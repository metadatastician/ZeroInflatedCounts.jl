#!/usr/bin/env bash
# SPDX-License-Identifier: MPL-2.0
# Regression coverage for the rebuilt actions.lock. Offline; requires GNU awk.
# Usage: bash tests/workflows/actions_lock_test.sh
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CHECK="$REPO/scripts/check-lock-sync.sh"
FIXTURE=$(mktemp -d)
trap 'rm -rf "$FIXTURE"' EXIT
PASS=0

reset_fixture() {
  rm -rf "$FIXTURE/workflows"
  cp -R "$REPO/.github/workflows" "$FIXTURE/workflows"
}

expect_failure() {
  local label="$1" diagnostic="$2" rc=0
  bash "$CHECK" "$FIXTURE/workflows" > "$FIXTURE/result" 2>&1 || rc=$?
  if [ "$rc" -ne 1 ] || ! grep -Fq -- "$diagnostic" "$FIXTURE/result"; then
    echo "FAIL: $label (expected exit 1 and '$diagnostic', got $rc)" >&2
    cat "$FIXTURE/result" >&2
    exit 1
  fi
  echo "PASS: $label"
  PASS=$((PASS + 1))
}

# Delete a complete record, retaining the generated lockfile's formatting so
# the production checker's parser is exercised exactly as it is in CI.
remove_record() {
  local key="$1"
  awk -v key="    '$key':" '
    /^    [^ ]/ { omit = index($0, key) == 1 }
    /^[^ #]/ { omit = 0 }
    !omit { print }
  ' "$FIXTURE/workflows/actions.lock" > "$FIXTURE/lock"
  mv "$FIXTURE/lock" "$FIXTURE/workflows/actions.lock"
}

# Positive control uses every real workflow, including job-level reusable
# workflows, subpath actions, mixed-case owner names and zero-uses workflows.
bash "$CHECK" "$REPO/.github/workflows"
echo 'PASS: repository lock covers all workflows and transitive dependencies'
PASS=$((PASS + 1))

reset_fixture
remove_record '.github/workflows/build-notification.yml'
expect_failure 'previously unlocked workflow cannot disappear' 'not onboarded'

reset_fixture
remove_record '.github/workflows/lock-sync-gate.yml'
expect_failure 'workflow with no uses still needs an empty lock entry' 'UNLISTED WORKFLOWS'

reset_fixture
sed -i 's|julia-actions/cache@v3|julia-actions/cache@v2|g' "$FIXTURE/workflows/actions.lock"
expect_failure 'stale cache v2 lock is rejected after the v3 bump' 'refs missing from the lockfile: julia-actions/cache@v3'

reset_fixture
sed -i '/^dependencies:/,$d' "$FIXTURE/workflows/actions.lock"
expect_failure 'workflow lists alone cannot replace dependency records' 'DANGLING EDGES'

# These records were all introduced by this PR. Losing any one must fail,
# even if every workflow still lists the right uses reference.
for ref in \
  'julia-actions/cache@v3' \
  'julia-actions/setup-julia@v3.0.2' \
  'julia-actions/julia-buildpkg@v1' \
  'julia-actions/julia-docdeploy@v1' \
  'haskell-actions/setup@v2.12.1' \
  'sonarsource/sonarqube-scan-action@v8.3.0'; do
  reset_fixture
  remove_record "$ref"
  expect_failure "missing dependency: $ref" 'DANGLING EDGES'
done

reset_fixture
remove_record 'actions/attest@508db95dd578ae2727ebd6217d5ba78e4fbda05d'
expect_failure 'transitive-only dependency must also resolve' 'DANGLING EDGES'

reset_fixture
remove_record '.github/workflows/governance.yml'
expect_failure 'job-level reusable workflow must remain locked' 'not onboarded'

reset_fixture
rm "$FIXTURE/workflows/julia-docs.yml"
expect_failure 'deleted workflow cannot leave an orphan lock entry' 'workflow file that does not exist'

# Leaf dependencies may omit uses entirely (all new Julia records do).
# Validate immutable pins and GitHub identity metadata, which the sync checker
# does not inspect. This intentionally reads the generator's record layout.
awk '
  function check_record() {
    if (key == "") return
    if (commit !~ /^sha1-[0-9a-f]{40}$/ || commit == "sha1-0000000000000000000000000000000000000000" ||
        ref == "" || owner !~ /^[1-9][0-9]*$/ || repo !~ /^[1-9][0-9]*$/) {
      print "FAIL: incomplete dependency metadata: " key > "/dev/stderr"
      bad = 1
    }
    split(key, parts, "@")
    if (parts[2] ~ /^[0-9a-f]{40}$/ && commit != "sha1-" parts[2]) {
      print "FAIL: SHA reference disagrees with commit: " key > "/dev/stderr"
      bad = 1
    }
    count++
  }
  /^dependencies:/ { dependencies = 1; next }
  dependencies && /^    [^ ]/ {
    check_record()
    key = $1; gsub(/\047|:$/, "", key)
    ref = commit = owner = repo = ""
  }
  dependencies && /^        ref:/ { ref = $2; gsub(/\047/, "", ref) }
  dependencies && /^        commit:/ { commit = $2; gsub(/\047/, "", commit) }
  dependencies && /^        owner_id:/ { owner = $2 }
  dependencies && /^        repo_id:/ { repo = $2 }
  END {
    check_record()
    if (!count) { print "FAIL: no dependency records" > "/dev/stderr"; bad = 1 }
    if (!bad) print "PASS: immutable pins and identity metadata for " count " dependencies"
    exit bad
  }
' "$REPO/.github/workflows/actions.lock"
PASS=$((PASS + 1))
echo "actions lock: $PASS passed"
