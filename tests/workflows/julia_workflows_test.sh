#!/usr/bin/env bash
# SPDX-License-Identifier: MPL-2.0
# Configuration regression tests for the Julia workflow security changes.
# Usage: bash tests/workflows/julia_workflows_test.sh
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
WF="$REPO/.github/workflows"

for workflow in julia-ci julia-docs; do
  # Match the security linter: actions-lock may prepend comments, but the
  # license must occur before the workflow body, not in a step or script.
  awk '
    /^# SPDX-License-Identifier: MPL-2.0[[:space:]]*$/ { found = 1 }
    /^[^#[:space:]]/ { exit }
    END { exit !found }
  ' "$WF/$workflow.yml" || { echo "FAIL: $workflow missing leading SPDX header" >&2; exit 1; }

  # Read permission blocks at the workflow/job indentation used in this repo.
  # Require exactly the intended grants; reject inline permission shorthands,
  # extra scopes, missing blocks, or accidental elevation of another job.
  awk -v workflow="$workflow" '
    function check_permissions() {
      if (!in_permissions) return
      expected = level == 0 ? "read" : "write"
      if (entries != 1 || grant != "contents: " expected ||
          (level != 0 && (workflow != "julia-docs" || job != "docs"))) bad = 1
      in_permissions = 0
    }
    /^[[:space:]]*(#|$)/ { next }
    {
      match($0, /[^ ]/); indent = RSTART - 1
      if (in_permissions && indent <= level) check_permissions()
      if (in_permissions) {
        line = $0; sub(/^[ ]+/, "", line); sub(/[ ]+#.*$/, "", line); sub(/[ ]+$/, "", line)
        grant = line; entries++
      }
      if (/^jobs:/) in_jobs = 1
      else if (indent == 0) in_jobs = 0
      if (in_jobs && indent == 2) { job = $1; sub(/:$/, "", job) }
      if ($0 ~ /^[ ]*permissions:/) {
        if ($0 !~ /^[ ]*permissions:[ ]*(#.*)?$/ || (indent != 0 && indent != 4)) bad = 1
        level = indent; entries = 0; grant = ""; in_permissions = 1
        if (indent == 0) top_blocks++; else job_blocks++
      }
    }
    END {
      check_permissions()
      if (top_blocks != 1 || job_blocks != (workflow == "julia-docs" ? 1 : 0)) bad = 1
      if (bad) print "FAIL: " workflow " must default to contents: read; only the docs job may write" > "/dev/stderr"
      exit bad
    }
  ' "$WF/$workflow.yml"
  echo "PASS: $workflow license and least-privilege permissions"
done

# The cache pin comment changed with the v3 lock: guard against restoring the
# stale v2 comment or documenting a different commit from the actual lock.
cache_commit=$(awk '
  /^    '\''julia-actions\/cache@v3'\'':/ { cache = 1; next }
  cache && /^        commit:/ { gsub(/\047|sha1-/, "", $2); print $2; exit }
  /^    [^ ]/ { cache = 0 }
' "$WF/actions.lock")
[ -n "$cache_commit" ]
grep -Eq "^# +julia-actions/cache@v3 += $cache_commit( |$)" "$WF/julia-ci.yml"
if grep -q 'julia-actions/cache@v2' "$WF/julia-ci.yml" "$WF/julia-docs.yml"; then
  echo 'FAIL: stale Julia cache v2 reference' >&2
  exit 1
fi
echo 'PASS: cache v3 documentation agrees with its locked commit'
