# TEMPORARY probe script — deleted before the PR lands.
#
# Answers, with a real pscl fit, the questions the package's design rests on:
#   1. what the zero-part coefficients mean in each model (measured, not assumed);
#   2. what pscl does in the degenerate cases the package refuses;
#   3. what pscl does at the size-parameter bounds.
# Nothing here is copied into the package: this is measurement.
options(width = 200)
set.seed(11)

section <- function(name) cat("\n== ", name, " ==\n", sep = "")
attempt <- function(tag, expr) {
    r <- tryCatch(expr, error = function(e) e, warning = function(w) w)
    if (inherits(r, "error")) {
        cat(tag, "ERROR:", conditionMessage(r), "\n")
    } else if (inherits(r, "warning")) {
        cat(tag, "WARNING:", conditionMessage(r), "\n")
    } else {
        cat(tag, "ok\n")
    }
    invisible(r)
}

n <- 60
g <- factor(rep(c("A", "B"), each = n / 2), levels = c("A", "B"))
s <- exp(rnorm(n, 0, 0.4))
p_true <- plogis(0.4 + 0.5 * (g == "B") + 0.8 * log(s))
mu_true <- exp(1.2 + 0.9 * (g == "B")) * s
theta_true <- 2.5
z <- rbinom(n, 1, p_true)
y <- ifelse(z == 1, rnbinom(n, mu = mu_true, size = theta_true), 0)
d <- data.frame(y = y, g = g, s = s)

# --- 1. hurdle -----------------------------------------------------------------
section("hurdle")
fit <- pscl::hurdle(y ~ g + offset(log(s)) | g + log(s), data = d,
                    dist = "negbin", zero.dist = "binomial")
cat("count coef:\n"); print(coef(fit, model = "count"))
cat("zero coef:\n"); print(coef(fit, model = "zero"))
cat("theta:", fit$theta, "\n")
cat("logLik:", as.numeric(logLik(fit)), " df:", attr(logLik(fit), "df"), "\n")
cat("converged:", fit$converged, "\n")
eta0 <- model.matrix(~ g + log(s), data = d) %*% coef(fit, model = "zero")
ph <- plogis(eta0)
cat("cor(plogis(eta0), y==0):", cor(ph, as.numeric(y == 0)),
    " cor(plogis(eta0), y>0):", cor(ph, as.numeric(y > 0)), "\n")
cat("mean(y==0):", mean(y == 0), " mean(plogis(eta0)):", mean(ph),
    " mean(y>0):", mean(y > 0), "\n")
cat("sd(plogis(eta0)):", sd(ph), "\n")
cat("zero-part SEs:\n"); print(summary(fit)$coefficients$zero[, 2, drop = FALSE])
fit0 <- pscl::hurdle(y ~ offset(log(s)) | log(s), data = d,
                     dist = "negbin", zero.dist = "binomial")
cat("null logLik:", as.numeric(logLik(fit0)), " df:", attr(logLik(fit0), "df"), "\n")
lr <- 2 * (as.numeric(logLik(fit)) - as.numeric(logLik(fit0)))
cat("LRT:", lr, " p(chi2_2):", pchisq(lr, df = 2, lower.tail = FALSE), "\n")

# --- 2. zeroinfl ---------------------------------------------------------------
section("zeroinfl")
fitz <- pscl::zeroinfl(y ~ g + offset(log(s)) | g + log(s), data = d,
                       dist = "negbin")
cat("count coef:\n"); print(coef(fitz, model = "count"))
cat("zero coef:\n"); print(coef(fitz, model = "zero"))
cat("theta:", fitz$theta, "\n")
cat("logLik:", as.numeric(logLik(fitz)), " df:", attr(logLik(fitz), "df"), "\n")
cat("converged:", fitz$converged, "\n")
etaz <- model.matrix(~ g + log(s), data = d) %*% coef(fitz, model = "zero")
pz <- plogis(etaz)
cat("cor(plogis(eta0), y==0):", cor(pz, as.numeric(y == 0)),
    " cor(plogis(eta0), y>0):", cor(pz, as.numeric(y > 0)), "\n")
cat("mean(y==0):", mean(y == 0), " mean(plogis(eta0)):", mean(pz), "\n")
cat("zero-part SEs:\n"); print(summary(fitz)$coefficients$zero[, 2, drop = FALSE])
fitz0 <- pscl::zeroinfl(y ~ offset(log(s)) | log(s), data = d, dist = "negbin")
cat("null logLik:", as.numeric(logLik(fitz0)), " df:", attr(logLik(fitz0), "df"), "\n")
lrz <- 2 * (as.numeric(logLik(fitz)) - as.numeric(logLik(fitz0)))
cat("LRT:", lrz, " p(chi2_2):", pchisq(lrz, df = 2, lower.tail = FALSE), "\n")

# --- 3. degenerate cases the package refuses -----------------------------------
section("degenerate cases")

d_sep <- d
d_sep$y[d_sep$g == "A"] <- 0
attempt("(sep) all-zero group A, hurdle:", pscl::hurdle(y ~ g + offset(log(s)) | g + log(s),
    data = d_sep, dist = "negbin", zero.dist = "binomial"))
attempt("(sep) all-zero group A, zinb:", pscl::zeroinfl(y ~ g + offset(log(s)) | g + log(s),
    data = d_sep, dist = "negbin"))

d_sep2 <- d
d_sep2$y[d_sep2$g == "A"] <- 1
attempt("(sep) no zeros in group A, hurdle:", pscl::hurdle(y ~ g + offset(log(s)) | g + log(s),
    data = d_sep2, dist = "negbin", zero.dist = "binomial"))
attempt("(sep) no zeros in group A, zinb:", pscl::zeroinfl(y ~ g + offset(log(s)) | g + log(s),
    data = d_sep2, dist = "negbin"))

d_nz <- d
d_nz$y <- d_nz$y + 1
attempt("(sep) no zeros at all, hurdle:", pscl::hurdle(y ~ g + offset(log(s)) | g + log(s),
    data = d_nz, dist = "negbin", zero.dist = "binomial"))
attempt("(sep) no zeros at all, zinb:", pscl::zeroinfl(y ~ g + offset(log(s)) | g + log(s),
    data = d_nz, dist = "negbin"))

# --- 4. the small-sample floors -------------------------------------------------
section("small-sample floors")
d6 <- d[1:6, ]
d8 <- d[1:8, ]
r6 <- attempt("(tiny) 6 samples, zinb:", pscl::zeroinfl(y ~ g + offset(log(s)) | g + log(s),
    data = d6, dist = "negbin"))
if (!inherits(r6, "error") && !inherits(r6, "warning")) {
    cat("  theta:", r6$theta, " logLik:", as.numeric(logLik(r6)), "\n")
    print(summary(r6)$coefficients$zero)
}
r8 <- attempt("(tiny) 8 samples, zinb:", pscl::zeroinfl(y ~ g + offset(log(s)) | g + log(s),
    data = d8, dist = "negbin"))
if (!inherits(r8, "error") && !inherits(r8, "warning")) {
    cat("  theta:", r8$theta, " logLik:", as.numeric(logLik(r8)), "\n")
    print(summary(r8)$coefficients$zero)
}
d3 <- d[1:3, ]
attempt("(tiny) 3 samples, hurdle:", pscl::hurdle(y ~ g + offset(log(s)) | g + log(s),
    data = d3, dist = "negbin", zero.dist = "binomial"))
d4 <- d[1:4, ]
attempt("(tiny) 4 samples, hurdle:", pscl::hurdle(y ~ g + offset(log(s)) | g + log(s),
    data = d4, dist = "negbin", zero.dist = "binomial"))

# --- 5. the size parameter at its bounds ----------------------------------------
section("theta at bound")
set.seed(3)
y_pois <- rpois(n, 3)
d_pois <- data.frame(y = y_pois, g = g, s = s)
rb <- attempt("(bound) poisson-ish data, zinb:", pscl::zeroinfl(y ~ g + offset(log(s)) | g + log(s),
    data = d_pois, dist = "negbin"))
if (!inherits(rb, "error") && !inherits(rb, "warning")) {
    cat("  theta:", rb$theta, " logLik:", as.numeric(logLik(rb)), "\n")
}
rh <- attempt("(bound) poisson-ish data, hurdle:", pscl::hurdle(y ~ g + offset(log(s)) | g + log(s),
    data = d_pois, dist = "negbin", zero.dist = "binomial"))
if (!inherits(rh, "error") && !inherits(rh, "warning")) {
    cat("  theta:", rh$theta, " logLik:", as.numeric(logLik(rh)), "\n")
}

section("versions")
cat("R.version:", R.version.string, "\n")
cat("pscl:", as.character(packageVersion("pscl")),
    " MASS:", as.character(packageVersion("MASS")), "\n")
