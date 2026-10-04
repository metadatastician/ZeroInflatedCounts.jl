# TEMPORARY probe script — deleted before the PR lands.
set.seed(11)
n <- 40
g <- factor(rep(c("A","B"), each = n/2))
s <- exp(rnorm(n, 0, 0.4))
logitp <- 0.4 + 0.5 * (g == "B") + 0.8 * log(s)
p <- plogis(logitp)
mu <- exp(1.2 + 0.9 * (g == "B")) * s
theta <- 2.5
z <- rbinom(n, 1, p)
y <- ifelse(z == 1, rnbinom(n, mu = mu, size = theta), 0)
d <- data.frame(y = y, g = g, s = s)

fit <- pscl::hurdle(y ~ g + offset(log(s)) | g + log(s), data = d,
                    dist = "negbin", zero.dist = "binomial")
cat("== hurdle ==\n")
print(coef(fit, model = "count"))
print(coef(fit, model = "zero"))
cat("theta:", fit$theta, "\n")
cat("logLik:", as.numeric(logLik(fit)), " df:", attr(logLik(fit), "df"), "\n")
eta0 <- model.matrix(~ g + log(s), data = d) %*% coef(fit, model = "zero")
cat("cor(plogis(eta0), y==0):", cor(plogis(eta0), as.numeric(d$y == 0)), "\n")
cat("cor(plogis(eta0), y>0):", cor(plogis(eta0), as.numeric(d$y > 0)), "\n")
cat("mean(y==0):", mean(d$y == 0), " mean(plogis(eta0)):", mean(plogis(eta0)), "\n")
cat("names(fit):", names(fit), "\n")
print(summary(fit)$coefficients$count)
print(summary(fit)$coefficients$zero)

fit0 <- pscl::hurdle(y ~ offset(log(s)) | log(s), data = d,
                     dist = "negbin", zero.dist = "binomial")
cat("null logLik:", as.numeric(logLik(fit0)), "\n")
lrt <- 2 * (as.numeric(logLik(fit)) - as.numeric(logLik(fit0)))
cat("LRT:", lrt, " p:", pchisq(lrt, df = 2, lower.tail = FALSE), "\n")
cat("lrtest:\n"); print(pscl::lrtest(fit, fit0))

fitz <- pscl::zeroinfl(y ~ g + offset(log(s)) | g + log(s), data = d,
                       dist = "negbin")
cat("== zeroinfl ==\n")
print(coef(fitz, model = "count"))
print(coef(fitz, model = "zero"))
cat("theta:", fitz$theta, "\n")
cat("logLik:", as.numeric(logLik(fitz)), " df:", attr(logLik(fitz), "df"), "\n")
cat("converged:", fitz$converged, "\n")
print(summary(fitz)$coefficients$count)
print(summary(fitz)$coefficients$zero)
fitz0 <- pscl::zeroinfl(y ~ offset(log(s)) | log(s), data = d, dist = "negbin")
lrtz <- 2 * (as.numeric(logLik(fitz)) - as.numeric(logLik(fitz0)))
cat("null logLik:", as.numeric(logLik(fitz0)), " LRT:", lrtz,
    " p:", pchisq(lrtz, df = 2, lower.tail = FALSE), "\n")
cat("names(fitz):", names(fitz), "\n")

cat("== degenerate cases ==\n")
d2 <- d; d2$y[d2$g == "A"] <- 0
r2 <- try(pscl::hurdle(y ~ g + offset(log(s)) | g + log(s), data = d2, dist = "negbin", zero.dist = "binomial"), silent = TRUE)
cat("(a) all-zero group A, hurdle:", class(r2), "\n")
r2z <- try(pscl::zeroinfl(y ~ g + offset(log(s)) | g + log(s), data = d2, dist = "negbin"), silent = TRUE)
cat("(a) all-zero group A, zinb:", class(r2z), if (!inherits(r2z, "try-error")) as.numeric(logLik(r2z)) else "", "\n")
d3 <- d; d3$y <- d3$y + 1
r3 <- try(pscl::hurdle(y ~ g + offset(log(s)) | g + log(s), data = d3, dist = "negbin", zero.dist = "binomial"), silent = TRUE)
cat("(b) no zeros at all, hurdle:", class(r3), if (!inherits(r3, "try-error")) paste(coef(r3, "zero"), collapse = ",") else "", "\n")
r3z <- try(pscl::zeroinfl(y ~ g + offset(log(s)) | g + log(s), data = d3, dist = "negbin"), silent = TRUE)
cat("(b) no zeros at all, zinb:", class(r3z), "\n")

cat("== tiny sample floor ==\n")
dt <- d[1:6, ]
rt <- try(pscl::zeroinfl(y ~ g + offset(log(s)) | g + log(s), data = dt, dist = "negbin"), silent = TRUE)
cat("(c) 6 samples:", class(rt), if (inherits(rt, "try-error")) as.character(rt) else "", "\n")
d6 <- d[1:8, ]
rt6 <- try(pscl::zeroinfl(y ~ g + offset(log(s)) | log(s), data = d6, dist = "negbin"), silent = TRUE)
cat("(c) 8 samples:", class(rt6), if (inherits(rt6, "try-error")) as.character(rt6) else as.numeric(logLik(rt6)), "\n")

cat("== theta at bound ==\n")
set.seed(3)
yb <- rpois(n, 3)
db <- data.frame(y = yb, g = g, s = s)
rb <- try(pscl::zeroinfl(y ~ g + offset(log(s)) | g + log(s), data = db, dist = "negbin"), silent = TRUE)
cat("poisson-ish:", class(rb), if (!inherits(rb, "try-error")) paste("theta", rb$theta) else "", "\n")
cat("R.version:", R.version.string, "\n")
cat("pscl:", as.character(packageVersion("pscl")), "MASS:", as.character(packageVersion("MASS")), "\n")
