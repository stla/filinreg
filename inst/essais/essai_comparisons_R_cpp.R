library(filinreg)

set.seed(666)
x <- rnorm(15)
y <- x + rnorm(15)

fi <- filinreg(y ~ x, L = 20L, distr = "normal", lucky = TRUE)
fiSummary(fi)
fiR <- filinregR(y ~ x, L = 20L, distr = "student", df = Inf, lucky = TRUE)
fiSummary(fiR) # idem


fi <- filinreg(y ~ x, L = 20L, distr = "student", df = 3, lucky = TRUE)
fiSummary(fi)
fiR <- filinregR(y ~ x, L = 20L, distr = "student", df = 3, lucky = TRUE)
fiSummary(fiR) # idem


fi <- filinreg(y ~ x, L = 20L, distr = "cauchy", lucky = TRUE)
fiSummary(fi)
fiR <- filinregR(y ~ x, L = 20L, distr = "student", df = 1, lucky = TRUE)
fiSummary(fiR) # idem


fi <- filinreg(y ~ x, L = 20L, distr = "logistic", lucky = TRUE)
fiSummary(fi)
fiR <- filinregR(y ~ x, L = 20L, distr = "logistic", lucky = TRUE)
fiSummary(fiR) # idem
