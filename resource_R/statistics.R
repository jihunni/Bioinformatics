# quantile
quantile(correlation_coefficient, probs=c(0.025, 1-0.025))


# normal distribution
num_of_samples = 100000
x <- seq(-6, 6, by = .01)
y <- dnorm(x, mean = 2.5, sd = 0.5)
#plot(x,y)
#y <- rgamma(num_of_samples, shape = 10, scale = 3)

qnorm(p=0.025, mean=0, sd=1, lower.tail = TRUE)
qnorm(p=0.025, mean=0, sd=1, lower.tail = FALSE)

#normalization
distribution = correlation_coefficient
distribution = (distribution - mean(distribution)) / sd(distribution) #normalization
summary(distribution)
hist(distribution)

# Goodness of fit test : Kolmogorov–Smirnov test
# ref: https://www.r-bloggers.com/2015/01/goodness-of-fit-test-in-r/
ks_test = ks.test(y_1, y_2) # y_1 and y_2 are probability density
