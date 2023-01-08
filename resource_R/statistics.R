# statistics
mean()
median()
sd()
## mode ()
get_mode = gene_score_anno$score_MinMaxNormalized
get_mode = unlist(lapply(get_mode, function(x) round(x, 3)))
get_mode = table(get_mode)
get_mode = as.numeric(names(get_mode)[which(get_mode==max(get_mode))])

## summary of dataframe (tidyverse)
base 참고

## summary of dataframe (data.table)
Ref: https://tariat.tistory.com/576
library(data.table)
summary_df = data.table(summary_df)
summary_df = summary_df[,.(mean=mean(colname_summary),median=median(colname_summary), sd=sd(colname_summary)),by = .(colname_group)]


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

#normalization (z-score)
## 1
distribution = correlation_coefficient
distribution = (distribution - mean(distribution)) / sd(distribution) #normalization
summary(distribution)
hist(distribution)

## 2
pathway_matrix = data.frame(lapply(pathway_matrix, as.numeric)) # convert list into dataframe
table(is.na(pathway_matrix)) #check NA
pathway_matrix = as.matrix(pathway_matrix) #convert dataframe into matrix
pathway_matrix_zscore = (pathway_matrix - mean(pathway_matrix))/sd(pathway_matrix) # transform into z-score

# Min-max normalization
distribution = (distribution - min(distribution)) / (max(distribution) - min(distribution)) #normalization                         
# Goodness of fit test : Kolmogorov–Smirnov test
# ref: https://www.r-bloggers.com/2015/01/goodness-of-fit-test-in-r/
ks_test = ks.test(y_1, y_2) # y_1 and y_2 are probability density

# Empirical Cumulative Distribution Function
Ref:https://stat.ethz.ch/R-manual/R-devel/library/stats/html/ecdf.html
                         
x = random_model_TotalFlux$total_flux
x_target = 31787.86967576691
Fn <- ecdf(x)
Fn()
plot(Fn)
lines(x_target, Fn(x_target), col = "red")
                         
# basic statistic by column in a dataframe
apply(data.frame, 1, median) #median by row
apply(data.frame, 2, median) #median by column

# Statistical Test
## non-parametric test
- Wilcoxon rank sum test
Ref (Youtube, concept): https://www.youtube.com/watch?v=la_M-OJvmI8&ab_channel=JaradNiemi  
Ref (R-code) : http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r  
```
wilcox.test(x, y, alternative = "two.sided")
> x,y: numeric vectors
> alternative: the alternative hypothesis. Allowed value is one of “two.sided” (default), “greater” or “less”.
```
