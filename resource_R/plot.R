# Plot
par("mar") # check margin
par(mar=c(4,4,4,2)) #button, left, up, right margin
plot(X_data, Y_data, pch=".", main="Title", ylab="Y_axis", xlab="X_axis",xlim=c(0,10),ylim=c(0,2.0))
hist(data, breaks= 1000, main="title", ,xlim=c(0,10), ylim=c(0,2.0))
#save plot
  # 1. Open jpeg file
  jpeg(filename= paste0(substr(x, start=3, stop=15), "_summits_peakDistribution.jpg"), width = 800, height = 400, pointsize = 20, quality = 100,)
  # 2. Create the plot
  hist(peak_data$score, xlim = c(0, 400), breaks = 1000, main = substr(x, start=3, stop=15))
  # 3. Close the file
  dev.off()

###heatmap###
library(pheatmap)
rownames(sampleDistMatrix) <- vsd$condition
colnames(sampleDistMatrix) <- colnames(countData)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#annotation
color_code=data.frame(sample = rownames(sampleDistMatrix))
rownames(color_code) = colnames(sampleDistMatrix)
#heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         #  cutree_rows = 2,
         #  cutree_cols = 2,
         col=colors,
         annotation_col = color_code,
         show_colnames = F,
         show_rownames = F,
         main = "Sample-to-sample distances(HCC without hepatitis))")


###ggplot###
#line plot
df = data.frame(
    method=c(rep('ideal',4), rep('reduction',4), rep('atomic',4)), 
    nprc=rep(c(2, 4, 6, 8), 3), 
    speedup=c(c(2, 4, 6, 8), c(4.27, 2.17, 1.45, 1.17, 51.53, 53.73, 47.11, 47.74)^(-1)*8.49942))

ggplot(data=df, aes(x=nprc, y=speedup, group=method)) +
    geom_line(aes(color=method)) +
    geom_point(aes(shape=method, color=method), size=2) +
    
#scatter plot
ggplot(dataframe, aes(x=colname1, y=colname2, color = padj_sex < 0.05 & padj_normal < 0.05)) +
    geom_point(size=2, shape=16) +
    scale_color_manual(values=c('green','red')) #color selection per group
#boxplot
plot_df = data.frame(score=x, target=y)
  button10 =  dplyr::top_n(plot_df, -dim(plot_df)[1]*0.1,target) #lowest
  button10['group'] = 'button'
  top10 = dplyr::top_n(plot_df, dim(plot_df)[1]*0.1,target) #highest
  top10['group'] = 'top'
  plot_df = rbind(button10, top10)
  rm(button10, top10)

ggplot(boxplot_data.frame, aes(x = sex, y=value)) +
    geom_boxplot(width=0.8, outlier.size = NULL, outlier.shape=16, outlier.colour='black')  +
    geom_signif(comparisons = list(c("button","top")),
              map_signif_level = TRUE, step_increase = 0.1)
	 ## x : categorical variable ; 
	 ## y : numerical variable

#violin plot with box plot
ggplot(data.frame, aes(x = colname_x_in_data.frame, y=colname_y, fill=colname)) +
    geom_violin() +
    stat_summary(fun=mean/median, geom="point", fill="red", shape=21, size=2.5) +
    geom_boxplot(width=0.1, outlier.size = NULL, outlier.shape=16, outlier.colour='black', fill='white') +
    stat_summary(fun="mean", geom="point", shape=21, size=3, fill="blue") +
    ggtitle("plot title")

# Bar chart
ggplot(plot_data.frame, aes(x=col_name1, y= col_name2)) +
    geom_col(aes(fill=col_name3), color=color_vector1, size=1.1) + #size indicated the border line of each bar, color indicates the color of each line.

    scale_fill_gradientn(colours = color_vector2) +
    scale_y_discrete(limits = name_inSequence) + 
        #name_inSequence contains the element of y-axis in sequence
    ggtitle("title") +
    theme(panel.grid.major = element_line(colour="#f0f0f0"))



    geom_hline(yintercept=0, size=.1) +
    geom_vline(xintercept=0, size=.1)

#additional thing
sp + labs(x='xlab', y='ylab', color="color-level legend name")
sp + ggtitle("graph title")
sp + geom_hline(yintercept=0, size=.1) +
sp + geom_hline(yintercept=1, size=.1, linetype="dotted")
sp + geom_vline(xintercept=0, size=.1, color="red")

ggsave(paste0("./figure/",".png"), width=25, height = 15, units='cm', limitsize = FALSE)

### regression

#linear regression
fit = lm(y ~ x, data=input_data.frame)
    #y: target variable : the variable whose values are to be modeled and predicted by other variables
    #x: predictor variable : the name given to an independent variable used in regression analyses
summary(fit)

##to draw a plot
a = 
b = 

##preduct function - versoin 1
y_predicted = a * x + b
plot(predicted_y, y) # predicted Y vs actual Y 
abline(a=0,b=1, col="red") # drawing diagnoal line

##preduct function - versoin 2
predicted_y = predict(fit)
plot(predicted_y, y)
abline(a=0,b=1, col="red") #diagnoal line

##drawing regression lines over Y ~ X
plot (x, y)
abline(lm(y ~ x, data=input_data.frame), col="red")
summary(lm(y~x7, data=input_data.frame))

#multiple regression
fit = lm(y ~ x1 + x2 + x3 + x4, data=input_data.frame)
AIC(fit.multiple)
    # Akaike information criterion = AIC 
    # Lower AIC is better in that model is simpler and better fit for the data
plot(predict(fit.multiple, input_data.frame), y)
abline(a=0,b=1, col="red")
summary(fit.multiple)

#model performance evaulation
eval_results <- function(true, predicted, df) {
    SSE <- sum((predicted - true)^2)
    SST <- sum((true - mean(true))^2)
    R_square <- 1 - SSE / SST
    RMSE = sqrt(SSE/nrow(df))
    
    
    # Model performance metrics
    data.frame(
        RMSE = RMSE,
        Rsquare = R_square)
}

# Lasso regression, regression that gives a penalty to complex model, choose best features (L2 regularization)
if (!require(glmnet)) install.packages("glmnet")
require(glmnet)

#Histogram (two group)
labs <- paste("common_string_", c("group1", "group2"))
data.frame %>%
    mutate(label = ifelse(predicate_for_group1 == 1, labs[1], labs[2])) %>%
    ggplot(aes(x = colName_for_Xaxis_score)) +
    geom_histogram(color = "white") +
    facet_wrap(~colName_for_grouping) +
    xlab("Probability of X") +
    theme_bw()

    e.g.
        labs <- paste("Actual school type attended:", c("Catholic", "Public"))
        prs_df %>%
            mutate(catholic = ifelse(catholic == 1, labs[1], labs[2])) %>%
            ggplot(aes(x = pr_score)) +
            geom_histogram(color = "white") +
            facet_wrap(~catholic) +
            xlab("Probability of going to Catholic school") +
            theme_bw()












### Dimensional reduction ###
# eigen vector
library(ggplot2)
set.seed(1)
base1 <- rnorm(20)
x1 <- base1 + rnorm(20) * 0.4
y1 <- base1 + rnorm(20) * 0.4

dat <- data.frame(id = 1:20, x1, y1)

ggplot(dat, aes(x = x1, y = y1)) + geom_point(size = 5) + 
    geom_text(aes(label = id), col = "white") + 
    geom_abline(slope = 1, intercept = 0, linetype = 2) + 
    geom_abline(slope = -1, intercept = 0, linetype = 2) + 
    xlim(-3.5, 3.5) + ylim(-3.5, 3.5)

### rotation by multiplications with eigen vectors
e <- eigen(cov(dat[, -1]))
rotated <- data.frame(id = dat$id, as.matrix(dat[, 2:3]) %*% e$vectors)

ggplot(rotated, aes(x = X1, y = X2)) + 
    geom_point(size = 5) + 
    geom_text(aes(label = id), col = "white") +
    geom_abline(intercept = 0, slope = 0)

#dataset
library(datasets)
data(iris)

# Dimensional reduction ###
#pca
pca.out = prcomp(iris[,1:4])
plot(pca.out$x[,1:2], col=iris$Species, pch=16)


# multidimensional scaling (MDS)
iris.dist = dist((iris[,1:4]), method="euclidean")
mds.out = cmdscale(iris.dist, eig = T, k = 2)
plot(mds.out$points[,1], mds.out$points[,2], col=iris$Species, pch=16, main = "MDS plot")

# non-metric multidimensional scaling (NMDS)
require(vegan)
set.seed(2)

nmds.out = metaMDS(comm = iris[,1:4],k=2, trymax=100)
plot(nmds.out$points, col=iris$Species, xlab="NMDS1", ylab="NMDS2", pch=16)
stressplot(nmds.out)


# tSNE analysis
require(Rtsne)
set.seed(1)
tsne.out = Rtsne(iris[,1:4], check_duplicates = F)
plot(tsne.out$Y[,1:2], col=iris$Species, xlab = "tSNE1", ylab = "tSNE2", pch=16 )

# UMAP analysis
require(umap)
umap.out = umap(iris[,1:4])
plot(umap.out$layout[,1:2], col=iris$Species, xlab = "UMAP1", ylab = "UMAP2", pch=16)


###Vendia gram
# install.packages('VennDiagram')
library(VennDiagram)
grid.newpage()
draw.triple.venn(
    area1 = 1000,
    area2 = 500,
    area3 = 100,
    n12 = 30,
    n23 = 20,
    n13 = 15,
    n123 = 10,
    category = c('Biggest', 'Bigger', 'Big'),
    fill = RColorBrewer::brewer.pal(3, 'Accent'),
    lty = 'blank',
    cex = rep(1.2, 7),
    cat.cex = rep(1.5, 3)
)
