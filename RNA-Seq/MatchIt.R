library(tidyverse)
library("MatchIt")

### example - HCC application ###
reference: https://cran.r-project.org/web/packages/MatchIt/vignettes/MatchIt.html

#load data
save(HCC_geneCount, file='HCC_geneCount.rda') #log transformed data (log+1) is not applied
save(PanCan_HCC_clinic2, file='HCC_clinic_data.rda')

##LIHC clinical data from TCGA Assembler2
duplicated(clinic_LIHC) #No duplication
substr(colnames(HCC_geneCount),start=1,stop=12) %in% clinic_LIHC$bcr_patient_barcode #all samples are included in clinical data

HCC_clinic = data.frame(sample=colnames(HCC_geneCount), patient = substr(colnames(HCC_geneCount),start=1,stop=12), type = substr(colnames(HCC_geneCount),start=14,stop=15))
HCC_clinic$type[HCC_clinic$type=='01'] = 'tumor'
HCC_clinic$type[HCC_clinic$type=='02'] = 'tumor'
HCC_clinic$type[HCC_clinic$type=='11'] = 'normal'
table(HCC_clinic$type)

HCC_clinic2 = left_join(HCC_clinic, clinic_LIHC, by=c("patient" = "bcr_patient_barcode"))

save(HCC_clinic2, file='HCC_clinic2_2021.06.06.rda')



#data pre-processing
colnames(HCC_clinic2)
HCC_clinic3 = HCC_clinic2[,c(1,2,3,8,11,12,13,29,30,31,32,33,58, 19, 20, 21)]
str(HCC_clinic3)

#gender
HCC_clinic3$sample == HCC_clinic$sample
HCC_clinic3$gender = HCC_clinic2$gender

HCC_clinic3$gender[HCC_clinic3$gender == 'MALE'] = 'male'
HCC_clinic3$gender[HCC_clinic3$gender == 'FEMALE'] = 'female'
HCC_clinic3$gender = as.factor(HCC_clinic3$gender)
##Diagnosis Age
table(HCC_clinic3$age_at_diagnosis)
which(HCC_clinic3$age_at_diagnosis == '[Not Available]')
HCC_clinic3$age_at_diagnosis[which(HCC_clinic3$age_at_diagnosis == '[Not Available]')] = NA


HCC_clinic3$age_at_diagnosis = unlist(lapply(HCC_clinic3$age_at_diagnosis, as.numeric))
    #unlist:convert list into vector
HCC_clinic3$age_at_diagnosis[134] = mean(HCC_clinic3$age_at_diagnosis[-134]) #fill mean value in NA row



##tumor stage
##convert TX, T0, Na into 0
colnames(HCC_clinic3)[9] ='tumor_pathologic'
table(HCC_clinic3$tumor_pathologic)

HCC_clinic3$tumor_pathologic[which(is.na(HCC_clinic3$tumor_pathologic))] = 0
HCC_clinic3$tumor_pathologic = gsub("TX",0, HCC_clinic3$tumor_pathologic)
HCC_clinic3$tumor_pathologic = gsub("T1",1, HCC_clinic3$tumor_pathologic)
HCC_clinic3$tumor_pathologic = gsub("T2",2, HCC_clinic3$tumor_pathologic)
HCC_clinic3$tumor_pathologic = gsub("T2a",2, HCC_clinic3$tumor_pathologic)
HCC_clinic3$tumor_pathologic = gsub("T2b",2, HCC_clinic3$tumor_pathologic)
HCC_clinic3$tumor_pathologic = gsub("T3",3, HCC_clinic3$tumor_pathologic)
HCC_clinic3$tumor_pathologic = gsub("T3a",3, HCC_clinic3$tumor_pathologic)
HCC_clinic3$tumor_pathologic = gsub("T3b",3, HCC_clinic3$tumor_pathologic)
HCC_clinic3$tumor_pathologic = gsub("T4",4, HCC_clinic3$tumor_pathologic)

HCC_clinic3$tumor_pathologic = unlist(lapply(HCC_clinic3$tumor_pathologic, as.numeric))


###chagned name
colnames(HCC_clinic3)[12] ='tumor_stage'
#HCC_clinic3$tumor_stage = HCC_clinic2$ajcc_pathologic_tumor_stage ##reverse
table(HCC_clinic3$tumor_stage)
is.na(HCC_clinic3$tumor_stage)

#Be careful! before chaning string, check whether it gives an impact to another strings.
HCC_clinic3$tumor_stage = gsub("Stage I",1, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("Stage II",2, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("Stage III",3, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("Stage IIIA",3, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("Stage IIIB",3, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("Stage IIIC",3, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("Stage IV",4, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("Stage IVA",4, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("Stage IVB",4, HCC_clinic3$tumor_stage)

HCC_clinic3$tumor_stage = gsub("1",1, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("1I",2, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("1II",3, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("1IIA",3, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("1IIB",3, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("1IIC",3, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("2I", 3, HCC_clinic3$tumor_stage) ? 
HCC_clinic3$tumor_stage = gsub("2IA", 3, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("2IB", 3, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("2IC", 3, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("3A", 3, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("3B", 3, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("3C", 3, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("1V",4, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("1VA",4, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("1VB",4, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("4A",4, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage = gsub("4B",4, HCC_clinic3$tumor_stage)

##randomly distribute NA value into stage
#HCC_clinic3$tumor_stage = gsub("[Discrepancy]",1, HCC_clinic3$tumor_stage)
HCC_clinic3$tumor_stage[HCC_clinic3$tumor_stage == '[Discrepancy]'] = '[Not Available]'
HCC_clinic3$tumor_stage[which(HCC_clinic3$tumor_stage == "[Not Available]")][1:16] = 1
HCC_clinic3$tumor_stage[which(HCC_clinic3$tumor_stage == "[Not Available]")][1:9] = 2
HCC_clinic3$tumor_stage[which(HCC_clinic3$tumor_stage == "[Not Available]")][1:6] = 3

HCC_clinic3$tumor_stage = unlist(lapply(HCC_clinic3$tumor_stage, as.numeric))

## Race Category
##HCC_clinic3$race = HCC_clinic2$race
table(HCC_clinic3$ethnicity)
table(HCC_clinic3$race)
HCC_clinic3$race[HCC_clinic3$ethnicity == 'HISPANIC OR LATINO'] = 'HISPANIC OR LATINO'
HCC_clinic3$race[HCC_clinic3$race == 'AMERICAN INDIAN OR ALASKA NATIVE'] = 'Hispanic/Native'
HCC_clinic3$race[HCC_clinic3$race == 'HISPANIC OR LATINO'] = 'Hispanic/Native'
HCC_clinic3$race[HCC_clinic3$race == '[Not Available]'] = 'Not'
HCC_clinic3$race[HCC_clinic3$race == '[Not Evaluated]'] = 'Not'
HCC_clinic3$race[HCC_clinic3$race == '[Unknown]'] = 'Not'

HCC_clinic3$race = as.factor(HCC_clinic3$race)

#metastasis_pathologic
colnames(HCC_clinic3)[11] = 'metastasis_pathologic'
table(HCC_clinic3$metastasis_pathologic)
HCC_clinic3$metastasis_pathologic = as.factor(HCC_clinic3$metastasis_pathologic)

#HBV and HCV
library(stringr)


#initialization
table(HCC_clinic3$history_hepato_carcinoma_risk_factors)
HCC_clinic3['HBV'] = NA
HCC_clinic3['HCV'] = NA
HCC_clinic3['alcohol'] = NA

##history_hepato_carcinoma_risk_factors
index = str_detect(HCC_clinic3$history_hepato_carcinoma_risk_factors, 'Hepatitis B')
HCC_clinic3$HBV[index] = TRUE

index = str_detect(HCC_clinic3$history_hepato_carcinoma_risk_factors, 'Hepatitis C')
HCC_clinic3$HCV[index] = TRUE

index = str_detect(HCC_clinic3$history_hepato_carcinoma_risk_factors, 'Alcohol consumption')
HCC_clinic3$alcohol[index] = TRUE
HCC_clinic3$alcohol[is.na(HCC_clinic3$alcohol)] = FALSE


##iral_hepatitis_serology
table(HCC_clinic3$viral_hepatitis_serology)

index = str_detect(HCC_clinic3$viral_hepatitis_serology, 'Hepatitis  C Antibody')
HCC_clinic3$HCV[index] = TRUE
index = str_detect(HCC_clinic3$viral_hepatitis_serology, 'Hepatitis C Virus RNA')
HCC_clinic3$HCV[index] = TRUE
HCC_clinic3$HCV[is.na(HCC_clinic3$HCV)] = FALSE



index = str_detect(HCC_clinic3$viral_hepatitis_serology, 'Hepatitis B Surface Antigen')
HCC_clinic3$HBV[index] = TRUE

index = str_detect(HCC_clinic3$viral_hepatitis_serology, 'HBV Core Antibody')
HCC_clinic3$HBV[index] = FALSE
HCC_clinic3$HBV[is.na(HCC_clinic3$HBV)] = FALSE

###HBV positive
Hepatitis  C Antibody|HBV Surface Antibody 
Hepatitis  C Antibody|Hepatitis B Surface Antigen|HBV DNA 
Hepatitis  C Antibody|Hepatitis B Surface Antigen|HBV Surface Antibody 
Hepatitis  C Antibody|Hepatitis B Surface Antigen|HBV Surface Antibody|HBV DNA 
Hepatitis  C Antibody|Hepatitis C Virus RNA|HCV Genotype|Hepatitis B Surface Antigen 
Hepatitis  C Antibody|Hepatitis C Virus RNA|Hepatitis B Surface Antigen|HBV Surface Antibody 
Hepatitis B Surface Antigen 
Hepatitis B Surface Antigen|HBV DNA 
Hepatitis B Surface Antigen|HBV Surface Antibody 
Hepatitis B Surface Antigen|HBV Surface Antibody|HBV DNA
Hepatitis C Virus RNA|HBV Surface Antibody|HBV DNA 
Hepatitis C Virus RNA|Hepatitis B Surface Antigen

### HBV negative
HBV Surface Antibody 
HCV Genotype 
Hepatitis  C Antibody 
Hepatitis  C Antibody|HBV Core Antibody 
Hepatitis  C Antibody|HBV Surface Antibody 
Hepatitis  C Antibody|Hepatitis B Surface Antigen 
Hepatitis  C Antibody|Hepatitis B Surface Antigen|HBV Core Antibody 
Hepatitis  C Antibody|Hepatitis B Surface Antigen|HBV Core Antibody|HBV DNA 
Hepatitis  C Antibody|Hepatitis B Surface Antigen|HBV DNA 
Hepatitis  C Antibody|Hepatitis B Surface Antigen|HBV Surface Antibody 
Hepatitis  C Antibody|Hepatitis B Surface Antigen|HBV Surface Antibody|HBV Core Antibody 
Hepatitis  C Antibody|Hepatitis B Surface Antigen|HBV Surface Antibody|HBV DNA
Hepatitis  C Antibody|Hepatitis C Virus RNA|HCV Genotype 
Hepatitis  C Antibody|Hepatitis C Virus RNA|HCV Genotype|Hepatitis B Surface Antigen
Hepatitis  C Antibody|Hepatitis C Virus RNA|Hepatitis B Surface Antigen|HBV Surface Antibody 
Hepatitis B Surface Antigen 
Hepatitis B Surface Antigen|HBV Core Antibody 
Hepatitis B Surface Antigen|HBV Core Antibody|HBV DNA 0
Hepatitis B Surface Antigen|HBV DNA 
Hepatitis B Surface Antigen|HBV Surface Antibody 0
Hepatitis B Surface Antigen|HBV Surface Antibody|HBV Core Antibody 
Hepatitis B Surface Antigen|HBV Surface Antibody|HBV Core Antibody|HBV DNA 
Hepatitis B Surface Antigen|HBV Surface Antibody|HBV DNA 
Hepatitis C Virus RNA 0
Hepatitis C Virus RNA|HBV Surface Antibody|HBV DNA 
Hepatitis C Virus RNA|HCV Genotype 
Hepatitis C Virus RNA|HCV Genotype|Hepatitis B Surface Antigen|HBV Surface Antibody|HBV Core Antibody 0
Hepatitis C Virus RNA|Hepatitis B Surface Antigen 
Hepatitis C Virus RNA|Hepatitis B Surface Antigen|HBV Core Antibody 
Hepatitis C Virus RNA|Hepatitis B Surface Antigen|HBV Core Antibody|HBV DNA
Hepatitis C Virus RNA|Hepatitis B Surface Antigen|HBV Surface Antibody|HBV Core Antibody 

table(HCC_clinic3$HBV)
TRUE
175

table(HCC_clinic3$HCV)
TRUE 
144 


#nodes_pathologic
colnames(HCC_clinic3)[10] = 'nodes_pathologic'
table(HCC_clinic3$nodes_pathologic)
HCC_clinic3$nodes_pathologic = gsub('N0', 0,HCC_clinic3$nodes_pathologic)
HCC_clinic3$nodes_pathologic = gsub('N1', 1,HCC_clinic3$nodes_pathologic)
HCC_clinic3$nodes_pathologic = gsub('[Not Available]','NX',HCC_clinic3$nodes_pathologic)
HCC_clinic3$nodes_pathologic[HCC_clinic3$nodes_pathologic =='[NXNXNXNXNXNXNXNXNXNXNXNXNX]'] = 'NX'

HCC_clinic3$nodes_pathologic = as.factor(HCC_clinic3$nodes_pathologic)

#wrap up the data parsing
save(HCC_clinic3, file='HCC_clinic3_2021.06.06.rda')

#normalization
norm_minmax <- function(x){(x- min(x)) /(max(x)-min(x))}
HCC_clinic4 = HCC_clinic3
HCC_clinic4$tumor_pathologic = norm_minmax(HCC_clinic3$tumor_pathologic)
HCC_clinic4$tumor_stage = norm_minmax(HCC_clinic3$tumor_stage)
HCC_clinic4$age_at_diagnosis = norm_minmax(HCC_clinic3$age_at_diagnosis)


save(HCC_clinic4, file='HCC_clinic4_2021.06.10.rda')

table(substr(HCC_clinic4$sample,start=14,stop=15))
HCC_clinic_tumor = HCC_clinic4[substr(HCC_clinic4$sample,start=14,stop=15) == '01' | substr(HCC_clinic4$sample,start=14,stop=15) == '02',]
HCC_clinic_normal = HCC_clinic4[substr(HCC_clinic4$sample,start=14,stop=15) == '11',]

table(HCC_clinic_tumor$HBV == FALSE & HCC_clinic_tumor$HCV == FALSE)
table(HCC_clinic_tumor$HBV == TRUE & HCC_clinic_tumor$HCV == FALSE)
table(HCC_clinic_tumor$HBV == FALSE & HCC_clinic_tumor$HCV == TRUE)
table(HCC_clinic_tumor$HBV == TRUE & HCC_clinic_tumor$HCV == TRUE)

table(HCC_clinic_tumor[HCC_clinic_tumor$HBV == FALSE & HCC_clinic_tumor$HCV == FALSE,]$gender == 'MALE')
table(HCC_clinic_tumor[HCC_clinic_tumor$HBV == TRUE & HCC_clinic_tumor$HCV == FALSE,]$gender == 'MALE')
table(HCC_clinic_tumor[HCC_clinic_tumor$HBV == FALSE & HCC_clinic_tumor$HCV == TRUE,]$gender == 'MALE')
table(HCC_clinic_tumor[HCC_clinic_tumor$HBV == TRUE & HCC_clinic_tumor$HCV == TRUE,]$gender == 'MALE')

HCC_noH = HCC_clinic_tumor[HCC_clinic_tumor$HBV == FALSE & HCC_clinic_tumor$HCV == FALSE,] 
HCC_HBV = HCC_clinic_tumor[HCC_clinic_tumor$HBV == TRUE & HCC_clinic_tumor$HCV == FALSE,]
HCC_HCV = HCC_clinic_tumor[HCC_clinic_tumor$HBV == FALSE & HCC_clinic_tumor$HCV == TRUE,]
HCC_both = HCC_clinic_tumor[HCC_clinic_tumor$HBV == TRUE & HCC_clinic_tumor$HCV == TRUE,]


### propensity score algorithm ###
1. load library and check data
library("MatchIt")
colnames(HCC_noH)

2. Check initial imbalance
# No matching; constructing a pre-match matchit object
m.noH_non <- matchit(gender ~ race+tumor_pathologic+nodes_pathologic+metastasis_pathologic+tumor_stage+age_at_diagnosis+HBV+HCV+alcohol, data = HCC_noH, method = NULL, distance = "glm")

# Checking balance prior to matching
summary(m.out0)

3. Matching with weight
##start
m.noH <- matchit(gender ~ race+tumor_pathologic+nodes_pathologic+metastasis_pathologic+tumor_stage+age_at_diagnosis+HBV+HCV+alcohol, data = HCC_noH, method = "full", distance = "glm", link = "probit", caliper = 0.2, estimand = "ATE")

##weight
m.noH <- matchit(gender ~ race+tumor_pathologic+nodes_pathologic+metastasis_pathologic+tumor_stage+age_at_diagnosis+HBV+HCV+alcohol, data = HCC_noH, method = "full", distance = "glm", link = "probit", caliper = 0.2, estimand = "ATE", s.weights = m.noH$weights)
m.noH
plot(summary(m.noH))

save(m.noH, file='m.noH_2021.06.10.rda')

4. Assessing the quality of matches
# Checking balance after NN matching
summary(m.noH)
    ##check std.mean diff. < 0.1

#visualize the distribution of propensity scores - jitter plot
plot(m.noH, type = "jitter", interactive = FALSE)

#visualize the distribution of propensity scores - jitter plot - QQ plot
plot(m.out1, type = "qq", interactive = FALSE,
     which.xs = c("neoplasm_stage", "tumor_stage", 'race'))
plot(summary(m.out1))


#visualization
reference:https://cran.r-project.org/web/packages/MatchIt/vignettes/assessing-balance.html#plot.matchit
library(cobalt)

unadjust.match =  matchit(sex ~ age + race + tumor_stage + metastatis_stage + neoplasm_stage, data = tumorFree_noHepatitis, method = NULL , distance = "glm", link = "probit", caliper = 0.2, estimand = "ATE")


##visualization option
match_for_visualization = m.noH_non
match_for_visualization = m.noH
#which = "unadjusted"
which = 'adjusted'

##visualizae
love.plot(match_for_visualization, binary = "std")
bal.plot(match_for_visualization, var.name = "age", which = which)
bal.plot(match_for_visualization, var.name = "neoplasm_stage", which = which)
bal.plot(match_for_visualization, var.name = "metastatis_stage", which = which)
bal.plot(match_for_visualization, var.name = "tumor_stage", which = which)
bal.plot(match_for_visualization, var.name = "race", which = which)
bal.plot(match_for_visualization, var.name = "distance", which = which,
         type = "histogram", mirror = TRUE)

m.noH.data <- match.data(m.noH)
m.noH.data = m.noH.data[m.noH.data$weights != 0,]
save(m.noH.data, file='m.noH.data_2021.06.10.rda')

# Checking balance after full matching
summary(m.out2, un = FALSE)

plot(m.out2, type = "qq", interactive = FALSE, which.xs = c("age","metastatis_stage", "neoplasm_stage", "tumor_stage", 'race'))


6. Estimating the Treatment effects
m.data2 <- match.data(m.out2)
head(m.data2)

library("lmtest") #coeftest
library("sandwich") #vcovCL

fit1 <- lm(sex ~ age   +status + race + metastatis_stage, data = HCC_clinic3, weights = m.data2$weights)
coeftest(fit1, vcov. = vcovCL, cluster = ~subclass)

m.data2 <- match.data(m.out2)
fit2 <- lm(re78 ~ treat + age + educ + race + married + nodegree + re74 + re75, data = m.data2, weights = weights)
coeftest(fit2, vcov. = vcovCL, cluster = ~subclass)

