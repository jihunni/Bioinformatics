library(tidyverse)
library("MatchIt")

### example - HCC application ###
reference: https://cran.r-project.org/web/packages/MatchIt/vignettes/MatchIt.html

save(HCC_geneCount, file='HCC_geneCount.rda') #log transformed data (log+1)?
save(PanCan_HCC_clinic2, file='HCC_clinic_data.rda')

#drop_na on column sex
clinic=drop_na(PanCan_HCC_clinic2, Sex) 
rm(PanCan_HCC_clinic2)

HCC_genecount_filter = HCC_geneCount[,colnames(HCC_geneCount) %in% clinic$case]
rownames(HCC_genecount_filter) = rownames(HCC_geneCount)

#check the order of data
colnames(HCC_genecount_filter) == clinic$case

save(HCC_genecount_filter, clinic, file='HCC_geneCount_and_clinic.rda')

#data pre-processing
clinic_process = clinic[,c(1,2,6,7,8,23,24,25,27)]
str(clinic_process)
clinic_process=cbind(clinic_process, clinic["American Joint Committee on Cancer Metastasis Stage Code"])

##Diagnosis Age
table(clinic_process$`Diagnosis Age`)
which(is.na(clinic_process$`Diagnosis Age`)) #133
clinic_process$`Diagnosis Age` = unlist(lapply(clinic_process$`Diagnosis Age`, as.numeric))
    #unlist:convert list into vector
clinic_process$`Diagnosis Age`[133] = mean(clinic_process$`Diagnosis Age`[-133]) #fill mean value in NA row
##sex
clinic_process$Sex = factor(clinic_process$Sex)
is.na(clinic_process$Sex)
##Neoplasm Disease Stage American Joint Committee on Cancer Code
clinic_process$'Neoplasm Disease Stage American Joint Committee on Cancer Code' = factor(clinic$'Neoplasm Disease Stage American Joint Committee on Cancer Code')
table(clinic_process$'Neoplasm Disease Stage American Joint Committee on Cancer Code')

###changed name
table(clinic_process$tumor_stage)
clinic_process$tumor_stage[clinic_process$tumor_stage=='T2A'|clinic_process$tumor_stage=='T2B'] = 'T2'
clinic_process$tumor_stage[clinic_process$tumor_stage=='T3A' | clinic_process$tumor_stage=='T3B'] = 'T3'
clinic_process$tumor_stage = factor(clinic_process$tumor_stage)

clinic_process$tumor_stage[is.na(clinic_process$tumor_stage)] = 'T2' #fill NA with most frequent value

## I want to change into int. 1,3,4 and fill mean in NA

##American Joint Committee on Cancer Tumor Stage Code
clinic_process$`American Joint Committee on Cancer Tumor Stage Code` = factor(clinic_process$`American Joint Committee on Cancer Tumor Stage Code`)
is.na(clinic_process$`American Joint Committee on Cancer Tumor Stage Code`)

##Person Neoplasm Cancer Status
clinic_process$`Person Neoplasm Cancer Status`  = factor(clinic_process$`Person Neoplasm Cancer Status`)
table(is.na(clinic_process$`Person Neoplasm Cancer Status`))
table(clinic_process$`Person Neoplasm Cancer Status`)
clinic_process$`Person Neoplasm Cancer Status`[is.na(clinic_process$`Person Neoplasm Cancer Status`)] = "Tumor Free" #fill the most frequent value

###chagned name
table(clinic_process$neoplasm_stage)
clinic_process$neoplasm_stage
clinic_process$neoplasm_stage[clinic_process$neoplasm_stage=='STAGE IIIA' | clinic_process$neoplasm_stage=="STAGE IIIB" | clinic_process$neoplasm_stage=="STAGE IIIC"] = "STAGE III"

clinic_process$neoplasm_stage[clinic_process$neoplasm_stage=='STAGE IVA' | clinic_process$neoplasm_stage=="STAGE IVB"] = "STAGE IV"
clinic_process$neoplasm_stage = factor(clinic_process$neoplasm_stage)
clinic_process$neoplasm_stage[is.na(clinic_process$neoplasm_stage)] = 'STAGE I'
## Race Category
clinic_process$`Race Category` = factor(clinic_process$`Race Category`)
table(clinic_process$`Race Category`)
is.na(clinic_process$`Race Category`)
clinic_process$`Race Category`[is.na(clinic_process$`Race Category`)] = "White" #fill the most frequent value
save(clinic_process, file='HCC_clinic_data_after_processing.rda')
rm(clinic, HCC_geneCount)

##American Joint Committee on Cancer Metastasis Stage Code
clinic_process$`American Joint Committee on Cancer Metastasis Stage Code` = factor(clinic_process$`American Joint Committee on Cancer Metastasis Stage Code`)
is.na(clinic_process$`American Joint Committee on Cancer Metastasis Stage Code`)
table(clinic_process$`American Joint Committee on Cancer Metastasis Stage Code`)

1. Planning and Data Load
library("MatchIt")

colnames(clinic_process) =c('case', 'patient_code', 'age', 'sex', 'neoplasm_stage', 'tumor_stage', 'status', 'race', 'metastatis_stage')

2. Check initial imbalance
# No matching; constructing a pre-match matchit object
m.out0 <- matchit(sex ~ age + neoplasm_stage + tumor_stage +status + race + metastatis_stage, data = clinic_process, method = NULL, distance = "glm")

# Checking balance prior to matching
summary(m.out0)

3. Matching
# 1:1 NN PS matching w/o replacement
m.out1 <- matchit(sex ~ age + neoplasm_stage + tumor_stage +status + race + metastatis_stage, data = clinic_process, method = "nearest", distance = "glm")

m.out1 <- matchit(sex ~ age   +status + race + metastatis_stage, data = clinic_process, method = "nearest", distance = "glm")

m.out1

4. Assessing the quality of matches
# Checking balance after NN matching
summary(m.out1, un = FALSE)
##n = FALSE to suppress display of the balance before matching for brevity and because we already saw it.

#visualize the distribution of propensity scores - jitter plot
plot(m.out1, type = "jitter", interactive = FALSE)

#visualize the distribution of propensity scores - jitter plot - QQ plot
plot(m.out1, type = "qq", interactive = FALSE,
     which.xs = c("neoplasm_stage", "tumor_stage", 'race'))
plot(summary(m.out1))

5. Trying a different matching specification
# Full matching on a probit PS
m.out2 <- matchit(sex ~ age + race + tumor_stage + metastatis_stage + neoplasm_stage, data = clinic_process, method = "full", distance = "glm", link = "probit")
m.out2

# Checking balance after full matching
summary(m.out2, un = FALSE)

plot(m.out2, type = "qq", interactive = FALSE, which.xs = c("age","metastatis_stage", "neoplasm_stage", "tumor_stage", 'race'))



# Love plot : summarize balance visually
plot(summary(m.out2))


6. Estimating the Treatment effects
m.data2 <- match.data(m.out2)
head(m.data2)

library("lmtest") #coeftest
library("sandwich") #vcovCL

fit1 <- lm(sex ~ age   +status + race + metastatis_stage, data = clinic_process, weights = m.data2$weights)
coeftest(fit1, vcov. = vcovCL, cluster = ~subclass)

m.data2 <- match.data(m.out2)
fit2 <- lm(re78 ~ treat + age + educ + race + married + nodegree + re74 + re75, data = m.data2, weights = weights)
coeftest(fit2, vcov. = vcovCL, cluster = ~subclass)

save(weight_vector, file='weight_vector.rda')


###moditifcation (normalization of clinic data) ###
normalized_clinic_process = clinic_process

#normalization function
norm_minmax <- function(x){(x- min(x)) /(max(x)-min(x))}

#convert tumor_stage
table(normalized_clinic_process$tumor_stage)
is.na(normalized_clinic_process$tumor_stage)
normalized_clinic_process$tumor_stage = as.character(normalized_clinic_process$tumor_stage)
normalized_clinic_process$tumor_stage[normalized_clinic_process$tumor_stage=='T0'] = 0
normalized_clinic_process$tumor_stage[normalized_clinic_process$tumor_stage=='T1'] = 1
normalized_clinic_process$tumor_stage[normalized_clinic_process$tumor_stage=='T2'] = 2
normalized_clinic_process$tumor_stage[normalized_clinic_process$tumor_stage=='T3'] = 3
normalized_clinic_process$tumor_stage[normalized_clinic_process$tumor_stage=='T4'] = 4
normalized_clinic_process$tumor_stage[normalized_clinic_process$tumor_stage=='TX'] = 1 #fill with the most frequent value
normalized_clinic_process$tumor_stage = as.integer(normalized_clinic_process$tumor_stage)
normalized_clinic_process$tumor_stage = norm_minmax(normalized_clinic_process$tumor_stage )

#neoplasm_stage
table(normalized_clinic_process$neoplasm_stage)
normalized_clinic_process$neoplasm_stage= as.character(normalized_clinic_process$neoplasm_stage)
normalized_clinic_process$neoplasm_stage[normalized_clinic_process$neoplasm_stage=='STAGE I'] = 1
normalized_clinic_process$neoplasm_stage[normalized_clinic_process$neoplasm_stage=='STAGE II'] = 2
normalized_clinic_process$neoplasm_stage[normalized_clinic_process$neoplasm_stage=='STAGE III'] = 3
normalized_clinic_process$neoplasm_stage[normalized_clinic_process$neoplasm_stage=='STAGE IV'] = 4
normalized_clinic_process$neoplasm_stage = as.integer(normalized_clinic_process$neoplasm_stage)
normalized_clinic_process$neoplasm_stage = norm_minmax(normalized_clinic_process$neoplasm_stage)
# metastatis_stage : leave it as it is
table(normalized_clinic_process$metastatis_stage)

# age
summary(normalized_clinic_process$age)
unique(is.na(normalized_clinic_process$age)) #check na value
normalized_clinic_process$age = norm_minmax(normalized_clinic_process$age)

#save data
save(clinic_process, normalized_clinic_process, file='clinic_normalized.rda')

# Propensity algorithm
m.normalized <- matchit(sex ~ age + race + tumor_stage + metastatis_stage + neoplasm_stage + status, data = normalized_clinic_process, method = "full", distance = "glm", link = "probit", estimand = "ATE", s.weights = ~SW)
m.normalized

# Checking balance after full matching
summary(m.normalized, un = FALSE)
plot(m.normalized, type = "qq", interactive = FALSE, which.xs = c("age","metastatis_stage", "neoplasm_stage", "tumor_stage", 'race', 'status')) ##QQ plot
plot(summary(m.normalized)) #Love plot

#output
m.normalized.data <- match.data(m.normalized)
save(m.normalized.data, file='m.normalized_weight.rda')

rm(m.normalized, m.normalized.data)


### divide clinic data into two group. (tumor and normal type) ### 2021.04.26
tumorFree_normalized = normalized_clinic_process[clinic_process$status == 'Tumor Free',]
tumor_normalized = normalized_clinic_process[clinic_process$status == 'With Tumor',]

# Propensity algorithm - tumor
m.tumor.normalized <- matchit(sex ~ age + race + tumor_stage + metastatis_stage + neoplasm_stage, data = tumor_normalized, method = "full", distance = "glm", link = "probit", caliper = 0.3, estimand = "ATE")
m.tumor.normalized

# Checking balance after full matching - tumor
summary(m.tumor.normalized, un = FALSE)
plot(m.tumor.normalized, type = "qq", interactive = FALSE, which.xs = c("age","metastatis_stage", "neoplasm_stage", "tumor_stage", 'race')) ##QQ plot
plot(summary(m.tumor.normalized)) #Love plot

# Propensity algorithm - tumor Free
m.tumorFree.normalized <- matchit(sex ~ age + race + tumor_stage + metastatis_stage + neoplasm_stage, data = tumorFree_normalized, method = "full", distance = "glm", link = "probit", caliper = 0.3, estimand = "ATE")
m.tumor.normalized

# Checking balance after full matching - tumor Free
summary(m.tumorFree.normalized, un = FALSE)
plot(m.tumorFree.normalized, type = "qq", interactive = FALSE, which.xs = c("age","metastatis_stage", "neoplasm_stage", "tumor_stage", 'race')) ##QQ plot
plot(summary(m.tumorFree.normalized)) #Love plot

#output
m.tumor.normalized.data <- match.data(m.tumor.normalized)
m.tumorFree.normalized.data <- match.data(m.tumorFree.normalized)
save(m.tumor.normalized.data, m.tumorFree.normalized.data, file='propensity_result_2021.04.26.rda')

rm(m.tumor.normalized.data, m.tumorFree.normalized.data)

#purity
clinic_purity = read_tsv('../data/TCGA_mastercalls.abs_tables_JSedit.fixed-purity.txt')
table(clinic_process$case %in% clinic_purity$array)
save(clinic_purity, file= 'PanCan-tumor purity.rda')

#clinic data - virus infection
clinic_PanCan = read_tsv('../data/lihc_tcga_pan_can_atlas_2018_clinical_data.tsv')
clinic_LIHC = read_tsv('../../../Resources/TCGA-Assembler-2-master/TCGA-Assembler/BiospecimenClinicalData/nationwidechildrens.org_clinical_patient_lihc.txt') 

table(clinic_process$patient_code %in% clinic_LIHC$bcr_patient_barcode)

table(clinic_LIHC$history_hepato_carcinoma_risk_factors)
clinic_hepatitis = select(clinic_LIHC,bcr_patient_barcode, history_hepato_carcinoma_risk_factors)
table(clinic_hepatitis$history_hepato_carcinoma_risk_factors)

##data pre-processing
colnames(clinic_hepatitis) = c('patient_barcode', 'hepatitis') 
table(clinic_hepatitis$hepatitis)
clinic_hepatitis$hepatitis[clinic_hepatitis$hepatitis == 'Hepatitis C|Hemochromatosis'] = 'Hepatitis C'
clinic_hepatitis$hepatitis[clinic_hepatitis$hepatitis == 'Alcohol consumption|Hepatitis B'] = 'Hepatitis B'
clinic_hepatitis$hepatitis[clinic_hepatitis$hepatitis == 'Hepatitis B|Hepatitis C'] = 'both'
clinic_hepatitis$hepatitis[clinic_hepatitis$hepatitis == '[Not Available]'] = NA
clinic_hepatitis_selected = clinic_hepatitis[clinic_hepatitis$patient_barcode %in% clinic_process$patient_code,]
table(clinic_hepatitis_selected$hepatitis)
both Hepatitis B Hepatitis C 
7          95          48
save(clinic_hepatitis_selected, file='clinic_hepatitis_selected.rda')

clinic_process = left_join(x=clinic_process, y=clinic_hepatitis, by = c('patient_code' ='bcr_patient_barcode'))
