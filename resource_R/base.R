#sound
library(beepr)
beep() # alarm function

#Basic command key in Rstudio
Comment command : Ctrl + shift + C
(consol) Ctrl + C

#library
library(readr)
library(ggplot2)
library(tidyverse)
library(readxl)
library("Hmisc") #correlation matrix with p-value
library(corrplot) #visualizing correlation matrix
library(beepr) #beep() : alarm function
library(DESeq2)
library("apeglm")
library(ggraph)
library(igraph)

#variable
class(variable)
dim()
"AsIs" : dataframe class, which does not allow to modify it
rm(variable)
length(list)
ls() : show all variable names
rm(list = ls())
rm()

a <-1 : assignment
a < -1 : boolian

# object identification
typeof(obj)
class(obj)
is.vector(object)

sapply(obj, class)
sapply(obj, attributes)
attributes(obj)
names(obj)

<boolian operation>
vector %in% vector
e.g. 1 %in% c(1,2,3) #TRUE

duplicated(vector, data.frame)
e.g. duplicated(data.frame)
#random selection
sample(vector, the number of selection)
e.g. sample(c(1,2,3,4,5,6), 2))

<factor>
#create vector containing factor type
vector = factor(c(v1, v2, ))

#change the name of factor
library(plyr)
new.vector = revalue(vector, c(name1_old = name1_new, name2_old = name2_new)) #change name of factor

<conversion of data type>
#list -> matrix -> data.frame
data.frame = data.frame(matrix = matrix(unlist(input_list), nrow=length(list), ncol=2, byrow=TRUE) 
                        
#transpose data.frame into matrix
matrix = t(data.frame)

#dataframe row -> vector
vector = as.numeric(dataframe)
#dataframe column -> vector
vector = as.vector(dataframe)
vector = dataframe$colname
                        
<string>
#example
LETTERS[1:5]

#select the substring
new string= substr(data.frame$column, start = 6, stop = 7)

#find and change in string
conversion_vector = gsub("[.]","-", input_vector) #change '.' into '-'
e.g. conversion = gsub("[.]","-", GTEx_sample_id$conversion)
e.g. TCGA id converter



<vector>
#random samples
e.g. sample(x=1:30, size = 5)
        
#merge vector
nerged_vectir = c(vector1, vector2)

#make a vector with null space
vector <- rep(NA, length)
e.g. a <- rep(NA, 10)

#count the number of TRUE element
length(vector[vector == condition])
e.g.length(z[z == TRUE])

#named numeric vector
names(vector) = ("name1", "name2")

# find a location of specific value
where(vector == "name")

<data frame>
#initialize new dataframe                        
new_df = data.frame(matrix(nrow=3, ncol=4))                        
new_df = data.frame(matrix(0, nrow=3, ncol=4))   
                        
#data frame indexing
dataframe[] #data frame
dataframe[[]] #vector
                       

#Explore data frame
head( dataframe )
str(datafrane) #structure of dataframe 
rowData( dataframe )
colData( sdataframe )

#drop a column
new.data.frame = subset(data.frame, select = -c(1))

#drop a NA in specified column
new.dataframe=drop_na(old.dataframe, specified column name) 
ex. new.dataframe=drop_na(old.dataframe, Sex) 
# merge by row & add single observation into data.frame
new_df = rbind(old_df, add_df)
merged_df = merge(dataframe1, dataframe2, by=0) #merge data according to row naems

#change name of column
rename.variable(df, old name of column, new name of column)
e.g. iris <- rename.variable(iris, "Species", "especes")

colnames(data.frame) = c("name1", "name2")
rownames(data.frame) = c("name1", "name2")

#add new column at the most front 
new_column = data.frame(sample = rownames(GTEx_liver))
merged_df = cbind(new_column, data.frame) #add sample column at first

#extract into new data.frame
new_data.frame=select(data.frame, colum name)
new_data.frame=filter(data.frame, boolian for all row)

#find out wheter the specific substring exists in data.frame (version 1)
    #detect mathces
function_name = function(input_char){
    for(i in 1:4){
        if(TRUE){
            execution line
        }
        else if{
            exeuction line
        }
        else{
            execution line
        }
    }
    
    while(TRUE){
        execution line
    }
    return(str_detect(input_char, "stiring needed to be detected"))
}
result = sapply(HMR$EQUATION, function_name) #HMR$EQUATION : vector

#find out wheter the specific substring exists in data.framee (version 2)
Boolian list = str_detect(data,frame, 'string need to be detected') #find the index (boolian)

#split the dataframe into two group by 'string'
list = strsplit(data.frame$column, 'string')  #convert dataframe into list 
vector = unlist(list) # convert list into vector

#change the order of column
new data.frame = old data.frame[,c(sequence)]
e.g. df2[,c(1,3,2,4)]

#merge dataframe
bind_rows(df1, df2)
bind_cols(df1,df2)

new.dataFrame = left_join(x=left_df, y=right_df, by = c('merged_left_df_column' ='merged_right_df_column'))

#convert char dataframe into numeric data.frame
numeric_data.frame = data.frame(lapply(char_data.frame, as.numeric))
data.frame$column = unlist(lapply(data.frame$column, as.numeric))

#the result of lapply is a list
vector = unlist(list) #convert vector into list
    #convert char dataframe into numeric dataframe.
rownames(numeric_data.frame) = rownames(char_data.frame)

#summary
table(vector or data.frame)


#multiply column by column with vector v
data.frame(mapply('*',df,v))

#summary with package 'dplyr'
#statistics
statistics = data.frame %>%
    group_by(colnames) %>%
    summarise(number=n())

Male tumor_free 94
male tumor 32
female tumor_free 60
female tumor 33

rm(statistics)

<matrix>
#the number of column and row
nrow(vector, array, dataframe, matrix)
ncol(vector, array, dataframe, matrix)

#merge two matrinx into one matrix
merged_matrix = merge(matrix_left, matrix_right, by = "row.names", all = TRUE)

<function>
??functioin_name
help(function_nam)






<File I/O>
#load the data
read_delim()
read_tsv(filename) #readr package
df <- read.csv(file="my.large.file.csv",nrows=2000)
    
    
#save and load in rda file
save(variable, file = "mydata.rda")
load(file = "mydata.rda")

#save and load of dataframe in txt format
write(vector, "filename.txt") #export vector into txt file

input_data = read.table("C:/R/default_working_directory/Secretome/BRCA Project/output/dataDEGs.txt", header =TRUE, sep="")

write.table(result, paste0(getwd(), "/conversion.txt"), quote=FALSE, append=FALSE, row.names = FALSE)

write.table(dataframe, file = "male_tumor_only.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote=FALSE, append=FALSE)

#sort
data.frame[sort(vector / df$colnames),] #sort row

#correlation
matrix = cor(correlation_prepared_data.frame) 
    row : gene
    column : sampple
    
    


