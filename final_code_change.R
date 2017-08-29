rm(list = ls())
setwd("D:/shahid/project")
getwd()

#header files
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(xgboost)
library(caret)
library(proxy)
library(logisticPCA)
library(rARPACK)
library(ica)
library(tibble)
library(stringr)
library(forcats)
library(corrplot)
library(lubridate)
library(mlr)
library(data.table)
library(MASS)
library(class)
library(e1071)
library(tm)
library(SnowballC)
library(Matrix)
library(syuzhet)
library(nnet)
library(gsubfn)
library(reshape2)
library(devtools)
library(qdapTools)

# LabelCount Encoding function
labelCountEncoding = function(column){
  return(match(column,levels(column)[order(summary(column,maxsum=nlevels(column)))]))
}



# Load CSV files
train_variants = read.csv("train_variants_csv.csv", header=T)
test_variants = read.csv("test_variants_csv.csv", header=T)

#load text file
"Read data"
train_text = do.call(rbind,strsplit(readLines('training_text'),'||',fixed=T))
train_text = as.data.table(train_text)
train_text = train_text[-1,]
colnames(train_text) = c("ID", "Text")
train_text$ID = as.numeric(train_text$ID)

test_text = do.call(rbind,strsplit(readLines('test_text'),'||',fixed=T))
test_text = as.data.table(test_text)
test_text = test_text[-1,]
colnames(test_text) = c("ID", "Text")
test_text$ID = as.numeric(test_text$ID)

#join train variants and train text by "ID"
train_data = merge(train_variants,train_text,by="ID")
test_data = merge(test_variants,test_text,by="ID")
rm(test_text,train_text,train_variants,test_variants)

train_data_class = train_data$Class
train_data$Class = NULL

#merge train data and test data
data = rbind(train_data, test_data)
rm(train_data,test_data)

#text features ##character/word ###
"Basic text features"
##create separate variable to count number of characters in a text
data$nchar = as.numeric(nchar(data$Text))
#create separate variable to count number of characters in a text
data$nwords = as.numeric(str_count(data$Text, "\\S+"))
#rm(test_txt,test_txt_dump,testfull,train_txt,train_txt_dump,trainfull,testing_variants,training_variants)

# for TF-IDF 
#create corpus
data_text = Corpus(VectorSource(data$Text))
#remove extra leading and trailing spaces
data_text = tm_map(data_text, stripWhitespace)
#convert to lower case
data_text = tm_map(data_text, content_transformer(tolower))
writeLines(as.character(data_text[1]))

#remove functuation marks
data_text = tm_map(data_text, removePunctuation)
writeLines(as.character(data_text[1]))
#remove stopwords
data_text = tm_map(data_text, removeWords, stopwords("english"))
writeLines(as.character(data_text[1]))
#convert to root form/base form
data_text = tm_map(data_text, stemDocument, language="english")
writeLines(as.character(data_text[1]))
#remove numbers
data_text = tm_map(data_text, removeNumbers)
writeLines(as.character(data_text[1]))
#create document term matrix
doc_term_mat = DocumentTermMatrix(data_text, control = list(weighting = weightTfIdf))
doc_term_mat = removeSparseTerms(doc_term_mat, 0.95)

# Create dataframe
data = cbind(data, as.matrix(doc_term_mat))

#symbol/letter extraction from variation
#extracting first letter from variaton for typical variations
data = extract(data,Variation, into=c("First_Letter", "Var_2"), 
          regex = "^(?=.{1,7}$)([a-zA-Z]+)([0-9].*)$", 
          remove = FALSE)

# extract number and last letter for typical variations
data = extract(data,Var_2, into=c("Gene_Location", "Last_Letter"),
          regex = "^([0-9]+)([a-zA-Z]|.*)$",
          remove = TRUE)

# Identify and encode deletions
data$is_deletion = ifelse(grepl("del", data$Variation, ignore.case = T), 1, 0)

# Identify and encode insertions
data$is_insertion = ifelse(grepl("ins", data$Variation, ignore.case = T), 1, 0)

# Identify and encode fusions
data$is_fusion = ifelse(grepl("fus", data$Variation, ignore.case = T), 1, 0)

# Identify and encode truncation
data$is_truncation = ifelse(grepl("trunc", data$Variation, ignore.case = T), 1, 0)

# Identify and encode methylation
data$is_methylation = ifelse(grepl("methyl", data$Variation, ignore.case = T), 1, 0)

# Identify and encode amplification
data$is_amplification = ifelse(grepl("amp", data$Variation, ignore.case = T), 1, 0)

# Identify and encode silencing
data$is_silencing = ifelse(grepl("sil", data$Variation, ignore.case = T), 1, 0)

# Identify and encode overexpression
data$is_ov_expression = ifelse(grepl("expr", data$Variation, ignore.case = T), 1, 0)

# Identify and encode splicing
data$is_splicing = ifelse(grepl("splice", data$Variation, ignore.case = T), 1, 0)

# Identify and encode exon variations
data$is_exon_var = ifelse(grepl("exon", data$Variation, ignore.case = T), 1, 0)

# One hot encode first letter variables
first_letter_pattern = mtabulate(data$First_Letter)
colnames(first_letter_pattern) = paste("first_", colnames(first_letter_pattern), sep="")
first_letter_pattern = subset(first_letter_pattern, select=-c(first_CASP, first_DNMT))

# One hot encode last letter variables
last_letter_pattern = mtabulate(data$Last_Letter)
colnames(last_letter_pattern) = paste("last_", colnames(last_letter_pattern), sep="")
last_letter_pattern = subset(last_letter_pattern, select=-c(last_BRAF, last_del, last_dup))

# Identify and encode normal variations 
last_letter_pattern$normal_vars = rowSums(last_letter_pattern[,])

# Bind first and last letter vars back in
data = cbind(first_letter_pattern, data)
data = cbind(last_letter_pattern, data)
rm(first_letter_pattern,last_letter_pattern)
# Identify and encode named variations insertion,deletion,truncation ....
data$named_variation = rowSums(data[,c("is_deletion", "is_insertion", "is_fusion", "is_truncation", "is_methylation", "is_amplification", "is_silencing", "is_ov_expression", "is_splicing", "is_exon_var")])

# Identify and encode wierd variations which is not above
data$weird_variation = rowSums(data[,c("normal_vars", "named_variation")])
data$weird_variation = ifelse(data$weird_variation == 0, 1, 0)

# Change  variation and Gene to factors for label encoding
data$Variation = as.factor(data$Variation)
data$Gene = as.factor(data$Gene)

# Label count encoding
data$Variation = labelCountEncoding(data$Variation)
data$Gene = labelCountEncoding(data$Gene)

# Sentiment Analysis
sentiment_rec = get_nrc_sentiment(data$Text) 
#combine sentiment with data
data = cbind(data,sentiment_rec)
rm(sentiment_rec)

# delete character variables
data = subset(data, select=-c(First_Letter, Last_Letter))

#if any n.a?? make it to 0
data[is.na(data)] = 0

# train/test data
train_data = data[1:3321,]
test_data = data[3322:8989,]

train_ID = train_data$ID
train_data_class = train_data_class - 1
train_data$ID = NULL
train_data$txt = NULL

test_ID = test_data$ID
test_data$ID = NULL
test_data$txt = NULL

# Make data numeric
train_data[] = lapply(train_data, as.numeric)
test_data[] = lapply(test_data, as.numeric)

# Create matrices
model_train_dat = xgb.DMatrix(Matrix(as.matrix(train_data), sparse = TRUE), label = train_data_class)
test_dat = xgb.DMatrix(Matrix(as.matrix(test_data), sparse = TRUE))

#model building using XGBoost parameters
#parameter list
parameters = list(booster = "gbtree",
                  eval_metric = "mlogloss", 
                  objective = "multi:softprob",
                  min_child_weight = 0.67,
                  colsample_bytree = 0.65,                       
                  subsample = 1,                              
                  max_depth =6,                           
                  eta = 0.04,                             
                  gamma = 1,
                  num_class = 9)


# Find optimal nrounds 
cv_para = xgb.cv(parameters, model_train_dat, 
                     early_stopping_rounds = 40, 
                     nfold = 10,
                     maximize = FALSE,
                     prediction = TRUE,
                     nrounds=2000)

# Train model with optimal number of rounds
model_train = xgb.train(parameters, model_train_dat, nrounds = 179)

# calculate importance of variables
importance_mat = xgb.importance(model = model_train)
x=importance_mat[1:20,]
xgb.plot.importance(importance_matrix = x)

# Predict on test data 
test_predict = predict(model_train, test_dat)
test_predict = t(matrix(test_predict, nrow=9, ncol=length(test_predict)/9))
test_predict = as.data.frame(test_predict)
names(test_predict)[1:9] = c("class1", "class2", "class3", "class4", "class5", "class6", "class7", "class8", "class9")
test_predict = cbind(test_ID, test_predict)
names(test_predict)[1] = "ID"
test_predict$ID = as.integer(test_predict$ID)

#write result
 write_csv(test_predict,'predictions_result.csv')
