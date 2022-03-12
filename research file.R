library(tidyverse)
library(CePa)
library(data.table)
library(matrixStats)
library(ggplot2)
library(caret)
library(rminer)
library(EnsDb.Hsapiens.v79)
library(doSNOW)
library(deepnet)
library(QCApro)

# register parallel
getDoParWorkers()
getDoParName()
# register parallel
registerDoSNOW(makeCluster(4, type = "SOCK"))
getDoParVersion()


#BiocManager::install('EnsDb.Hsapiens.v79')

set.seed(67)
############ change the data file
read.gct('CCLE_RNAseq_genes_rpkm_20180929.gct')->CCLE
as.data.frame(CCLE)->CCLE
write.csv(CCLE,"ccle-full.csv")


read.csv("GDSC.csv")->GDSC

########### Read files ###########

read.csv("GDSC.csv")->GDSC
fread("ccle-full.csv")->CCLE

########### Check for sensitive and insensitive #############


gdsc.median<- median(GDSC$IC50)
GDSC$Category<-if_else(GDSC$IC50<gdsc.median,'Sensitive','Insensitive')

############ Extract Cell lines ############

names(CCLE)->cellcols
str_extract(cellcols,"(\\w+?)_")->cellcols
str_extract(cellcols,"([0-9a-zA-Z]+)")->cellcols
cellcols[1]<-"ID"
cellcols[2]<-"Gene"

names(CCLE)<-cellcols


################################################################################

########### Prune Columns ###########

# Extract common cell lines from GDSC and CCLE
keep<-intersect(GDSC$Cell_line,names(CCLE))

# Obtain 
CCLE[,..keep]->CCLE.prune

# Add the Gene ID
cbind(CCLE[,"Gene"],CCLE.prune)->gene.cell

########### T-TEST ###########

# Store only IC50 and celline from GDSC
rho.gdsc<-GDSC %>% filter(Cell_line %in% keep) %>% select(Cell_line,IC50)

# Store only Cellline and Category from GDSC
cat.gdsc<-GDSC %>% filter(Cell_line %in% keep) %>% select(Cell_line,Category)

################################################################################

# create function for t-test pvalue and Correlation values 

ttest.rho<-function(x){
  
  # transpose the row
  conv.tab<-melt(x,id.vars = "Gene")
  
  # convert into data frame
  as.data.frame(conv.tab)->conv.tab
  
  # filter 0 values
  conv.tab %>% filter(value!=0)->conv.tab
  
  # merge and fetch category information for T-test
  conv.tab.T<-merge(conv.tab,cat.gdsc,by.x="variable",by.y="Cell_line")
  
  # sensitive group
  sens<-conv.tab.T %>% filter(Category == 'Sensitive') %>% pull(value)
  
  # insensitive group
  insens<-conv.tab.T %>% filter(Category == 'Insensitive') %>% pull(value)
  
  # result
  result<-c()
  
  # T-Test
  if(length(sens)==0|| length(insens)==0){
    
    result[1]<-NA
    
  }else{
    
    # Catch the ttest error of low columns.
    tryCatch(
      {
        # obtain the pvalue 
        value<-t.test(sens,insens,var.equal=T)
        # return pvalue
        result[1]<-value$p.value
      },
      error=function(cond) {
        result[1]<-NA
        }
      
    )
    
  }
  
  # Correlation coefficient for gene
  # Merge two dataframes
  merge.tab<-merge(conv.tab,rho.gdsc,by.x="variable",by.y="Cell_line")
  
  # filter the 0 values from correlation
  merge.tab %>% filter(value!=0)->merge.tab
  
  # Spearman correlation or RHO
  coeff<-cor(merge.tab$value, merge.tab$IC50, method = "pearson")
  
  result[2]<-coeff
  
  return(result)
  
}
################################################################################

########### Evaluate the ttest and correlation values ###########

# create temporary datatable for pvalues.
temp.p<-c()
temp.r<-c()

# calculate the Pvalue using ttest
for(i in 1:dim(gene.cell)[1]){
  temp<-gene.cell[i,ttest.rho(.SD)]
  temp.p[i]<-temp[1]
  temp.r[i]<-temp[2]
}

# assign the p-values back to the gene.cell
gene.cell[,Pvalue:=temp.p]

# assign the RHO back to the gene.cell
gene.cell[,Rho:=temp.r]

# Remove Na values and store to new Data.table
# Set Threshold for Pvalue as 0.001
GCT<-gene.cell[!is.na(Pvalue)]

# Store the file for future reference.
fwrite(GCT,"Gene_Cell_Ttest_Pvalue.csv")


################################################################################

###########  Volcano plot. ########### 

as.data.frame(GCT)->GCT

ggplot(data=GCT, aes(y=-log10(Pvalue), x= abs(Rho))) + 
  geom_point() + 
  theme_minimal()+
  geom_vline(xintercept=0.2, col="red") +
  geom_hline(yintercept=-log10(0.01), col="red")

################################################################################

# Prune the data based on the threshold for correlation and -log10(Pvalue)

GCT %>%  filter(abs(Rho)>0.2 & Pvalue<0.01) %>% select(-Pvalue,-Rho) -> final
rownames(final)<-final$Gene
final %>% select(-Gene)->final

# check dimension of final data.
dim(final)

write.csv(final,"GDSC_CCLE_FINAL.csv")

################################################################################

############ Transpose the dataframe ############

t_final <- transpose(final)
colnames(t_final) <- rownames(final)
rownames(t_final) <- colnames(final)

# store Gene as a column to merge and obtain category
t_final$Gene<-rownames(t_final)
# perform merge
merge(t_final,cat.gdsc,by.x="Gene",by.y="Cell_line",all.x = TRUE)->mt_final

# Remove duplicated Cell Lines 
mt_final[!duplicated(mt_final$Gene),]->mt_final

# Set Rownames as gene and remove gene column
rownames(mt_final)<-mt_final$Gene

# drop gene column
mt_final %>% select(-Gene)->mt_final
mt_final$Category<-as.factor(mt_final$Category)

write.csv(mt_final,"Train.csv")

################################################################################
#=======================>   START HERE FOR MODELS   <===========================

# Load data just for model.
mt_final<-read.csv("Train.csv")
rownames(mt_final)<-mt_final$X
mt_final<-mt_final[,c(-1,-2)]
mt_final$Category<-as.factor(mt_final$Category)

# For Prostate cancer.
#mt_final<- read.csv("ProstateCancer.csv")
#mt_final<-mt_final[,c(-1)]
#mt_final$Category<-as.factor(mt_final$Category)


#remove zero variance columns and correlated variables
# Skip this for prostate cancer 
mt_final_T<-mt_final[,-nearZeroVar(mt_final)] 
mt_final_c<-mt_final_T[,-findCorrelation(cor(mt_final_T[,-length(mt_final_T)]), .8)]
mt_final<- mt_final_c

################################################################################


########### Build SVM model ###########

# Create Train validation 

TrainIndex<-createDataPartition(mt_final$Category, p = 0.8, list = FALSE, times = 1)
Train<-mt_final[TrainIndex,]
Validation<-mt_final[-TrainIndex,]

# Traincontrol
train_control_5 <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 3
)
# Train control 10
train_control_10 <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 3
)

# SVM model
svm1 <- train(Category ~., data = Train, method = "svmLinear", trControl = train_control_5,metric="Accuracy")
svm1
svm2 <- train(Category ~., data = Train, method = "svmPoly", trControl = train_control_5,metric="Accuracy")
svm2
svm3 <- train(Category ~., data = Train, method = "svmRadial", trControl = train_control_5,metric="Accuracy")
svm3

svmt<- train(Category ~., data = Train, method = "svmRadial", trControl = train_control_5,metric="Accuracy",
             tuneGrid = expand.grid(C = c(0.25,0.5,1),sigma=1/seq(20,50,length=20)^2))
svmt
besttune<-svmt$bestTune
# Random Forest model

rf <- train(Category ~., data = Train, method = "rf", trControl = train_control_5,metric='Accuracy')
rf
# variable importance
#varImp(rf)

# XGBTree model

xgb <- train(Category ~., data = Train, method = "xgbTree", trControl = train_control_5,
             metric='Accuracy')
xgb

pred<-predict(svm3, Validation[,1:554])
confusionMatrix(pred,Validation$Category)

# Naive Bayes
nb <- train(Category ~., data = Train, method = "naive_bayes", trControl = train_control_5,metric='Accuracy')
nb

# neural net

nnet <- train(Category ~.-Category, data = Train, method = "nnet", trControl = train_control_5,metric='Accuracy')
nnet

# logistic regression
glm <- train(Category ~., data = Train, method = "glm", trControl = train_control_5,metric='Accuracy')
glm

# Gradient booseted Tree
gbm<- train(Category ~., data = Train, method = "gbm", trControl = train_control_5,metric='Accuracy')
gbm

# deepnet
dnn<- train(Category ~., data = Train, method = "dnn", trControl = train_control_5,metric='Accuracy')


.######################## Gene Importance #############################################

## Identify Gene Importance using Looping ##
#nb_gene <- train(Category ~., data = Train, method = "naive_bayes", trControl = train_control_5,metric='Accuracy')
#mean(nb_gene$results$Accuracy)

############################# using loop to identify inference score #############################
j<-1
res<-vector("list", length(Train))
for(i in 1:length(Train)){
  
  nb_gene <- train(Category ~., data = Train[-i], method = "svmRadial", trControl = train_control_5,metric='Accuracy')
  res[[j]]<-max(nb_gene$results$Accuracy)
  print(j)
  j<-j+1
}

#unlist the result and store to data frame
unlist(res)->res
data.frame(Gene=names(Train[-length(Train)]),Result=res)->Gene.Importance

# plot the importance
Gene.Importance %>% 
  ggplot(aes(x=Gene,y=Result))+ 
  geom_point()+
  geom_hline(yintercept=max(nb$results$Accuracy), col="red")

#percentage contribution

per_cont<-function(x) ((x-mean(nb$results$Accuracy))/mean(nb$results$Accuracy))*100

unlist(lapply(Gene.Importance$Result, per_cont))->Gene.Importance$Percentage

# filter the genes that don't contribute

Gene.Importance %>% filter(round(Result,3)> 0.64 | round(Result,3)< 0.63)->Gene.percent

# using variable importance
varImp(nb)$importance %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  arrange(desc(Sensitive)) %>% 
  mutate(Gene=forcats::fct_inorder(rowname ),Importance=Sensitive) %>% 
  select(Gene,Importance)->Gene.var.imp

# check the correlation
merge(Gene.percent,Gene.var.imp,by.x = "Gene",by.y = "Gene")->imp.vs.inf

# plot graph

imp.vs.inf %>% ggplot(aes(x=Percentage,y=Importance))+ geom_point()

# filter important genes 
Train %>% select(Gene.percent$Gene,Category) ->Train.cut


################################# Step 1 ###########################
# need to add genes one by one and see if they increase or decrease . store the accuracy and then plot.

j<-1
#res<-vector("list", 667)
res<-vector("list", length(Train))
#for(i in 1:667){
for(i in 1:length(Train)-1){  
  res[[j]] <- max(train(Category ~., data = Train[c(1:i,length(Train))], method = "svmRadial",
                        trControl = train_control_5,metric='Accuracy')$results$Accuracy)
  print(j)
  j<-j+1
}

# transform the accuracy data.
unlist(res)->res
#data.frame(Gene=names(Train[-668]),Accuracy=res)->Gene.Importance
data.frame(Gene=names(Train[-length(Train)]),Accuracy=res[1:length(res)-1])->Gene.Importance
#percentage contribution

per_cont<-function(x) ((x-max(svm3$results$Accuracy))/max(svm3$results$Accuracy))*100

unlist(lapply(Gene.Importance$Accuracy, per_cont))->Gene.Importance$Percentage


# plot without sorting.
ggplot(Gene.Importance,aes(x=1:nrow(Gene.Importance),y=Accuracy))+ geom_line(color="orange")+
  geom_point(size=0.75,color="blue")+
  geom_hline(yintercept=max(svm3$results$Accuracy), col="red")+
  xlab("Genes")+
  labs(title = "Train Data Metric Evaluation",x="Genes",y="Accuracy",color="Metrics\n")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


Gene.Importance %>% arrange(Accuracy)->reorder.gene

# plot the accuracy and show increase and decrease with addition of genes
  ggplot(reorder.gene %>% arrange(Accuracy),aes(x=1:nrow(Gene.Importance),y=Accuracy))+ geom_line(color="blue")+
  geom_point(size=0.5)+
  geom_hline(yintercept=max(svm1$results$Accuracy), col="red")+
  labs(title = "Train Data Metric Evaluation",x="Genes",y="Accuracy",color="Metrics\n")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  
write.csv(reorder.gene,"Gene_Importance_Step1.csv")

############################ Step 2 ###########################

reorder.gene<-read.csv("Gene_Importance_Step1.csv")
reorder.gene %>% arrange(desc(Accuracy)) %>% select(-X)->reorder.gene

set.seed(56)

# add gene from the ascending order. add one by one and see if 
# there is a 2% increase to accuracy only then include that gene into the model.
# if there is no increase remove the gene and move on to the next one.

thershold<-0
cols<-c()
acc<-c()
#for(i in 1:667){
for(i in 1:(length(Train)-1)){  
  cols<-c(cols,reorder.gene$Gene[i]) 
  temp <- max(train(Category ~., data = Train[c(cols,"Category")], method = "svmRadial",
                         trControl = train_control_5,metric='Accuracy')$results$Accuracy)
  
  # check for 2% increase in the accuracy
  if(temp>((0.01*thershold)+thershold))
  {
    thershold<-temp
    acc<-c(acc,thershold)
    print(thershold)
    
  }else{
    cols<-cols[1:length(cols)-1]
  }

}

########################### Retrain with important genes ###########################

write(cols,"Important Genes.txt")

# Important Genes
cols

# Plot the data
plot.data<-data.frame(Genes=cols,Accuracy=acc)
plot.data$Genes<-with(plot.data, reorder(Genes, Accuracy))

# plot the graph
ggplot(plot.data,aes(x=Genes,y=Accuracy,group=1))+ geom_line(color="orange")+
  geom_point(size=1,color="blue")+
  geom_hline(yintercept=max(svm3$results$Accuracy), col="red")+
  labs(title = "Train Data Metric Evaluation with important genes",x="Genes",y="Accuracy")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# prune the train
Train.cut<-Train[c(cols,"Category")]
Validation.cut<-Validation[c(cols,"Category")]


# retrain the model with important cols.
svm.imp<-train(Category ~., data = Train.cut, method = "svmRadial",
      trControl = train_control_5,metric='Accuracy')
svm.imp



# Test data prediction
pred<-predict(svm.imp, Validation[,1:552])
confusionMatrix(pred,Validation$Category)

######### Gene names #########

col_changed<-str_extract(cols,"([0-9a-zA-Z]+)")

data.frame(ensembldb::select(EnsDb.Hsapiens.v79, keys= cols,
                              keytype = "GENEID", columns = c("SYMBOL","GENEID")))->GTemp

GTemp

####################################################################################

Truth<-mt_final[c(cols,"Category")]
Truth %>% dplyr::filter(Category=="Sensitive") %>% dplyr::select(-Category)->Truth.Sen
Truth %>% dplyr::filter(Category=="Insensitive") %>% dplyr::select(-Category)->Truth.Insen

med.s<-reshape::melt(lapply(Truth.Sen,median))$value
med.is<-reshape::melt(lapply(Truth.Insen,median))$value


# create truth table
kcalc<-function(r,flag){
  if(flag=="Sen"){
    temp<-if_else(Truth.Sen[r]<med.s,0,1)    
  }else if(flag=="Insen"){
    temp<-if_else(Truth.Insen[r]<med.is,0,1)    
  } 
  
  return(temp)
}

# define truth table 
Truth.Table.sensitive<-data.frame(C=NA)
Truth.Table.insensitive<-data.frame(C=NA)

# Calculate Sentive truth table
for(i in 1:length(Truth.Sen)){
  temp<-data.frame(kcalc(i,flag = "Sen"))
  Truth.Table.sensitive<-cbind(Truth.Table.sensitive,temp)
}
# update colnames
Truth.Table.sensitive<-Truth.Table.sensitive[-1]
colnames(Truth.Table.sensitive)<-names(Truth.Sen)
Truth.Table.sensitive$outcome<-1


# Calculate InSentive truth table
for(i in 1:length(Truth.Insen)){
  temp<-data.frame(kcalc(i,flag = "Insen"))
  Truth.Table.insensitive<-cbind(Truth.Table.insensitive,temp)
}

Truth.Table.insensitive<-Truth.Table.insensitive[-1]
colnames(Truth.Table.insensitive)<-names(Truth.Insen)
Truth.Table.insensitive$outcome<-0

# final table with sensitive 1 and insenstive 0.
Truth_table<-rbind(Truth.Table.sensitive,Truth.Table.insensitive)


# Use Enchanced Quine-McCluskey Algorithm
QMO<-eQMC(Truth_table,outcome = 'outcome')
# display the truth table
print(QMO$tt)

# display equation
QMO

# write to file
capture.output(QMO,file = "truthtable_equation.txt")
capture.output(QMO$tt,file = "truthtable_output.txt")


################################################################################

