library(tidyverse)
#library(CePa)
library(data.table)
library(matrixStats)
library(ggplot2)
library(caret)
library(rminer)
#library(EnsDb.Hsapiens.v79)
#library(doSNOW)
library(deepnet)
#library(QCApro)
library(PRROC)

# register parallel
getDoParWorkers()
getDoParName()
# register parallel
registerDoSNOW(makeCluster(4, type = "SOCK"))
getDoParVersion()


#BiocManager::install('EnsDb.Hsapiens.v79')

set.seed(67)


################################################################################
#=======================>   START HERE FOR MODELS   <===========================

# Load data just for model.
mt_final<-read.csv("ProstateCancer.csv")
rownames(mt_final)<-mt_final$X
mt_final<-mt_final[,c(-1)]
mt_final$Category<-as.factor(mt_final$Category)

# For Prostate cancer.
mt_final<- read.csv("ProstateCancer.csv")
mt_final<-mt_final[,c(-1)]
mt_final$Category<-as.factor(mt_final$Category)


#remove zero variance columns and correlated variables
# Skip this for prostate cancer 
#mt_final_T<-mt_final[,-nearZeroVar(mt_final)] 
mt_final_c<-mt_final[,-findCorrelation(cor(mt_final[,-length(mt_final)]), .8)]
mt_final<- mt_final_c

################################################################################


########### Build SVM model ###########

# Create Train validation 

TrainIndex<-createDataPartition(mt_final$Category, p = 0.7, list = FALSE, times = 1)
Train<-mt_final[TrainIndex,]
Validation<-mt_final[-TrainIndex,]

# Traincontrol
train_control_5 <- trainControl(
  method = "cv",
  number = 6
)
# Train control 10
train_control_10 <- trainControl(
  method = "cv",
  number = 6,
  classProbs = TRUE,
  savePredictions = T,
  summaryFunction = twoClassSummary
)

#### 



####
# SVM model
svm1 <- train(Category ~., data = Train, method = "svmLinear", trControl = trainControl(
  method = "cv",
  number = i
),metric="Accuracy")
svm1
svm2 <- train(Category ~., data = Train, method = "svmPoly", trControl = train_control_5,metric="Accuracy")
svm2
svm3 <- train(Category ~., data = Train, method = "svmRadial", trControl = train_control_5,metric="Accuracy")
svm3

svmt<- train(Category ~., data = Train, method = "svmPoly", trControl = train_control_5,metric="Accuracy",
             tuneGrid = expand.grid(
                                    degree= c(1,2),
                                    scale=c(0.1,0.01,0.001),
                                    C=seq(0, 1, by = 0.1)
                                    ))
svmt
besttune<-svmt$bestTune
# Random Forest model

rf <- train(Category ~., data = Train, method = "rf", trControl = train_control_5,metric='Accuracy')
rf
# variable importance
#varImp(rf)

# XGBTree model

xgb <- train(Category ~., data = Train, method = "xgbTree", trControl = train_control_5,
             metric='Accuracy',preProcess = c('center', 'scale'))
xgb

#pred<-predict(xgb, Validation[,1:554])
#confusionMatrix(pred,Validation$Category)


# Naive Bayes
nb <- train(Category ~., data = Train, method = "naive_bayes", trControl = train_control_5,metric='Accuracy')
nb

# neural net

nnet <- train(Category ~., data = Train, method = "nnet", trControl = train_control_5,preProcess = c('center', 'scale'))
nnet

# logistic regression
glm <- train(as.factor(Category) ~., data = Train, method = "glm", trControl = train_control_5,metric='Accuracy')
glm

# Gradient booseted Tree
gbm<- train(Category ~., data = Train, method = "gbm", trControl = train_control_5,metric='Accuracy')
gbm

# deepnet
dnn<- train(Category ~., data = Train, method = "dnn", trControl = train_control_5,metric='Accuracy')
dnn

# Show results of all models
acc<-c()
model<-list(svm1,svm2,svm3,xgb,rf,nb,glm,dnn)
m_names<-c('SVM_Linear','SVM_Poly','SVM_Radial','XGBoost','RF','NB','GLM','DNN')
for(i in model){
  acc<-c(acc,max(na.rm = T,i$results$Accuracy))
}

# create the plot with accuracy
pl<-data.frame(Model=m_names,Accuracy=acc)
pl %>% ggplot(aes(x=Model,y=Accuracy,fill=Model))+geom_bar(stat = 'identity')+theme_light()


######################## Gene Importance #############################################
##Variable Importance to filter important genes from 3400

varImp(svm1)->vimp
as.data.frame(vimp$importance)->vimp
vimp$gene<-rownames(vimp)
vimp %>% dplyr::filter(NR>50) %>% dplyr::select(gene)->gene_cols
rownames(gene_cols)<-NULL
as.vector(gene_cols)->gene_cols

Train %>% dplyr::select(gene_cols$gene,Category)->Train
Validation %>% dplyr::select(gene_cols$gene,Category)->Validation

# Identify the CV number
cv_num<-c()
cv_acc<-c()
for(i in 2:nrow(Train)){
  cv_acc <- c(cv_acc,max(train(Category ~., data = Train, method = "svmLinear", trControl = trainControl(
    method = "cv",
    number = i
  ),metric="Accuracy")$results$Accuracy))
  cv_num<-c(cv_num,i)
  
}

data.frame(CV=cv_num,Accuracy=cv_acc)->pp
pp %>% ggplot(aes(x=CV,y=Accuracy))+geom_line(color='blue')+geom_point(color='black')+theme_light()

svm1<-train(Category ~., data = Train, method = "svmLinear", trControl = train_control_5,metric="Accuracy")

################################# Step 1 ###########################
# need to add genes one by one and see if they increase or decrease . store the accuracy and then plot.

j<-1
#res<-vector("list", 667)
res<-vector("list", length(Train))
#for(i in 1:667){
for(i in 1:length(Train)-1){  
  res[[j]] <- max(train(Category ~., data = Train[c(1:i,length(Train))], method = "svmLinear",
                        trControl = train_control_5,metric='Accuracy')$results$Accuracy)
  print(j)
  j<-j+1
}

# transform the accuracy data.
unlist(res)->res
#data.frame(Gene=names(Train[-668]),Accuracy=res)->Gene.Importance
data.frame(Gene=names(Train[-length(Train)]),Accuracy=res[1:length(res)-1])->Gene.Importance
#percentage contribution

per_cont<-function(x) ((x-max(svm3$results$Accuracy))/max(svm3$results$Accuracy))
diff_cont<-function(x) ((x-mean(Gene.Importance$Accuracy))/sd(Gene.Importance$Accuracy))

unlist(lapply(Gene.Importance$Accuracy, per_cont))->Gene.Importance$Percentage
unlist(lapply(Gene.Importance$Accuracy, diff_cont))->Gene.Importance$Difference

# plot without sorting.
ggplot(Gene.Importance,aes(x=1:nrow(Gene.Importance),y=Percentage))+ geom_line(color="orange")+
  geom_point(size=0.75,color="blue")+
  geom_hline(yintercept=0, col="red")+
  xlab("Genes")+
  labs(title = "Train Data Metric Evaluation",x="Genes",y="Accuracy",color="Metrics\n")+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))




Gene.Importance %>% arrange(desc(Difference))->reorder.gene

# plot the accuracy and show increase and decrease with addition of genes
ggplot(reorder.gene %>% arrange(desc(Difference)) %>% head(20),aes(x=1:20,y=Difference))+ 
  geom_point(aes(label=Gene),size=2)+
  geom_line(color="blue",size=1)+
  geom_text(aes(label=Gene),hjust=0, vjust=0,angle = 25)+
  labs(x="Input features",
       y="Influence score",color="Metrics\n")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


top_20<-reorder.gene %>% arrange(desc(Difference)) %>% head(20)
res<-c()
for(i in 1:20){
  Train %>% dplyr::select(top_20$Gene[1:i],Category)->temp
  res<-c(res,max(train(Category ~., data = temp, method = "svmLinear",
            trControl = train_control_5,metric='Accuracy')$results$Accuracy))
}

data.frame(Accuracy=res,Names=top_20$Gene)->pp

pp %>% mutate(Names = fct_reorder(Names, Accuracy, .fun='median')) ->pp
ggplot(pp,aes(x=Names,y=Accuracy))+geom_bar(stat='identity',fill='darkblue')+theme_classic()+
  theme(axis.text.x = element_text(angle = 25,vjust = 0.45))

Train %>% dplyr::select(top_20$Gene[1:20],Category)->temp
Validation %>% dplyr::select(top_20$Gene[1:20],Category)->tv
topSvm<-train(Category ~., data = temp, method = "svmLinear",trControl = train_control_10,metric='Accuracy')
evalm(topSvm,gnames = c('RP','NR'))

pred<-predict(topSvm, tv[,1:20])
confusionMatrix(pred,tv$Category)

# ROC curve
evalm(svm1,gnames = c('RP','NR'))

write.csv(reorder.gene,"Gene_Importance_Step1.csv")

############################ Step 2 ###########################

reorder.gene<-read.csv("Gene_Importance_Step1.csv")
reorder.gene %>% arrange(desc(Accuracy))->reorder.gene

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
  temp <- max(train(Category ~., data = Train[c(cols,"Category")], method = "svmLinear",
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

write(cols,"Important Genes_PC.txt")

# Important Genes
cols

# Plot the data
plot.data<-data.frame(Genes=cols,Accuracy=acc)
plot.data$Genes<-with(plot.data, reorder(Genes, Accuracy))

# plot the graph
ggplot(plot.data,aes(x=Genes,y=Accuracy,group=1))+ geom_line(color="black")+
  geom_point()+
  geom_bar(stat='identity',fill="darkblue")+
  geom_hline(yintercept=max(svm3$results$Accuracy), col="red")+
  labs(title = "Train Data Metric Evaluation with important genes",x="Genes",y="Accuracy")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5,size=12))

# plot to show inference score.
res<-c()
for(i in 1:10){
  Train %>% dplyr::select(plot.data$Gene[1:i],Category)->temp
  res<-c(res,max(train(Category ~., data = temp, method = "svmLinear",
                       trControl = train_control_5,metric='Accuracy')$results$Accuracy))
}

data.frame(Accuracy=res,Names=plot.data$Gene)->pp
pp$Names<-with(pp, reorder(Names,typical))
#pp$typical<-acc

ggplot(pp,aes(x=Names,y=typical,group=1))+geom_bar(stat='identity',fill='darkblue',width = 0.5)+
  geom_line(aes(x=Names,y=))+
  geom_point(aes(y=Accuracy),shape=2)
  theme_classic()+
  theme(axis.text.x = element_text(angle = 25,vjust = 0.45))

# Heatmaps
t(Train[c(cols,'Category')])->heat
colnames(heat)<-heat['Category',]
heat[-8,]->heat
as.data.frame(heat)->heat

heatmap(as.matrix(heat))

# prune the train
Train.cut<-Train[c(cols,"Category")]
Validation.cut<-Validation[c(cols,"Category")]


# retrain the model with important cols.
svm.imp<-train(Category ~., data = Train.cut, method = "svmLinear",
               trControl = train_control_10,metric='Accuracy')
svm.imp

svm.imp1<-train(Category ~., data = Train.cut, method = "svmLinear",
               trControl = train_control_10,metric='Accuracy')
svm.imp1

evalm(svm.imp1,gnames = c('NR','RP'))

# Test data prediction
pred<-predict(svm.imp, Validation.cut[,-11])
confusionMatrix(pred,Validation.cut$Category)


# plot to identify best number

Gene.Importance %>% dplyr::filter(Gene %in% plot.data$Genes)->plot.data
plot.data$Gene<-with(plot.data, reorder(Gene, Accuracy))

ggplot(plot.data,aes(x=Gene,y=Accuracy,group=1))+ geom_line(color="black")+
  geom_point()+
  geom_bar(stat='identity',fill="darkblue")+
  geom_hline(yintercept=max(svm3$results$Accuracy), col="red")+
  labs(title = "Train Data Metric Evaluation with important genes",x="Genes",y="Accuracy")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5,size=12))


######### Gene names #########

col_changed<-str_extract(cols,"([0-9a-zA-Z]+)")

data.frame(ensembldb::select(EnsDb.Hsapiens.v79, keys= cols,
                             keytype = "GENEID", columns = c("SYMBOL","GENEID")))->GTemp

GTemp

####################################################################################

write.csv(plot.data,'PC_plot_data.csv')

Truth<-mt_final[c(cols,"Category")]
Truth %>% dplyr::filter(Category=="RP") %>% dplyr::select(-Category)->Truth.Sen
Truth %>% dplyr::filter(Category=="NR") %>% dplyr::select(-Category)->Truth.Insen

med.s<-reshape::melt(lapply(Truth.Sen,median))$value
med.is<-reshape::melt(lapply(Truth.Insen,median))$value


# create truth table
kcalc<-function(r,flag){
  if(flag=="RP"){
    temp<-if_else(Truth.Sen[r]<med.s,0,1)    
  }else if(flag=="NR"){
    temp<-if_else(Truth.Insen[r]<med.is,0,1)    
  } 
  
  return(temp)
}

# define truth table 
Truth.Table.sensitive<-data.frame(C=NA)
Truth.Table.insensitive<-data.frame(C=NA)

# Calculate Sentive truth table
for(i in 1:length(Truth.Sen)){
  temp<-data.frame(kcalc(i,flag = "RP"))
  Truth.Table.sensitive<-cbind(Truth.Table.sensitive,temp)
}
# update colnames
Truth.Table.sensitive<-Truth.Table.sensitive[-1]
colnames(Truth.Table.sensitive)<-names(Truth.Sen)
Truth.Table.sensitive$outcome<-'RP'


# Calculate InSentive truth table
for(i in 1:length(Truth.Insen)){
  temp<-data.frame(kcalc(i,flag = "NR"))
  Truth.Table.insensitive<-cbind(Truth.Table.insensitive,temp)
}

Truth.Table.insensitive<-Truth.Table.insensitive[-1]
colnames(Truth.Table.insensitive)<-names(Truth.Insen)
Truth.Table.insensitive$outcome<-'NR'

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


library(ComplexHeatmap)

ha = rowAnnotation(
  df = data.frame(Outcome=Truth_table$outcome),
  annotation_height = unit(4, "mm"),
  show_annotation_name = TRUE,
  col= list(Outcome = c('RP' = "Green", 'NR' = "darkorange"))
)

Heatmap(Truth_table[cols],
        border = T,rect_gp = gpar(col='Black'),
        name = "Expression",
        right_annotation=ha,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10,font=3)
        )

Heatmap(Truth[cols],
        border = T,rect_gp = gpar(col='Black'),
        name = "Expression",
        right_annotation=ha,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10,font=3)
)



# all the data.
Heatmap(mt_final[-3793],km=2,
        show_column_names = FALSE,
        name = "Expression",
        right_annotation=ha,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10,font=3)
)




################################################################################

