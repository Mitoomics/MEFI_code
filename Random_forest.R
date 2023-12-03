rm(list=ls())
path<-'G:Model construction\\your_file_path'
setwd(path)

library(randomForest)
train<-read.table(file='train.txt',sep='\t',header=T,fileEncoding="GBK")
train$label<-as.factor(train$label)


#ģÐ͵Ĺٽ¨ 
rf_model<-randomForest(label ~ . , data = train,
                                   #ntree = 500,
                                   #mtry = 3,
                                   importance = TRUE,
                                   proximity = TRUE)

importance<-rf_model$importance
importance<-cbind(rownames(importance),importance)
write.table(importance,"randomForest importance.txt",sep="\t",quote = F,row.names = F,col.names=T)

varImpPlot(rf_model, main = "variable importance")
savePlot(filename='randomForest importance.jpg',type='jpg')
dev.off()

test<-read.table("test.txt",header = T,sep = "\t")
predicted.result<-as.character(predict(rf_model,test))
predicted.score<-predict(rf_model,test,type = "prob")
res<-cbind(test,predicted.result,predicted.score)
colnames(res)<-c(colnames(test),'predicted.result',paste0(colnames(predicted.score),'(predicted.prob)'))
write.table(res,"randomForest predicted_result.txt",sep="\t",quote = F,row.names = F)
