# library(glmnet)
#
# rm(list=ls())
# data<-read.csv("C:/Projects/R packages/corteva/users/vahe/data.csv")
# dim(data)
# #Remove the missing observations
# data<-data[!is.na(data[,350]),]
#
# #Set of environmental covariates
# V<-data[, -c(1,2,3)]
# V[is.na(V)]<-0
# p<-dim(V)[2]
# n<-dim(V)[1]
# V<-V[,1:p]
# dim(V)
#
# #Response Variable
# Y<-data$yield
#
# #Indexes of Genotypes
# A1<-matrix(rep(0, n*22), ncol=22)
# t<-1
# for(i in data$gen)
# {
#   A1[t,i]=1
#   t<-t+1
# }
#
# #Indexes of the environments
# A2<-matrix(rep(0, n*383), ncol=383)
# t<-1
# for(i in data$env)
# {
#   A2[t,i]=1
#   t<-t+1
# }
#
#
# #A2<-A2[,as.numeric(names(table(data$env)))]
# #Remove the missing environments
# A2<-A2[,-c(168,213,255,256,258,261,262,343)]
# VV<-data.frame(matrix(rep(0,n*22*p), nrow= n))
# t<-1
# for (i in data$gen)
# {
#   b<-1+(i-1)*p
#   e<-i*p
#   VV[t,b:e]<-V[t,]
#   t<-t+1
# }
#
# VV<-as.matrix(VV)
# p1<-dim(A1)[2]
# p2<-dim(A2)[2]
# p3<-dim(VV)[2]
# X<-cbind(A1[,-p1], A2[,-p2], VV)
#
#
# #A Model with only the main effects
# XN<-cbind(A1[,-p1], A2[,-p2])
# set.seed(1)
# simple.model<-lm(Y~XN)
# sum(is.na(coef(simple.model)))
# Yfit.simple<-predict(simple.model, newx = XN)
#
# #Ridge Regression
# set.seed(1)
# cv.model<-cv.glmnet(X, Y, alpha=0, family = "gaussian" )
# ridge.model<-glmnet(X, Y, alpha=0, family = "gaussian" , lambda = cv.model$lambda.min)
#
# #Lasso Regression
# set.seed(1)
# cv.model<-cv.glmnet(X, Y, alpha=1, family = "gaussian" )
# lasso.model<-glmnet(X, Y, alpha=1, family = "gaussian" , lambda = cv.model$lambda.min)
#
#
# #Lasso Regression (with main effects unpenalized)
# w<-rep(1, dim(X)[2])
# w[1:(p1+p2-2)]<-0
# set.seed(1)
# cv.model<-cv.glmnet(X, Y, alpha=1, family = "gaussian", penalty.factor=w  )
# lasso.model.unp<-glmnet(X, Y, alpha=1, family = "gaussian" ,
#                     lambda = cv.model$lambda.min, penalty.factor=w )
#
# #Adaptive Lasso with prespecified weights
# weights<-1/(abs(coef(lasso.model)[-1]))
# weights[1:(p1+p2-2)]<-rep(1,(p1+p2-2) )
# weights[is.infinite(weights)]<-1000
# set.seed(1)
# cv.model<-cv.glmnet(X, Y, alpha=1, family = "gaussian", penalty.factor=weights )
# weighted.lasso.model<-glmnet(X, Y, alpha=1, family = "gaussian" ,
#                     lambda = cv.model$lambda.min , penalty.factor=weights )
#
# #Lasso Regression with squared transformations
# set.seed(1)
# XBig<-cbind(A1[,-p1], A2[,-p2], VV, VV^2)
# cv.model<-cv.glmnet(XBig, Y, alpha=1, family = "gaussian" )
# lasso.big.model<-glmnet(XBig, Y, alpha=1, family = "gaussian" ,
#                         lambda = cv.model$lambda.min)
#
# #Post-Selection, with selected significant variables: Refitting
# #Select the non-zero coefficients of env. variables
# coeffic<-coef(lasso.model)[-1]
# coeffic<-coeffic[-c(1:(p2+p1-2))]
# #Keep the subsample of the significant
# positions<-(1:p3)[coeffic!=0]
# #Refitting with the selected variables, without penalty
# XPost<-cbind(A1[,-p1], A2[,-p2], VV[, positions])
# post.lasso.model<-glmnet(XPost, Y, lambda = 0, family = "gaussian")
#
#
# #Post-Selection, with selected significant variables
# #and Squared transfations
# XPost.big<-cbind(A1[,-p1], A2[,-p2], VV[, positions], VV[, positions]^2)
# post.lasso.big.model<-glmnet(XPost.big, Y, lambda = 0, family = "gaussian")
#
# #Fitted values
# Yfit.simple<-predict(simple.model, newx = XN)
# Yfit.ridge<-predict(ridge.model, newx = X)
# Yfit.lasso<-predict(lasso.model, newx = X)
# Yfit.lasso.unp<-predict(lasso.model.unp, newx = X)
# Yfit.w.lasso<-predict(weighted.lasso.model, newx = X)
# Yfit.big.lasso<-predict(lasso.big.model, newx = XBig)
# Yfit.post.lasso<-predict(post.lasso.model, newx = XPost)
# Yfit.post.big.lasso<-predict(post.lasso.big.model, newx = XPost.big)
#
#
# Results<-matrix(rep(0, 7*4), ncol=4)
#
# #RMSE for overall data
# Results[1,1]<-sqrt(mean((Yfit.simple-Y)^2))
# Results[2,1]<-sqrt(mean((Yfit.ridge-Y)^2))
# Results[3,1]<-sqrt(mean((Yfit.lasso-Y)^2))
# Results[4,1]<-sqrt(mean((Yfit.w.lasso-Y)^2))
# Results[5,1]<-sqrt(mean((Yfit.big.lasso-Y)^2))
# Results[6,1]<-sqrt(mean((Yfit.post.lasso-Y)^2))
# Results[7,1]<-sqrt(mean((Yfit.post.big.lasso-Y)^2))
#
#
#
# #Average RMSE over the environments
# Accuracy.rmse<-matrix(rep(0,7*383), ncol = 7)
# t<-1
# for(i in seq(1,383, by=1))
# {
#   Accuracy.rmse[t,1]<-sqrt(mean((Y[data$env==i]- Yfit.simple[data$env==i])^2))
#   Accuracy.rmse[t,2]<-sqrt(mean((Y[data$env==i]- Yfit.ridge[data$env==i])^2))
#   Accuracy.rmse[t,3]<-sqrt(mean((Y[data$env==i]- Yfit.lasso[data$env==i])^2))
#   Accuracy.rmse[t,4]<-sqrt(mean((Y[data$env==i]- Yfit.w.lasso[data$env==i])^2))
#   Accuracy.rmse[t,5]<-sqrt(mean((Y[data$env==i]- Yfit.big.lasso[data$env==i])^2))
#   Accuracy.rmse[t,6]<-sqrt(mean((Y[data$env==i]- Yfit.post.lasso[data$env==i])^2))
#   Accuracy.rmse[t,7]<-sqrt(mean((Y[data$env==i]- Yfit.post.big.lasso[data$env==i])^2))
#   t<-t+1
# }
#
# Accuracy.rmse<-Accuracy.rmse[as.numeric(names(table(data$env))),]
# Results[,2]<-colMeans(Accuracy.rmse)
#
# #Average Pearson correlation over the environments
# Accuracy.cor<-matrix(rep(0,7*383), ncol = 7)
# t<-1
# for(i in seq(1,383, by=1))
# {
#   Accuracy.cor[t,1]<-cor(Y[data$env==i], Yfit.simple[data$env==i])
#   Accuracy.cor[t,2]<-cor(Y[data$env==i], Yfit.ridge[data$env==i])
#   Accuracy.cor[t,3]<-cor(Y[data$env==i], Yfit.lasso[data$env==i])
#   Accuracy.cor[t,4]<-cor(Y[data$env==i], Yfit.w.lasso[data$env==i])
#   Accuracy.cor[t,5]<-cor(Y[data$env==i], Yfit.big.lasso[data$env==i])
#   Accuracy.cor[t,6]<-cor(Y[data$env==i], Yfit.post.lasso[data$env==i])
#   Accuracy.cor[t,7]<-cor(Y[data$env==i], Yfit.post.big.lasso[data$env==i])
#   t<-t+1
# }
#
# Accuracy.cor<-Accuracy.cor[as.numeric(names(table(data$env))),]
# Results[,3]<-colMeans(Accuracy.cor)
#
#
#
# Accuracy.cor.sp<-matrix(rep(0,7*383), ncol = 7)
# t<-1
# for(i in seq(1,383, by=1))
# {
#   Accuracy.cor.sp[t,1]<-cor(Y[data$env==i], Yfit.simple[data$env==i], method="spearman")
#   Accuracy.cor.sp[t,2]<-cor(Y[data$env==i], Yfit.ridge[data$env==i], method="spearman")
#   Accuracy.cor.sp[t,3]<-cor(Y[data$env==i], Yfit.lasso[data$env==i], method="spearman")
#   Accuracy.cor.sp[t,4]<-cor(Y[data$env==i], Yfit.w.lasso[data$env==i], method="spearman")
#   Accuracy.cor.sp[t,5]<-cor(Y[data$env==i], Yfit.big.lasso[data$env==i], method="spearman")
#   Accuracy.cor.sp[t,6]<-cor(Y[data$env==i], Yfit.post.lasso[data$env==i], method="spearman")
#   Accuracy.cor.sp[t,7]<-cor(Y[data$env==i], Yfit.post.big.lasso[data$env==i], method="spearman")
#   t<-t+1
# }
#
# Accuracy.cor.sp<-Accuracy.cor.sp[as.numeric(names(table(data$env))),]
# Results[,4]<-colMeans(Accuracy.cor.sp)
#
# rownames(Results)<-c("Simple","Ridge", "Lasso",
#                      "W. Lasso", "Sq. Lasso",
#                      "Post-Lasso","Sq. Post-Lasso" )
# colnames(Results)<-c("RMSE", "Aver. RMSE",
#                      "Aver. P. Corr.","Aver. Sp. Corr.")
#
#
# # Selected coefficients based on the lasso regression
# Betas<-coef(lasso.model)[-(1:(p1+p2-1))]
# #Genotypes X Coefficients
# Betas<-matrix(Betas, nrow = 22)
# colnames(Betas)<-colnames(data)[-(1:3)]
# sel.betas<-colnames(Betas)[colSums(Betas)!=0]
# length(sel.betas)
#
# Beta.selected<-Betas[,colSums(Betas)!=0]
#
# GxC<-data.frame(matrix(NA, ncol=10, nrow=dim(Beta.selected)[2]))
# for (i in seq(1, dim(Beta.selected)[2], by=1))
# {
#   w<-(which(Beta.selected[,i]!=0))
#   GxC[i,1:length(w)]<-w
# }
# rownames(GxC)<-colnames(Beta.selected)
#
# #Beta.all<-Beta.selected
# #Beta.all[abs(Beta.selected)>0]=1
# #Betas.num<-Betas
#
# #Beta.selected<-sort(colMeans(Betas.num)[colSums(Betas.num)!=0])
# #write.csv(data.frame(GxC), "Betas.csv")
#
#
#
#
# Acc<-Accuracy.cor.sp[,c(1,3,6,7)]
# colnames(Acc)<-c("Simple", "Lasso", "PLasso", "PLasso sq.")
# heatmap(t(Acc), Colv=NA, Rowv=NA, scale='none',
#         main="Spearman Cor")
#
# #
# # dd<-t(Accuracy.rmse[,c(1,3,6,7)])
# # ddd<-data.frame(dd)
# # barplot(as.matrix(ddd[,1:15]), beside = TRUE,
# #         xaxt="n", ylim=c(0,max(ddd[,1:15])),
# #         col=c(1,2,3,4))
# # axis(1, at=seq(1,75,by=5), labels= 1:15, las=3)
#
