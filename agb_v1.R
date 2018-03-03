library(rpart)
library(tidyverse)
library(reshape2)

agb <- function(formula=formula(data),data,nesterov=TRUE,train.fraction=0.75,n.trees=100,shrinkage=0.01,new.data=data,distribution="gaussian",depth.tree=1,n.minobsinnode = 10){
  n <- nrow(data)
  set.seed(12345)
  perm <- sample(n);n.train <- round(n*train.fraction); n.valid <- n-n.train
  perm <- 1:n
  data2 <- model.frame(formula,data)
  train <- data2[perm[1:n.train],]
  valid <- data2[-perm[1:n.train],]
  new.data2 <- model.frame(formula,new.data)
  fitted <- matrix(0,nrow=n.train,ncol=n.trees)
  prev_valid <- matrix(0,nrow=n.valid,ncol=n.trees)
  prev_new <- matrix(0,nrow=nrow(new.data2),ncol=n.trees)
  g_fitted <- fitted
  g_prev_valid <- prev_valid
  g_prev_new <- prev_new
  err_train <- rep(0,n.trees)
  err_valid <- rep(0,n.trees)
  err_new <- rep(0,n.trees)
  if (distribution=="gaussian"){
    fitted[,1] <- mean(train[,1])
    prev_valid[,1] <- mean(train[,1])
    prev_new[,1] <- mean(train[,1])
    g_fitted[,1] <- fitted[,1]
    g_prev_valid[,1] <- prev_valid[,1]
    g_prev_new[,1] <- prev_new[,1]
    err_train[1] <- mean((train[,1]-fitted[,1])^2)
    err_valid[1] <- mean((valid[,1]-prev_valid[,1])^2)
    err_new[1] <- mean((new.data2[,1]-prev_new[,1])^2)
  }
  if (distribution=="adaboost"){
    fitted[,1] <- 0.5*log(sum(train[,1])/sum(1-train[,1]))
    prev_valid[,1] <- fitted[1,1]
    prev_new[,1] <- fitted[1,1]
    g_fitted[,1] <- fitted[,1]
    g_prev_valid[,1] <- prev_valid[,1]
    g_prev_new[,1] <- prev_new[,1]
    err_train[1] <- mean(exp(-(2*train[,1]-1)*fitted[,1]))
    err_valid[1] <- mean(exp(-(2*valid[,1]-1)*prev_valid[,1]))
    err_new[1] <- mean(exp(-(2*new.data2[,1]-1)*prev_new[,1]))
  }
  if (distribution=="bernoulli"){
    fitted[,1] <- log(mean(train[,1])/mean(1-train[,1]))
    prev_valid[,1] <- fitted[1,1]
    prev_new[,1] <- fitted[1,1]
    g_fitted[,1] <- fitted[,1]
    g_prev_valid[,1] <- prev_valid[,1]
    g_prev_new[,1] <- prev_new[,1]
    err_train[1] <- -2*mean(train[,1]*fitted[,1]-log(1+exp(fitted[,1])))
    err_valid[1] <- -2*mean(valid[,1]*prev_valid[,1]-log(1+exp(prev_valid[,1])))
    err_new[1] <- -2*mean(new.data2[,1]*prev_new[,1]-log(1+exp(prev_new[,1])))
    p_fitted <- fitted
    p_fitted[,1] <- 1/(1+exp(-fitted[,1]))
  }
  data_boucle <- data.frame(train[,-1])
  names(data_boucle) <- names(train)[-1]
  data_boucle$U <- 0
  tree.ctrl <- rpart.control(xval=1,maxdepth=depth.tree,minbucket=n.minobsinnode,cp=0.000001,maxcompete = 0,maxsurrogate = 0,usesurrogate = 0)
  #  tree.ctrl <- rpart.control(xval=1,maxdepth=depth.tree,minbucket=n.minobsinnode,minsplit=0,cp=0.000001,maxcompete = 0,maxsurrogate = 0,usesurrogate = 0)
  #  tree.ctrl <- rpart.control(xval=10,maxdepth=depth.tree,minsplit=2,cp=0.000001)
  lambda <- rep(1,n.trees)
  gamma <- rep(1,n.trees)
  for (i in 2:n.trees){
    lambda[i] <- 0.5*(1+sqrt(1+4*lambda[i-1]^2))
  }
  if (distribution=="gaussian"){
    for (i in 2:(n.trees)){
      if (nesterov){
        gamma[i] <- (1-lambda[i])/lambda[i+1]
        U <- train[,1]-g_fitted[,i-1]
        data_boucle$U <- U
        arbre <- rpart(U~.,data=data_boucle,control=tree.ctrl)
        data_boucle1 <- data_boucle
        data_boucle1$node.term <- arbre$where
        toto <- data_boucle1 %>% group_by(node.term) %>% summarize(prev.node=mean(U))
        arbre$frame[arbre$frame$var=="<leaf>","yval"] <- toto$prev.node
        fitted[,i] <- g_fitted[,i-1]+shrinkage*predict(arbre)
        prev_valid[,i] <- g_prev_valid[,i-1]+shrinkage*predict(arbre,newdata=valid)
        prev_new[,i] <- g_prev_new[,i-1]+shrinkage*predict(arbre,newdata=new.data2)
        g_fitted[,i] <- (1-gamma[i-1])*fitted[,i]+gamma[i-1]*fitted[,i-1]
        g_prev_valid[,i] <- (1-gamma[i-1])*prev_valid[,i]+gamma[i-1]*prev_valid[,i-1]
        g_prev_new[,i] <- (1-gamma[i-1])*prev_new[,i]+gamma[i-1]*prev_new[,i-1]
      } else {
        U <- train[,1]-fitted[,i-1]
        data_boucle$U <- U
        arbre <- rpart(U~.,data=data_boucle,control=tree.ctrl)
        data_boucle1 <- data_boucle
        data_boucle1$node.term <- arbre$where
        toto <- data_boucle1 %>% group_by(node.term) %>% summarize(prev.node=mean(U))
        arbre$frame[arbre$frame$var=="<leaf>","yval"] <- toto$prev.node
        fitted[,i] <- fitted[,i-1]+shrinkage*predict(arbre)
        prev_valid[,i] <- prev_valid[,i-1]+shrinkage*predict(arbre,newdata=valid)
        prev_new[,i] <- prev_new[,i-1]+shrinkage*predict(arbre,newdata=new.data2)
      }
      err_train[i] <- mean((train[,1]-fitted[,i])^2)
      err_valid[i] <- mean((valid[,1]-prev_valid[,i])^2)
      err_new[i] <- mean((new.data2[,1]-prev_new[,i])^2)
      if (err_valid[i]>=100){
        n.trees <- i
        err_train <- err_train[1:i]
        err_valid <- err_valid[1:i]
        err_new <- err_new[1:i]
        prev_new <- prev_new[,1:i]
        break
      }
    }
  }
  if (distribution=="adaboost"){
    for (i in 2:(n.trees)){
      if (nesterov){
        gamma[i] <- (1-lambda[i])/lambda[i+1]
        U <- -(2*train[,1]-1)*exp(-(2*train[,1]-1)*g_fitted[,i-1])
        U1 <- exp(-(2*train[,1]-1)*g_fitted[,i-1])
        data_boucle$U <- U
        arbre <- rpart(U~.,data=data_boucle,control=tree.ctrl)
        #      arbre <- prune(arbre,arbre$cptable[which(arbre$cptable[,"nsplit"]==2),"CP"])
        data_boucle1 <- data_boucle
        data_boucle1$U1 <- U1
        data_boucle1$node.term <- arbre$where
        toto <- data_boucle1 %>% group_by(node.term) %>% summarize(prev.node=mean(-U)/mean(U1))
        arbre$frame[arbre$frame$var=="<leaf>","yval"] <- toto$prev.node
        fitted[,i] <- g_fitted[,i-1]+shrinkage*predict(arbre)
        prev_valid[,i] <- g_prev_valid[,i-1]+shrinkage*predict(arbre,newdata=valid)
        prev_new[,i] <- g_prev_new[,i-1]+shrinkage*predict(arbre,newdata=new.data2)
        g_fitted[,i] <- (1-gamma[i-1])*fitted[,i]+gamma[i-1]*fitted[,i-1]
        g_prev_valid[,i] <- (1-gamma[i-1])*prev_valid[,i]+gamma[i-1]*prev_valid[,i-1]
        g_prev_new[,i] <- (1-gamma[i-1])*prev_new[,i]+gamma[i-1]*prev_new[,i-1]
      } else {
        U <- -(2*train[,1]-1)*exp(-(2*train[,1]-1)*fitted[,i-1])
        U1 <- exp(-(2*train[,1]-1)*fitted[,i-1])
        data_boucle$U <- U
        arbre <- rpart(U~.,data=data_boucle,control=tree.ctrl)
        #      arbre <- prune(arbre,arbre$cptable[which(arbre$cptable[,"nsplit"]==2),"CP"])
        data_boucle1 <- data_boucle
        data_boucle1$U1 <- U1
        data_boucle1$node.term <- arbre$where
        toto <- data_boucle1 %>% group_by(node.term) %>% summarize(prev.node=mean(-U)/mean(U1))
        arbre$frame[arbre$frame$var=="<leaf>","yval"] <- toto$prev.node
        fitted[,i] <- fitted[,i-1]+shrinkage*predict(arbre)
        prev_valid[,i] <- prev_valid[,i-1]+shrinkage*predict(arbre,newdata=valid)
        prev_new[,i] <- prev_new[,i-1]+shrinkage*predict(arbre,newdata=new.data2)
      }
      err_train[i] <- mean(exp(-(2*train[,1]-1)*fitted[,i]))
      err_valid[i] <- mean(exp(-(2*valid[,1]-1)*prev_valid[,i]))
      err_new[i] <- mean(exp(-(2*new.data2[,1]-1)*prev_new[,i]))
      if (err_valid[i]>=100){
        n.trees <- i
        err_train <- err_train[1:i]
        err_valid <- err_valid[1:i]
        err_new <- err_new[1:i]
        prev_new <- prev_new[,1:i]
        break
      }
    }
  }
  if (distribution=="bernoulli"){
    for (i in 2:(n.trees)){
      if (nesterov){
        gamma[i] <- (1-lambda[i])/lambda[i+1]
        U <- train[,1]-1/(1+exp(-g_fitted[,i-1]))
        U1 <- p_fitted[,i-1]*(1-p_fitted[,i-1])
        U2 <- train[,1]-p_fitted[,i-1]
        data_boucle$U <- U
        arbre <- rpart(U~.,data=data_boucle,control=tree.ctrl)
        data_boucle1 <- data_boucle
        data_boucle1$U1 <- U1
        data_boucle1$U2 <- U2
        data_boucle1$node.term <- arbre$where
        toto <- data_boucle1 %>% group_by(node.term) %>% summarize(prev.node=mean(U2)/mean(U1))
        arbre$frame[arbre$frame$var=="<leaf>","yval"] <- toto$prev.node
        fitted[,i] <- g_fitted[,i-1]+shrinkage*predict(arbre)
        prev_valid[,i] <- g_prev_valid[,i-1]+shrinkage*predict(arbre,newdata=valid)
        prev_new[,i] <- g_prev_new[,i-1]+shrinkage*predict(arbre,newdata=new.data2)
        g_fitted[,i] <- (1-gamma[i-1])*fitted[,i]+gamma[i-1]*fitted[,i-1]
        p_fitted[,i] <- 1/(1+exp(-g_fitted[,i]))
        g_prev_valid[,i] <- (1-gamma[i-1])*prev_valid[,i]+gamma[i-1]*prev_valid[,i-1]
        g_prev_new[,i] <- (1-gamma[i-1])*prev_new[,i]+gamma[i-1]*prev_new[,i-1]
      } else {
        U <- train[,1]-1/(1+exp(-fitted[,i-1]))
        U1 <- p_fitted[,i-1]*(1-p_fitted[,i-1])
        U2 <- train[,1]-p_fitted[,i-1]
        data_boucle$U <- U
        arbre <- rpart(U~.,data=data_boucle,control=tree.ctrl)
        data_boucle1 <- data_boucle
        data_boucle1$U1 <- U1
        data_boucle1$U2 <- U2
        data_boucle1$node.term <- arbre$where
        toto <- data_boucle1 %>% group_by(node.term) %>% summarize(prev.node=mean(U2)/mean(U1))
        arbre$frame[arbre$frame$var=="<leaf>","yval"] <- toto$prev.node
        fitted[,i] <- fitted[,i-1]+shrinkage*predict(arbre)
        p_fitted[,i] <- 1/(1+exp(-fitted[,i]))
        prev_valid[,i] <- prev_valid[,i-1]+shrinkage*predict(arbre,newdata=valid)
        prev_new[,i] <- prev_new[,i-1]+shrinkage*predict(arbre,newdata=new.data2)
      }
      err_train[i] <- -2*mean(train[,1]*fitted[,i]-log(1+exp(fitted[,i])))
      err_valid[i] <- -2*mean(valid[,1]*prev_valid[,i]-log(1+exp(prev_valid[,i])))
      err_new[i] <- -2*mean(new.data2[,1]*prev_new[,i]-log(1+exp(prev_new[,i])))
      if (err_valid[i]>=100){
        n.trees <- i
        err_train <- err_train[1:i]
        err_valid <- err_valid[1:i]
        err_new <- err_new[1:i]
        prev_new <- prev_new[,1:i]
        break
      }
    }
  }
  
  df <- data.frame(x=1:n.trees,train=err_train,valid=err_valid)
  df1 <- melt(df,id.vars = "x")
  names(df1) <- c("Iterations","Dataset","Error")
  graph_err <- ggplot(df1)+aes(x=Iterations,y=Error,color=Dataset)+geom_line(size=1)+geom_vline(xintercept = which.min(df$valid),color="blue",size=1)+theme_classic()
  return(list(graph=graph_err,error=df1,Mopt=which.min(err_valid),prev_new=prev_new,err_new=err_new))
}
