netopen<-function(TG_list,TF_list,enhancerOpen_trn,GExp_trn,GExp_tst)
{
  # this procedure is to predict enhancer's openness value from their neighboring genes' expression
  # TG_list is the list for enhancer's neighboring gene 
  # TF_list is the list for TFs that are binded in enhancer region
  # enhancerOpen_trn is enhancer openness value in training set (ENCODE 167 cells)
  # GExp_trn is gene expression in training set (ENCODE 167 cells)
  # GExp_tst is the expression in testing set
  
  # only enhancers with more two neighboring gene and binding TFs are remaining for training
  Flist<-lapply(1:length(enhancer.withOpen.TF.list),function(x) intersect(rownames(RNAdata.train),intersect(c(enhancer.withOpen.TF.list[[x]],enhancer.withOpen.TG.list[[x]]),rownames(RNAdata.test))))
  len<-sapply(1:length(Flist),function(x) length(Flist[[x]]))
  use_lab<-which(len>=2)
  enhancerOpen_trn1<-enhancerOpen_trn[use_lab,]
  Flist1<-Flist[use_lab]
  library(glmnet)
  prdEopen<-NULL
  for (k in 1:nrow(enhancerOpen_trn1)) {
    Xrn=t(GExp_trn[Flist1[[k]],])
    Yrn=enhancerOpen_trn1[k,]
    Xst=t(GExp_tst[Flist1[[k]],])
    cv.Ri <- cv.glmnet(Xrn, Yrn, alpha = 0)
    model.Ri <- glmnet(Xrn, Yrn, alpha = 0, lambda = cv.Ri$lambda.min)
    pre<-predict(model.Ri,newx = Xst)
    prdEopen<-cbind(prdEopen,pre)
    rm(Xrn,Yrn,Xst,cv.Ri,model.Ri,pre)
    cat(k,"\n")
  }
  return(netopen=t(prdEopen))
}
