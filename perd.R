perd<-function(TG_list,TF_list,enhancerOpen_trn,GExp_trn,drugGExp,drug_instance_info)
{
  # this procedure is to predict drug responsive regulatory elements (enhancers) by NetOpen
  # TG_list is the list for enhancer's neighboring gene 
  # TF_list is the list for TFs that are binded in enhancer region
  # enhancerOpen_trn is enhancer openness value in training set (ENCODE 167 cells)
  # GExp_trn is gene expression in training set (ENCODE 167 cells)
  # drugGExp is the expression under drug different treatment condition
  # drug_instance_info contains the drug treatment conditions
  
  # only enhancers with more two neighboring gene and binding TFs are remaining for training
  
  library(glmnet)
  Flist<-lapply(1:length(TG_list),function(x) intersect(rownames(GExp_trn),intersect(c(TG_list[[x]],TF_list[[x]]),rownames(drugGExp[[1]]))))
  len<-sapply(1:length(Flist),function(x) length(Flist[[x]]))
  use_lab<-which(len>=2)
  enhancerOpen_trn1<-enhancerOpen_trn[use_lab,]
  Flist1<-Flist[use_lab]
  predictedEnhancer.openness<-list()
  for (k in 1:length(drugGExp)) {
    A<-drugGExp[[k]]
    prediction<-NULL
    for (i in 1:nrow(enhancerOpen_trn1)) {
      Xrn<-t(GExp_trn[Flist1[[i]],])
      Yrn<-enhancerOpen_trn1[i,]
      Xst<-t(A[Flist1[[i]],])
      cv.Ri <- cv.glmnet(Xrn, Yrn, alpha = 0)
      model.Ri <- glmnet(Xrn, Yrn, alpha = 0, lambda = cv.Ri$lambda.min)
      pre<-predict(model.Ri,newx = Xst)
      prediction<-cbind(prediction,pre)
      rm(Xrn,Yrn,Xst,cv.Ri,model.Ri,pre)
      # cat(i,"\n")
    }
    colnames(prediction)<-rownames(enhancerOpen_trn1)
    predictedEnhancer.openness[[k]]<-prediction
    rm(A)
    cat(k,"\n")
  }
  names(predictedEnhancer.openness)<-names(drugGExp)
  
  # identify the enhancers with significant openness after given drug treatment
  library(limma)
  
  drug_diffEopen<-list()
  for (i in 1:length(predictedEnhancer.openness)) {
    A<-predictedEnhancer.openness[[i]]
    B<-drug_instance_info[[i]]
    ind<-B$ind
    design <- model.matrix(~ 0 + factor(ind))
    rownames(design) <- rownames(A)
    colnames(design)<-c("before","after")
    fit<-lmFit(t(A),design)
    contrast.matrix <- makeContrasts(before-after,levels=design)
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    fit3 <- eBayes(fit2)
    Results <- topTable(fit3,number=nrow(fit3))
    drug_diffEopen[[i]]<-Results
    rm(A,B,ind,design,fit,contrast.matrix,fit2,fit3,Results)
    cat(i,"\n")
  }
  names(drug_diffEopen)<-names(drugGExp)
  drug_diffEopen_list<-lapply(1:length(drug_diffEopen),function(x) {A<-drug_diffEopen[[x]];rownames(A)[which(abs(A$logFC)>0.5 & A$P.Value<0.001)]})
  names(drug_diffEopen_list)<-names(drugGExp)
  return(perd=drug_diffEopen_list)
}
