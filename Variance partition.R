##################################
load('imac data.Rdata') # Data used are probit transformed. Ranking of genes are scaled by 1/(1 + n.gene) not 1/n.gene.

n.sample <- ncol(raw.imac)

n.gene <- nrow(raw.imac)

celltype <- annotation.imac$celltype

batch <- annotation.imac$Platform_Category


################ sample weight estimation (optional, I didn't use this in your paper, delete the 'weights' arguement in 'lmer' if you don't want it)##################
require(mgcv)

mu <- rowMeans(raw.imac)

sigma <- c()

for(i in 1:n.gene){
  
  sigma[i] <- summary(lm(raw.imac[i,] ~ celltype * batch))$sigma  
  
  print(i)
  
}


mod <- gam(sigma^0.5 ~ s(mu,bs="cr"))  # gene wise mean-variance trend estimation

weight <- matrix(ncol = n.sample, nrow = n.gene)

for(i in 1:n.sample){
  
  weight[,i] <- predict(mod, data.frame(mu = raw.imac[,i]))
  
  print(i)
  
}




################ fit mixture model ##################
require(lme4)
require(variancePartition)

varProp <- matrix(ncol =3, nrow= n.gene)

colnames(varProp) <- c('Batch', 'Celltype', 'Residuals')

rownames(varProp) <- rownames(raw.imac)
  
for(i in 1:n.gene){
  
  mixed.mod <-lmer(raw.imac[i,] ~ (1|celltype) + (1|batch),control=lmerControl(optimizer="bobyqa",
                                                               optCtrl=list(maxfun=2e5)), weights = weight[i,])
  
  # bad fit are not included
  
  if(!isSingular(mixed.mod)){ 
    
    varProp[i,] <- calcVarPart(mixed.mod)
    
  }
  
  print(i)

}



varRatio <- varProp[,'Celltype']/varProp[,'Batch'] # Calculate variance ratio

varRatio[is.na(varRatio)] <- 0 

varRatio <- sort(varRatio, decreasing = T)  # Sort genes according to variance ratio

varProp <- varProp[names(varRatio),]  

ribbonBound1 <- varProp[,'Residuals'] # Ribbon plot index

ribbonBound2 <- ribbonBound1 + varProp[,'Batch']






################ Plot ##################


MyColour <- c("#a0a3fa", "#F15D22", "#7C51A1") 
names(MyColour) <- c("Celltype", "Batch", "Residuals")


ggplot() + geom_ribbon(aes(x= 1:n.gene, ymin = 0, ymax=ribbonBound1, fill = 'Residuals')) + 
           geom_ribbon(aes(x= 1:n.gene, ymin = ribbonBound1, ymax=ribbonBound2, fill = 'Batch')) + 
           geom_ribbon(aes(x= 1:n.gene, ymin = ribbonBound2, ymax= 1, fill = 'Celltype')) +
           geom_vline(xintercept = which(names(varRatio) %in% rownames(imac)), alpha = 0.02) +
           scale_x_continuous(limits = c(1,sum(varRatio != 0)), expand = c(0,0)) + 
           scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
           scale_fill_manual('', values = MyColour) + 
           xlab('Genes ordered according to Progenitor Type / Platform variance ratio')+
           ylab('Proportion of variance explained')



# varProp is the requested table
write.table(varProp, 'variance parition result.txt')
