
ngenes <- 2000
baselineprop <- qAbundanceDist( (1:ngenes)/(ngenes+1) )
baselineprop <- baselineprop/sum(baselineprop)



ngroup <- 2
nlibs <- 6
group <- rep(1:2,each=nlibs/2)
nbatch <- 4
p <- runif(4)
batch <- sample(1:nbatch, nlibs, p =p/sum(p), replace = T)
expected.lib.size <- 11e6* runif(nlibs,min = 0.1, max = 1)


limma.df <- data.frame(matrix(nrow = 100, ncol = 3))
colnames(limma.df) <- c('type1','type2', 'fc')




for(j in 1:100){

ndiff <- c(100,100)#round(runif(ngroup, ngenes/200, ngenes/50)) #how many genes are differentially expressed in each  group
dif.express <- list()
fc <- list()
for(i in 1:ngroup){
  dif.express[[i]] <- sample(1:ngenes, ndiff[i]) #sample which genes are differentially expressed in each group
  fc[[i]] <- runif(ndiff[i], 1, 3)
}


prop.matrix <- matrix(rep(baselineprop,nlibs),ncol = nlibs) #set up the baseline proportion matrix
for(i in 1:ngroup){
  prop.matrix[dif.express[[i]],group == i] <- prop.matrix[dif.express[[i]],group == i] * fc[[i]]
}



mu0 <- t(t(prop.matrix) * expected.lib.size)



BCV0 <- 0.2+1/sqrt(mu0) #trend of dispersion
df.BCV <- 40
BCV <- BCV0*sqrt(df.BCV/rchisq(ngenes,df=df.BCV))
shape <- 1/BCV^2
scale <- mu0/shape
mu <- matrix(rgamma(ngenes*nlibs,shape=shape,scale=scale),ngenes,nlibs)
counts <- matrix(rpois(ngenes*nlibs,lambda=mu),ngenes,nlibs)


design <- model.matrix(~0 + as.factor(group)) #design matrix 
colnames(design) <- c('cellA','cellB') #rename design matrix
contrastsMatrix <- makeContrasts(c1vsc2 = cellA - cellB, levels = colnames(design)) #contrast matrix between cellA and cellB, test H0: c1 = c2
dgevoom <- voom(counts, design, plot = T) #weights the residual
dgeFit <- lmFit(dgevoom, design) #Fit the model with weighted least square
dgeFit <- contrasts.fit(dgeFit, contrasts = contrastsMatrix) #Estimate contrast with fitted model
dgeEFit <- eBayes(dgeFit) 
dgeResults <- topTable(dgeEFit, coef=1, n=Inf)
dgeResults$error <- NA
sum(dgeResults$adj.P.Val <= 0.05)

limma.difgene  <- rownames(dgeResults)[dgeResults$adj.P.Val <= 0.05]  #recovered differntially expressed gene
true.difgene <- as.character(do.call(c,dif.express))
type2.gene <- true.difgene[!(true.difgene %in% limma.difgene)] #true differntial expressed gene that are not found by limma (type 2 error)
type1.gene <- limma.difgene[!(limma.difgene %in% true.difgene)] #limma recovered gene that are not truly differntially expressed (type 1 error)


udiscovered.fc <- c(fc[[1]][which(dif.express[[1]] %in% type2.gene)], fc[[2]][which(dif.express[[2]] %in% type2.gene)])
limma.df[j, 'type1'] <- length(type1.gene)
limma.df[j, 'type2'] <- length(type2.gene)
limma.df[j, 'fc'] <- mean(udiscovered.fc)
}




require(mgcv)
qnorm.df <- data.frame(matrix(nrow = 100, ncol = 3))
colnames(qnorm.df) <- c('type1','type2', 'fc')

for(j in 1:100){
  
  ndiff <- c(100,100)#round(runif(ngroup, ngenes/200, ngenes/50)) #how many genes are differentially expressed in each  group
  dif.express <- list()
  fc <- list()
  for(i in 1:ngroup){
    dif.express[[i]] <- sample(1:ngenes, ndiff[i]) #sample which genes are differentially expressed in each group
    fc[[i]] <- runif(ndiff[i], 1, 3)
  }
  
  
  prop.matrix <- matrix(rep(baselineprop,nlibs),ncol = nlibs) #set up the baseline proportion matrix
  for(i in 1:ngroup){
    prop.matrix[dif.express[[i]],group == i] <- prop.matrix[dif.express[[i]],group == i] * fc[[i]]
  }
  
  
  
  mu0 <- t(t(prop.matrix) * expected.lib.size)
  
  
  
  BCV0 <- 0.2+1/sqrt(mu0) #trend of dispersion
  df.BCV <- 40
  BCV <- BCV0*sqrt(df.BCV/rchisq(ngenes,df=df.BCV))
  shape <- 1/BCV^2
  scale <- mu0/shape
  mu <- matrix(rgamma(ngenes*nlibs,shape=shape,scale=scale),ngenes,nlibs)
  counts <- matrix(rpois(ngenes*nlibs,lambda=mu),ngenes,nlibs)
  rank.count <- qnorm(apply(t(counts),1,rank)/(ngenes + 1))


  weights <- data.frame(matrix(nrow = ngenes*2,ncol = 2))
  colnames(weights) <- c('mean','sd')
  weights[1:ngenes,'mean'] <- apply(rank.count[,1:3],1,mean)
  weights[1:ngenes,'sd'] <- apply(rank.count[,1:3],1,sd)
  weights[(ngenes+1):(ngenes*2),'mean'] <- apply(rank.count[,4:6],1,mean)
  weights[(ngenes+1):(ngenes*2),'sd'] <- apply(rank.count[,4:6],1,sd)
  weights <- arrange(weights, mean)

  
  gam<- gam(sd ~ s(mean),data = weights)
  weights<- 1/matrix(predict(gam, data.frame(mean = as.numeric(rank.count))), ncol = nlibs)^4
  plot(gam)
  
  
  fit <- lmFit(rank.count, design, weights = weights)
  fit <- contrasts.fit(fit, contrasts = contrastsMatrix)
  fit <- topTable(eBayes(fit), coef=1, n=Inf)
  rownames(fit)[fit$adj.P.Val < 0.05]
  
  
  
  qnorm.difgene  <- rownames(fit)[fit$adj.P.Val < 0.05]  #recovered differntially expressed gene
  true.difgene <- as.character(do.call(c,dif.express))
  type2.gene <- true.difgene[!(true.difgene %in% qnorm.difgene)] #true differntial expressed gene that are not found by limma (type 2 error)
  type1.gene <- qnorm.difgene[!(qnorm.difgene %in% true.difgene)] #limma recovered gene that are not truly differntially expressed (type 1 error)

  udiscovered.fc <- c(fc[[1]][which(dif.express[[1]] %in% type2.gene)], fc[[2]][which(dif.express[[2]] %in% type2.gene)])

  qnorm.df[j, 'type1'] <- length(type1.gene)
  qnorm.df[j, 'type2'] <- length(type2.gene)
  qnorm.df[j, 'fc'] <- mean(udiscovered.fc)
  print(j)
}



rank.df <- data.frame(matrix(nrow = 100, ncol = 3))
colnames(rank.df) <- c('type1','type2', 'fc')

for(j in 1:100){
  
  ndiff <- c(100,100)#round(runif(ngroup, ngenes/200, ngenes/50)) #how many genes are differentially expressed in each  group
  dif.express <- list()
  fc <- list()
  for(i in 1:ngroup){
    dif.express[[i]] <- sample(1:ngenes, ndiff[i]) #sample which genes are differentially expressed in each group
    fc[[i]] <- runif(ndiff[i], 1, 3)
  }
  
  
  prop.matrix <- matrix(rep(baselineprop,nlibs),ncol = nlibs) #set up the baseline proportion matrix
  for(i in 1:ngroup){
    prop.matrix[dif.express[[i]],group == i] <- prop.matrix[dif.express[[i]],group == i] * fc[[i]]
  }
  
  
  
  mu0 <- t(t(prop.matrix) * expected.lib.size)
  
  
  
  BCV0 <- 0.2+1/sqrt(mu0) #trend of dispersion
  df.BCV <- 40
  BCV <- BCV0*sqrt(df.BCV/rchisq(ngenes,df=df.BCV))
  shape <- 1/BCV^2
  scale <- mu0/shape
  mu <- matrix(rgamma(ngenes*nlibs,shape=shape,scale=scale),ngenes,nlibs)
  counts <- matrix(rpois(ngenes*nlibs,lambda=mu),ngenes,nlibs)
  rank.count <- apply(t(counts),1,rank)
  

  

  fit <- lmFit(rank.count, design)
  fit <- contrasts.fit(fit, contrasts = contrastsMatrix)
  fit <- topTable(eBayes(fit), coef=1, n=Inf)



  difgene  <- rownames(fit)[fit$adj.P.Val < 0.05]  #recovered differntially expressed gene
  true.difgene <- as.character(do.call(c,dif.express))
  type2.gene <- true.difgene[!(true.difgene %in% difgene)] #true differntial expressed gene that are not found by limma (type 2 error)
  type1.gene <- difgene[!(difgene %in% true.difgene)] #limma recovered gene that are not truly differntially expressed (type 1 error)
  
  udiscovered.fc <- c(fc[[1]][which(dif.express[[1]] %in% type2.gene)], fc[[2]][which(dif.express[[2]] %in% type2.gene)])
  
  rank.df[j, 'type1'] <- length(type1.gene)
  rank.df[j, 'type2'] <- length(type2.gene)
  rank.df[j, 'fc'] <- mean(udiscovered.fc)
  print(j)
}


noweight.df <- data.frame(matrix(nrow = 100, ncol = 3))
colnames(noweight.df) <- c('type1','type2', 'fc')

for(j in 1:100){
  
  ndiff <- c(100,100)#round(runif(ngroup, ngenes/200, ngenes/50)) #how many genes are differentially expressed in each  group
  dif.express <- list()
  fc <- list()
  for(i in 1:ngroup){
    dif.express[[i]] <- sample(1:ngenes, ndiff[i]) #sample which genes are differentially expressed in each group
    fc[[i]] <- runif(ndiff[i], 1, 3)
  }
  
  
  prop.matrix <- matrix(rep(baselineprop,nlibs),ncol = nlibs) #set up the baseline proportion matrix
  for(i in 1:ngroup){
    prop.matrix[dif.express[[i]],group == i] <- prop.matrix[dif.express[[i]],group == i] * fc[[i]]
  }
  
  
  
  mu0 <- t(t(prop.matrix) * expected.lib.size)
  
  
  
  BCV0 <- 0.2+1/sqrt(mu0) #trend of dispersion
  df.BCV <- 40
  BCV <- BCV0*sqrt(df.BCV/rchisq(ngenes,df=df.BCV))
  shape <- 1/BCV^2
  scale <- mu0/shape
  mu <- matrix(rgamma(ngenes*nlibs,shape=shape,scale=scale),ngenes,nlibs)
  counts <- matrix(rpois(ngenes*nlibs,lambda=mu),ngenes,nlibs)
  rank.count <- qnorm(apply(t(counts),1,rank)/(ngenes+1))
  
  


  fit <- lmFit(rank.count, design)
  fit <- contrasts.fit(fit, contrasts = contrastsMatrix)
  fit <- topTable(eBayes(fit), coef=1, n=Inf)
  
  
  
  difgene  <- rownames(fit)[fit$adj.P.Val < 0.05]  #recovered differntially expressed gene
  true.difgene <- as.character(do.call(c,dif.express))
  type2.gene <- true.difgene[!(true.difgene %in% difgene)] #true differntial expressed gene that are not found by limma (type 2 error)
  type1.gene <- difgene[!(difgene %in% true.difgene)] #limma recovered gene that are not truly differntially expressed (type 1 error)
  
  udiscovered.fc <- c(fc[[1]][which(dif.express[[1]] %in% type2.gene)], fc[[2]][which(dif.express[[2]] %in% type2.gene)])
  
  noweight.df[j, 'type1'] <- length(type1.gene)
  noweight.df[j, 'type2'] <- length(type2.gene)
  noweight.df[j, 'fc'] <- mean(udiscovered.fc)
  print(j)
}




colMeans(limma.df)
colMeans(qnorm.df)
colMeans(rank.df)
colMeans(noweight.df)



df <- rbind(limma.df, qnorm.df, rank.df, noweight.df)
df$mod <- rep(c('limma','probit', 'rank', 'no.weight.probit'),each = 100)
df$true.discov <- 200 - df$type2
df$FDR <- df$type1 / (df$true.discov + df$type1)


ggplot() +geom_boxplot(data = df,aes(x=mod,y = type1, fill = mod)) + theme_bw()
ggplot() +geom_boxplot(data = df,aes(x=mod,y = type2, fill = mod)) + theme_bw()
ggplot() +geom_boxplot(data = df,aes(x=mod,y = fc, fill = mod))+ theme_bw()
ggplot() + geom_boxplot(data = df,aes(x=mod,y = FDR, fill = mod)) + theme_bw() + geom_hline(yintercept = 0.05, linetype = 'dashed')

voom(counts, design, plot = T, save.plot = T)$voom.xy
hist(voom(counts, design, plot = T)$E)
hist(rank.count)

