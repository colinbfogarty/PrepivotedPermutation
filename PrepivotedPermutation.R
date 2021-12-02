#univariate, greater than alternative
PermAllOS<- function(data, ind, mvec=0, nrep = 1000)
{
  require(mvtnorm)
  require(boot)
  require(emulator)
  require(parallel)
  require(CompQuadForm)
  require(multcomp)
  DiMStuOS <- function(data, ind, mvec=0, perm = F)
  {
    require(mvtnorm)
    require(boot)
    require(emulator)
    require(parallel)
    require(CompQuadForm)
    y = data[ind,2]
    strat = data[ind,1]
    if(perm == T)
    {
      strat = data[,1]
    }
    y1 = y[strat==1]
    y2 = y[strat==2]
    n1 = length(y1)
    n2 = length(y2)
    n = n1+n2
    stu = (mean(y1) - mean(y2) - mvec)/sqrt(var(y1)/length(y1) + var(y2)/length(y2))
    sighat = sqrt(var(y1)/n1 + var(y2)/n2)*sqrt(n)
    mu3hat = ((n/n1)^3*sum((y1 - mean(y1))^3) - (n/n2)^3*sum((y2 - mean(y2))^3))/n
    edge = pnorm(stu) + ((mu3hat)/(6*sqrt(n)*sighat^3))*dnorm(stu)*(2*stu^2+1)
    c((mean(y1) - mean(y2) - mvec), (mean(y1) - mean(y2) - mvec)/sqrt(var(y1)/length(y1) + var(y2)/length(y2)), edge)
  }

  y = data[ind,2]
  strat = data[,1]
  y1 = y[strat==1]
  y2 = y[strat==2]
  n1 = length(y1)
  n2 = length(y2)
  n = n1+n2
  dim = (mean(y1) - mean(y2) - mvec)
  stu = (mean(y1) - mean(y2) - mvec)/sqrt(var(y1)/n1 + var(y2)/n2)
  sighat = sqrt(var(y1)/n1 + var(y2)/n2)*sqrt(n)
  mu3hat = ((n/n1)^3*sum((y1 - mean(y1))^3) - (n/n2)^3*sum((y2 - mean(y2))^3))/n
  edge = pnorm(stu) + ((mu3hat)/(6*sqrt(n)*sighat^3))*dnorm(stu)*(2*stu^2+1)
  datnew = cbind(strat, y)
  predim = 0
  preedge = 0
  prestu = 0
  if(nrep > 0)
  {
    aall = boot(data = datnew, statistic = DiMStuOS, R = nrep, strata = datnew[,1],  mvec = mean(y1)-mean(y2))$t
    adim = aall[,1]
    astu = aall[,2]
    aedge = aall[,3]
    predim = mean(adim <= dim)
    prestu = mean(astu <= stu)
    preedge = mean(aedge <= edge)
  }
  c(dim = dim, stu = stu, edge = edge, predim = predim, prestu = prestu, preedge = preedge)
}



#multivariate
PermAllMult <- function(data, ind, mvec=0, nrep = 1000)
{
  require(mvtnorm)
  require(boot)
  require(emulator)
  require(parallel)
  require(CompQuadForm)
  require(multcomp)
  MaxPre <- function(data, ind, mvec=0, perm = F)
  {
    require(mvtnorm)
    require(boot)
    require(emulator)
    require(parallel)
    require(CompQuadForm)
    require(multcomp)
    y = data[ind,-1]
    strat = data[ind,1]
    if(perm == T)
    {
      strat = data[,1]
    }
    y1 = y[strat==1,]
    y2 = y[strat==2,]
    n1 = nrow(y1)
    n2 = nrow(y2)
    cov1 = (cov(y1)/n1 + cov(y2)/n2)
    covpool = ((cov(y1)*(n1-1)+ cov(y2)*(n2-1))/(n1+n2-2))*(1/n1+1/n2)
    
    ses = sqrt(diag(cov1))
    covstu = cov2cor(cov1)
    dimm = (colMeans(y1) - colMeans(y2) - mvec)
    stu = abs(dimm)/ses
    mm = max(stu)
    hot <- quad.form.inv(cov1, dimm)
    pool <- quad.form.inv(covpool, dimm)
    datnew = cbind(strat, y)
    #this would use pmvnorm to prepivot, but it is computationally intensive
    #we replace it with 1000 monte-carlo simulations from the relevant Gaussian
    #premm  = as.numeric(pmvnorm(lower = -rep(mm, ncol(y1)), upper = rep(mm, ncol(y1)), sigma = covstu))
    BB = rmvnorm(1000, sigma = covstu)
    premm = mean(apply(abs(BB), 1, max) <= mm)
    vchalf =(t(chol(cov1)))
    A1 = quad.form.inv(covpool, vchalf)
    eig = eigen(A1)
    lambda = eig$values
    prepool= as.numeric((1 - davies(pool, lambda, h = rep(1, length(lambda)), delta = rep(0, nrow(cov1)), lim = 100000, acc = .000005)$Qq))
    c(hot, pool, mm, prepool, premm)
  }
  y = data[ind,-1]
  strat = data[,1]
  y1 = y[strat==1,]
  y2 = y[strat==2,]
  n1 = nrow(y1)
  n2 = nrow(y2)
  cov1 = (cov(y1)/n1 + cov(y2)/n2)
  covpool = ((cov(y1)*(n1-1)+ cov(y2)*(n2-1))/(n1+n2-2))*(1/n1+1/n2)
  
  ses = sqrt(diag(cov1))
  covstu = cov2cor(cov1)
  dimm = (colMeans(y1) - colMeans(y2) - mvec)
  stu = abs(dimm)/ses
  mm = max(stu)
  hot <- quad.form.inv(cov1, dimm)
  pool <- quad.form.inv(covpool, dimm)
  #this would use pmvnorm to prepivot, but it is computationally intensive
  #we replace it with 1000 monte-carlo simulations from the relevant Gaussian
  #premm = as.numeric(pmvnorm(lower = rep(-mm, ncol(y1)), upper = rep(mm, ncol(y1)), sigma = covstu))
  BB = rmvnorm(1000, sigma = covstu)
  premm = mean(apply(abs(BB), 1, max) <= mm)
  vchalf =(t(chol(cov1)))
  A1 = quad.form.inv(covpool, vchalf)
  eig = eigen(A1)
  lambda = eig$values
  prepool= as.numeric((1 - davies(pool, lambda, h = rep(1, length(lambda)), delta = rep(0, nrow(cov1)), lim = 100000, acc = .000005)$Qq))
  
  datnew = cbind(strat, y)
  boothot = 0
  bootpool = 0
  bootmm = 0
  bootprepool = 0
  bootpremm = 0
  if(nrep > 0)
  {
    aall = boot(data = datnew, statistic = MaxPre, R = nrep, strata = datnew[,1], mvec = colMeans(y1)-colMeans(y2))$t
    ahot= aall[,1]
    apool = aall[,2]
    amm = aall[,3]
    aprepool = aall[,4]
    apremm = aall[,5]
    boothot = mean(ahot <= hot)
    bootpool = mean(apool <= pool)
    bootmm = mean(amm <= mm)
    bootprepool = mean(aprepool <= prepool)
    bootpremm = mean(apremm <= premm)
  }
  
  
  c(hot = hot, pool = pool, mm = mm, prepool = prepool, premm = premm, boothot = boothot, bootpool = bootpool, bootmm = bootmm, bootprepool = bootprepool, bootpremm = bootpremm)
}




#ANOVA
PermAnova <- function(data, ind, mvec=rep(0, length(unique(data[,1]))), C, nrep = 1000)
{
  require(mvtnorm)
  require(boot)
  require(emulator)
  require(parallel)
  require(CompQuadForm)
  require(multcomp)
  AnovaPre <- function(data, ind, mvec=rep(0, length(unique(data[,1]))), C, perm = F)
  {
    require(mvtnorm)
    require(boot)
    require(emulator)
    require(parallel)
    require(CompQuadForm)
    require(multcomp)
    y = data[ind,-1]
    strat = data[ind,1]
    if(perm == T)
    {
      strat = data[,1]
    }
    y = y - mvec[strat]
    k = length(unique(strat))
    n = table(strat)
    N = sum(n)
    SSE = sum((n-1)*tapply(y, strat, var))
    MSE = SSE/(N-k)
    muhat = tapply(y, strat, mean)
    Gammahat = diag(MSE/n)
    SSR = sum(n*(tapply(y, strat, mean) - mean(y))^2)
    Vhat = quad.form(Gammahat, C)
    mhat = muhat%*%C
    tukey = max(abs(mhat)/diag(Vhat)^(1/2))
    fstat = SSR/(k-1)/MSE
    sigma2 = tapply(y, strat, var)
    CR = sum((n/sigma2)*(muhat - sum((n/sigma2)*muhat)/sum((n/sigma2)))^2)
    c(fstat, tukey, CR)
  }
  y = data[ind,-1]
  strat = data[,1]
  k = length(unique(strat))
  n = table(strat)
  N = sum(n)
  SSE = sum((n-1)*tapply(y, strat, var))
  MSE = SSE/(N-k)
  muhat = tapply(y, strat, mean)
  Gammahat = diag(MSE/n)
  SSR = sum(n*(tapply(y, strat, mean) - mean(y))^2)
  Vhat = quad.form(Gammahat, C)
  mhat = muhat%*%C
  tukey = max(abs(mhat)/diag(Vhat)^(1/2))
  fstat = SSR/(k-1)/MSE
  sigma2 = tapply(y, strat, var)
  CR = sum((n/sigma2)*(muhat - sum((n/sigma2)*muhat)/sum((n/sigma2)))^2)
  mvec = muhat
  bootfstat = 0
  boottukey = 0
  bootCR = 0
  datnew = cbind(strat, y)
  if(nrep > 0)
  {
    aall = boot(data = datnew, statistic = AnovaPre, R = nrep, strata = datnew[,1], mvec = mvec, C = C)$t
    afstat= aall[,1]
    atukey = aall[,2]
    aCR = aall[,3]
    bootfstat = mean(afstat<= fstat)
    boottukey  = mean(atukey <= tukey)
    bootCR = mean(aCR <= CR)
  }
  c(fstat=fstat, tukey = tukey,CR = CR, bootfstat = bootfstat, boottukey = boottukey,  bootCR = bootCR)
}



#MANOVA
PermManova <- function(data, ind, mvec=matrix(0, length(unique(data[,1]), ncol(data[,-1]))), nrep = 1000)
{
  require(mvtnorm)
  require(boot)
  require(emulator)
  require(parallel)
  require(CompQuadForm)
  require(multcomp)
  ManovaPre <- function(data, ind, mvec=matrix(0, length(unique(data[,1]), ncol(data[,-1]))), perm = F)
  {
    require(mvtnorm)
    require(boot)
    require(emulator)
    require(parallel)
    require(CompQuadForm)
    require(multcomp)
    y = data[ind,-1]
    strat = data[ind,1]
    if(perm == T)
    {
      strat = data[,1]
    }
    y = y - mvec[strat,]
    fstrat = factor(strat)
    linmod = lm(y~fstrat)
    mm = manova(y~fstrat)
    eigs = summary(mm)$Eigenvalues
    tP = sum(eigs/(1+eigs))
    tL = sum(eigs)
    tW = -prod(1/(1+eigs))
    tR = max(eigs)
    c(tP = tP, tL = tL, tW = tW, tR = tR)
  }
  y = data[ind,-1]
  strat = data[,1]
  fstrat = factor(strat)
  linmod = lm(y~fstrat)
  mm = manova(y~fstrat)
  eigs = summary(mm)$Eigenvalues
  tP = sum(eigs/(1+eigs))
  tL = sum(eigs)
  tW = -prod(1/(1+eigs))
  tR = max(eigs)
  
  dd = data.frame(fstrat, y)
  mvec = as.matrix(aggregate(y~fstrat, data = dd, mean)[,-1])
  boottP = 0
  boottW = 0
  boottL= 0
  boottR = 0
  datnew = cbind(strat, y)
  if(nrep > 0)
  {
    aall = boot(data = datnew, statistic = ManovaPre, R = nrep, strata = datnew[,1],  mvec = mvec)$t
    atP= aall[,1]
    atL = aall[,2]
    atW= aall[,3]
    atR = aall[,4]
    boottP = mean(atP <= tP)
    boottL = mean(atL <= tL)
    boottW = mean(atW <= tW)
    boottR = mean(atR <= tR)
  }
  c(tP = tP, tL = tL, tW = tW, tR = tR, boottP = boottP, boottL = boottL, boottW = boottW, boottR = boottR)
}




#Median

PermMedian <- function(data, ind, mvec=rep(0, length(unique(data[,1]))), nrep = 1000)
{
  require(mvtnorm)
  require(boot)
  require(emulator)
  require(parallel)
  require(CompQuadForm)
  require(multcomp)
  MedPre <- function(data, ind, mvec = rep(0, length(unique(data[,1]))), perm = F)
  {
    require(mvtnorm)
    require(boot)
    require(emulator)
    require(parallel)
    require(CompQuadForm)
    require(multcomp)
    y = data[ind,-1]
    strat = data[ind,1]
    if(perm == T)
    {
      strat = data[,1]
    }
    y = y - mvec[strat]
    meds = tapply(y, strat, median)
    dim = meds[1] - meds[2]
    
    y1s = sort(y[strat==1])
    y2s = sort(y[strat==2])
    
    n1 = length(y1s)
    n2 = length(y2s)
    l1 = 1:n1
    l2 = 1:n2
    
    p1 = pbinom((n1-1)/2, n1, (l1-1)/n1) - pbinom((n1-1)/2, n1, l1/n1)
    p2 = pbinom((n2-1)/2, n2, (l2-1)/n2) - pbinom((n2-1)/2, n2, l2/n2)
    
    v1 = n1*sum((y1s - meds[1])^2*p1)
    v2 = n2*sum((y2s-meds[2])^2*p2)
    stu = dim/sqrt(v1/n1 +v2/n2)
    c(dim = dim, stu = stu)
  }
  y = data[ind,-1]
  strat = data[,1]
  meds = tapply(y, strat, median)
  dim = meds[1] - meds[2]
  y1s = sort(y[strat==1])
  y2s = sort(y[strat==2])
  
  n1 = length(y1s)
  n2 = length(y2s)
  l1 = 1:n1
  l2 = 1:n2
  
  p1 = pbinom((n1-1)/2, n1, (l1-1)/n1) - pbinom((n1-1)/2, n1, l1/n1)
  p2 = pbinom((n2-1)/2, n2, (l2-1)/n2) - pbinom((n2-1)/2, n2, l2/n2)
  
  v1 = n1*sum((y1s - meds[1])^2*p1)
  v2 = n2*sum((y2s-meds[2])^2*p2)
  stu = dim/sqrt(v1/n1 +v2/n2)
  mvec = meds
  datnew = cbind(strat, y)
  bootdim = 0
  bootstu = 0
  if(nrep > 0)
  {
    aall = boot(data = datnew, statistic = MedPre, R = nrep, strata = datnew[,1],  mvec = mvec)$t
    adim= aall[,1]
    astu = aall[,2]
    bootdim = mean(adim<= dim)
    bootstu  = mean(astu <= stu)
  }
  c(dim = dim, stu = stu, bootdim = bootdim, bootstu = bootstu)
}





