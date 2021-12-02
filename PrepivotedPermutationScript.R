require(mvtnorm)
require(boot)
require(emulator)
require(parallel)
require(CompQuadForm)
require(multcomp)

source("PrepivotedPermutation.R")

#Note: multicore (used within the boot function to run permutations in parallel) 
#may not work on Windows machines. You can either say parallel = "no" in the 
#function call, or use the option parallel == "snow"

#if you don't want to use the bootstrap prepivoted test statistics,
#set nbootstrap = 0 and ignore the bootstrap prepivoted output.

###########################
#Examples of prepivoted permutation tests 
#(Following Sections 5.2-5.4 of the manuscript;
#Sections A1-A2 of the supplement.)
############################

#############################
#Section 5.2: Univariate two-sample problem
#############################
#generate data, define parameters for conducting permutation test
nperm = 999
nbootstrap = 500
N = 50
n1 = floor(.3*N)
n2 = N - n1
SdT = 
X1 =  -(rexp(n1, 1/8)-8)
X2 =  ((rexp(n2, 1/5)-5))
dat = cbind(c(rep(c(1,2), c(n1, n2))), c(X1, X2))

#conduct prepivoted permutation test
prepivoted_permutation = boot(data = dat, statistic = PermAllOS, R = nperm, sim = "permutation",  parallel = "multicore", ncpus =  detectCores(), mvec = 0, nrep = nbootstrap)

#extract observed statistics (pstat) and permutation distributions, then calculate
#permutation p-values
#dim = difference in means
#stu = studentized difference in means
#edge = Gaussian prepivoted version of stu with an Edgeworth correction
#boot = bootstrap prepivoted transforms

#pvaldim will NOT provide an asymptotically valid 
#p-value for testing equality of means. The rest will.

pstat = prepivoted_permutation$t0
permdistdim = prepivoted_permutation$t[,1]
permdiststu = prepivoted_permutation$t[,2]
permdistedge = prepivoted_permutation$t[,3]
permdistbootdim = prepivoted_permutation$t[,4]
permdistbootstu = prepivoted_permutation$t[,5]
permdistbootedge = prepivoted_permutation$t[,6]

statdim = pstat[1]
statstu = pstat[2]
statedge = pstat[3]
statbootdim = pstat[4]
statbootstu= pstat[5]
statbootedge= pstat[6]

pvaldim = (1+sum(permdistdim >= statdim))/(1+nperm)
pvalstu = (1+sum(permdiststu >= statstu))/(1+nperm)
pvaledge = (1+sum(permdistedge >= statedge))/(1+nperm)
pvalbootdim = (1+sum(permdistbootdim >= statbootdim))/(1+nperm)
pvalbootstu = (1+sum(permdistbootstu >= statbootstu))/(1+nperm)
pvalbootedge = (1+sum(permdistbootedge >= statbootedge))/(1+nperm)


############
#Multivariate
############

#generate data, define parameters for conducting permutation test
nperm = 999
nbootstrap = 200
N = 150
n1 = floor(.3*N)
n2 = N - n1
d=15
Sigma1 = diag(1)
Sigma2 = matrix(.8, d, d)
diag(Sigma2) = 1
X1 = exp(rmvnorm(n1, mean= rep(0, d)))
X2 = rmvnorm(n2, mean = rep(exp(.5), d), sigma = Sigma2)

dat = cbind(c(rep(c(1,2), c(n1, n2))), rbind(X1, X2))

#conduct prepivoted permutation test
prepivoted_permutation = boot(data = dat, statistic = PermAllMult, R = nperm, sim = "permutation",  parallel = "multicore", ncpus =  detectCores(), mvec = 0, nrep = nbootstrap)

#extract observed statistics (pstat) and permutation distributions, then calculate
#permutation p-values

#unpool = unpooled Hotelling
#pool = pooled Hotelling
#max = max absolute t statistic
#pre = Gaussian prepivoted
#boot = bootstrap prepivoted
#(for instance, "bootprepool" takes the bootstrap prepivoted transform of the 
#Gaussian prepivoted transform of the pooled Hotelling Test)

#pvalpool and pvalmax will NOT provide an asymptotically valid 
#p-value for testing equality of means. The rest will.

pstat = prepivoted_permutation$t0
permdistunpool = prepivoted_permutation$t[,1]
permdistpool = prepivoted_permutation$t[,2]
permdistmax = prepivoted_permutation$t[,3]
permdistprepool = prepivoted_permutation$t[,4]
permdistpremax = prepivoted_permutation$t[,5]
permdistbootunpool = prepivoted_permutation$t[,6]
permdistbootpool = prepivoted_permutation$t[,7]
permdistbootmax = prepivoted_permutation$t[,8]
permdistbootprepool = prepivoted_permutation$t[,9]
permdistbootpremax = prepivoted_permutation$t[,10]

statunpool = pstat[1]
statpool = pstat[2]
statmax = pstat[3]
statprepool = pstat[4]
statpremax = pstat[5]
statbootunpool = pstat[6]
statbootpool = pstat[7]
statbootmax = pstat[8]
statbootprepool = pstat[9]
statbootpremax = pstat[10]

pvalunpool = (1+sum(permdistunpool >= statunpool))/(1+nperm)
pvalpool = (1+sum(permdistpool >= statpool))/(1+nperm)
pvalmax = (1+sum(permdistmax >= statmax))/(1+nperm)
pvalprepool = (1+sum(permdistprepool >= statprepool))/(1+nperm)
pvalpremax = (1+sum(permdistpremax >= statpremax))/(1+nperm)

pvalbootunpool = (1+sum(permdistbootunpool >= statbootunpool))/(1+nperm)
pvalbootpool = (1+sum(permdistbootpool >= statbootpool))/(1+nperm)
pvalbootmax = (1+sum(permdistbootmax >= statbootmax))/(1+nperm)
pvalbootprepool = (1+sum(permdistbootprepool >= statbootprepool))/(1+nperm)
pvalbootpremax = (1+sum(permdistbootpremax >= statbootpremax))/(1+nperm)

##################
#Analysis of variance (ANOVA)
##################

#generate data, define parameters for conducting permutation test

nperm = 999
nbootstrap = 200
N = 100 
k = 4
#The matrix C will contain contrasts for pairwise differences
#(used for Tukey-Kramer statistic)
C = matrix(0, k, k*(k-1)/2)
ii = 1
for(j in 1:(k-1))
{
  for(jj in (j+1):k)
  {
    C[j,ii] = 1
    C[jj, ii] = -1
    ii = ii+1
  }
}
n = floor(c(.1*N, .2*N, .3*N, .4*N))
N = sum(n)
gg = c(rep(1:k, n))
sigma = c(.7,.55, .4, .25)
Z = exp(rnorm(N, sd = sigma[gg])) - exp(sigma[gg]^2/2)
dat = cbind(gg, Z)

#conduct prepivoted permutation test
prepivoted_permutation = boot(data = dat, statistic = PermAnova, R = nperm, sim = "permutation",  parallel = "multicore", ncpus =  detectCores(), mvec = rep(0, k), C = C, nrep = nbootstrap)

#extract observed statistics (pstat) and permutation distributions, then calculate
#permutation p-values

#fstat = F statistic
#tukey = Tukey-Kramer
#CR = Robust Statistic in equation (3.2) of Chung and Romano (2013)
#boot = bootstrap prepivoted transforms

#pvalfstat and pvaltukey will NOT provide an asymptotically valid 
#p-value for testing equality of means. The rest will.


pstat = prepivoted_permutation$t0
permdistfstat = prepivoted_permutation$t[,1]
permdisttukey = prepivoted_permutation$t[,2]
permdistCR = prepivoted_permutation$t[,3]
permdistbootfstat = prepivoted_permutation$t[,4]
permdistboottukey = prepivoted_permutation$t[,5]
permdistbootCR = prepivoted_permutation$t[,6]

statfstat = pstat[1]
stattukey = pstat[2]
statCR = pstat[3]
statbootfstat = pstat[4]
statboottukey = pstat[5]
statbootCR = pstat[6]

pvalfstat = (1+sum(permdistfstat >= statfstat))/(1+nperm)
pvaltukey = (1+sum(permdisttukey >= stattukey))/(1+nperm)
pvalCR = (1+sum(permdistCR >= statCR))/(1+nperm)
pvalbootfstat = (1+sum(permdistbootfstat >= statbootfstat))/(1+nperm)
pvalboottukey = (1+sum(permdistboottukey >= statboottukey))/(1+nperm)
pvalbootCR = (1+sum(permdistbootCR >= statbootCR))/(1+nperm)

##################
#Multivariate analysis of variance (MANOVA)
##################

#generate data, define parameters for conducting permutation test

nperm = 999
nbootstrap = 200
N = 200
k = 4
d=10
n = floor(c(.1*N, .2*N, .3*N, .4*N))
N = sum(n)
gg = c(rep(1:k, n))
sigma2 = (c(1,.8,.6,.4))
corre = c(.3,.5,.7,.9)
Y = matrix(0, N, d)
for(kk in 1:k)
{
  Sigma1 = matrix(corre[kk], d, d)
  diag(Sigma1) = 1
  Sig = Sigma1*sigma2[kk]
  Y[gg==kk,] = exp(rmvnorm(n[kk], sigma = Sig)) - matrix(exp(sigma[kk]/2), n[kk], d)
}
dat = cbind(gg, Y)

#conduct prepivoted permutation test
prepivoted_permutation = boot(data = dat, statistic = PermManova, R = nperm, sim = "permutation",  parallel = "multicore", ncpus =  detectCores(), mvec = 0, nrep = nbootstrap)

#extract observed statistics (pstat) and permutation distributions, then calculate
#permutation p-values

#tP = Pillai-Bartlett Trace
#tL = Lawley-Hotelling Trace
#tW = Wilk's Lambda (we take the negative internally to look in the right tail)
#tR = Roy's Largest Root
#boot = bootstrap prepivoted transform

#pval[tP, tL, tW, tR] will NOT provide an asymptotically valid 
#p-value for testing equality of means. The bootstrap
#prepivoted transforms will


pstat = prepivoted_permutation$t0
permdisttP = prepivoted_permutation$t[,1]
permdisttL = prepivoted_permutation$t[,2]
permdisttW = prepivoted_permutation$t[,3]
permdisttR = prepivoted_permutation$t[,4]
permdistboottP = prepivoted_permutation$t[,5]
permdistboottL = prepivoted_permutation$t[,6]
permdistboottW = prepivoted_permutation$t[,7]
permdistboottR = prepivoted_permutation$t[,8]

stattP = pstat[1]
stattL = pstat[2]
stattW = pstat[3]
stattR = pstat[4]
statboottP = pstat[5]
statboottL = pstat[6]
statboottW = pstat[7]
statboottR = pstat[8]

pvaltP = (1+sum(permdisttP >= stattP))/(1+nperm)
pvaltL = (1+sum(permdisttL >= stattL))/(1+nperm)
pvaltW = (1+sum(permdisttW >= stattW))/(1+nperm)
pvaltR = (1+sum(permdisttR >= stattR))/(1+nperm)

pvalboottP = (1+sum(permdistboottP >= statboottP))/(1+nperm)
pvalboottL = (1+sum(permdistboottL >= statboottL))/(1+nperm)
pvalboottW = (1+sum(permdistboottW >= statboottW))/(1+nperm)
pvalboottR = (1+sum(permdistboottR >= statboottR))/(1+nperm)

################
#Median
###############

#generate data, define parameters for conducting permutation test

nperm = 999
nbootstrap = 500
n1 = 13
n2 = 13
N = n1 + n2
X1 = rnorm(n1)
X2 = rnorm(n2, sd = 5)

dat = cbind(c(rep(c(1,2), c(n1, n2))), c(X1, X2))

#conduct prepivoted permutation test
prepivoted_permutation = boot(data = dat, statistic = PermMedian, R = nperm, sim = "permutation",  parallel = "multicore", ncpus =  detectCores(), mvec = rep(0, length(unique(dat[,1]))), nrep = nbootstrap)

#extract observed statistics (pstat) and permutation distributions, then calculate
#permutation p-values

#dim = difference in medians
#stu = studentized difference in medians
#boot = bootstrap prepivoted transformations

pstat = prepivoted_permutation$t0
permdistdim = prepivoted_permutation$t[,1]
permdiststu = prepivoted_permutation$t[,2]
permdistbootdim = prepivoted_permutation$t[,3]
permdistbootstu = prepivoted_permutation$t[,4]

statdim = pstat[1]
statstu = pstat[2]
statbootdim = pstat[3]
statbootstu = pstat[4]

pvaldim = (1+sum(permdistdim >= statdim))/(1+nperm)
pvalstu = (1+sum(permdiststu >= statstu))/(1+nperm)
pvalbootdim = (1+sum(permdistbootdim >= statbootdim))/(1+nperm)
pvalbootstu = (1+sum(permdistbootstu >= statbootstu))/(1+nperm)




