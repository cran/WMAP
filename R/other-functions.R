#' @importFrom caret confusionMatrix
#' @importFrom randomForest randomForest
#' @importFrom stats cor optim quantile rgamma
#' @importFrom utils head



is.natural <-
  function(x, tol = .Machine$double.eps^0.5)  return((abs(x - round(x)) < tol) | (x>0))


is.probvector <- function(x, tol = .Machine$double.eps^0.5)  (abs(sum(x) - 1) < tol) & (sum(x<0)==0)


fn.rdirichlet <- function(alpha)
{

  vec = rgamma(n=length(alpha), shape=alpha, rate=1)
  vec/sum(vec)

}



fn.PMF_quantile <- function(val, q, mass)
{
  pmf = mass/sum(mass)

  indx = order(val)

  val2 = val[indx]
  pmf2 = pmf[indx]

  cdf = cumsum(pmf2)

  indx.q = min(which(cdf > q))

  if (indx.q==1)
  {ret = min(val2)
  }
  if (indx.q>1)
  {ret = mean(c(val2[indx.q-1],val2[indx.q]))
  }

  ret
}


##########################


preStage1 <- function(data.local)
{
  data.local$N.sz = array(,c(data.local$J, data.local$K))
  data.local$indx.sz.lst = array(list(NA),c(data.local$J,data.local$K))

  for (s in 1:data.local$J)
    for (z in 1:data.local$K)
    {data.local$indx.sz.lst[[s,z]] = which((data.local$S==s)&(data.local$Z==z))
    data.local$N.sz[s,z] = length(data.local$indx.sz.lst[[s,z]])
    }

  data.local$N.z = colSums(data.local$N.sz)
  data.local$N.s = rowSums(data.local$N.sz)

  data.local$which.binary = as.numeric(which(unlist(lapply(apply(data.local$X,2,unique),length))==2))

  data.local

}


check.factor <- function(factor.vec, num)
{
  err = FALSE

  numeric.vec = as.numeric(levels(factor.vec))

  if ((min(numeric.vec) != 1) | (max(numeric.vec) != num) | (length(numeric.vec) != num))
  {
    err = TRUE
  }

  err
}


fn.userInputChecks_1 <- function(method, flexorFlag, gammaMin, gammaMax, num.random, naturalGroupProp)
{
  if (!(method %in% c("IC", "IGO", "FLEXOR", "ic", "igo", "flexor")))
    stop("Not an implemented method.")

  if (flexorFlag){

    if ((gammaMin<=0) | (gammaMax>=1) | (gammaMin>=gammaMax))
      stop("gammaMin and gammaMax must be proportions with gammaMin<gammaMax.")

    if ((!is.null(num.random)) & (!is.natural(num.random)))
      stop("num.random must be a natural number.")

    if (!is.probvector(naturalGroupProp))
      stop("naturalGroupProp must be a probability vector with no zeros.")

  } # if (flexorFlag)



}


fn.userInputChecks_2 <- function(data, B=NULL)
{
  text = NULL

  if (!is.null(B))
  {
    if (!is.natural(B))
    {addText = "B cannot be negative or a non-integer. "
    text = paste(text, addText, sep="")
    }
  }


  if ((length(data$S) != length(data$Z)) | (length(data$S) != nrow(data$X)) | (!is.matrix(data$X))){
    addText = "(S,Z,X) have incompatible dimensions and/or X is not a data matrix. "
    text = paste(text, addText, sep="")
  }


  if ((sum(is.na(data$S)) + sum(is.na(data$Z)) + sum(is.na(data$X)) + sum(is.na(data$Y)))>0){
    addText = "(S,Z,X,Y) has missing values. "
    text = paste(text, addText, sep="")
  }

  if (check.factor(factor.vec=data$S, num=data$J)){
    addText = paste("S should have factor levels 1,...", data$J, ". ", sep="")
    text = paste(text, addText, sep="")
  }

  if (check.factor(factor.vec=data$Z, num=data$K))
  {addText = paste("Z should have factor levels 1,...", data$K, ". ", sep="")
  text = paste(text, addText, sep="")
  }

  if (!is.null(text))
    stop(text)
}



fn.userInputChecks_3 <- function(data)
{

  if (length(data$S) != nrow(data$Y))
    stop("Y has incompatible dimensions")

}



fn.clip <- function(vec, gammaMin, gammaMax)
{
  vec[vec < gammaMin] = gammaMin
  vec[vec > gammaMax] = gammaMax

  vec/sum(vec)
}



fn.createListData <- function(S, Z, X, Y=NULL)
{

  data = NULL

  S ->  data$S

  Z ->  data$Z

  X ->  data$X

  if (!is.null(Y))
  {Y ->  data$Y}

  data$N = length(data$S)

  data$J = max(as.numeric(data$S))

  data$K = max(as.numeric(data$Z))

  data$p = ncol(data$X)

  data = preStage1(data)

  data

}



fn.SVD <- function(data, cutoff)
  {

    matX = data.matrix(data$X)
    SVD = svd(matX); eigen.v = SVD$d
    maxx= min(which(cumsum(eigen.v)/sum(eigen.v)>cutoff))

    matX2 = array(0, dim(matX))
    for (q in 1:maxx)
      matX2 = matX2 + SVD$d[q] * (SVD$u[,q] %*% t(SVD$v[,q]))

    cor(as.vector(matX), as.vector(matX2))

    matX2
  }


fn.RF_ZFirst <- function(data, svdCutoff)
{
  parm = NULL

  parm$e$z = array(,c(data$N, data$K))

  df = data.frame(cbind(data$Z, data$X))
  names(df) = c("Z", paste("X", 1:data$p, sep=""))

  matX2 = fn.SVD(data, cutoff=svdCutoff)

  df = data.frame(cbind(data$Z, matX2))
  names(df) = c("Z", paste("X", 1:data$p, sep=""))

  tmp1 = randomForest(as.factor(Z)~., data=df)

  parm$rf_z = tmp1

  parm$e$z = tmp1$votes

  confusionMatrix(data=as.factor(apply(parm$e$z,1,which.max)), reference=data$Z)
  apply(parm$e$z,2,summary)

  parm$rf_sGiven_z = rep(list(NA), data$K)

  parm$e$sGiven_z = array(0,c(data$N, data$J, data$K))

  for (z in 1:data$K)
    {indx.z = which(data$Z == z)

    df = data.frame(cbind(data$S[indx.z], matX2[indx.z,]))
    names(df) = c("S", paste("X", 1:data$p, sep=""))

    tmp1 = randomForest(as.factor(S)~., data=df, xtest=data.frame(data$X[-indx.z,]))

    parm$rf_sGiven_z[[z]] = tmp1

    cols = as.numeric(dimnames(tmp1$test$votes)[[2]])

    tmp1$votes[which(!is.finite(tmp1$votes))] = 1/length(cols)
    tmp1$test$votes[which(!is.finite(tmp1$test$votes))] = 1/length(cols)

    parm$e$sGiven_z[-indx.z,cols,z] = tmp1$test$votes
    parm$e$sGiven_z[indx.z,cols,z] = tmp1$votes

    parm$e$sGiven_z[,,z] = parm$e$sGiven_z[,,z]/rowSums(parm$e$sGiven_z[,,z])

    unique(apply(parm$e$sGiven_z[,,z],1,sum))

    confusionMatrix(data=as.factor(apply(parm$e$sGiven_z[indx.z,,z],1,which.max)), reference=data$S[indx.z])
    apply(parm$e$sGiven_z[,,z],2,summary)
    }

  parm$e$sz = array(,c(data$N, data$K, data$J))

  for (z in 1:data$K)
  {
    parm$e$sz[,z,] = parm$e$sGiven_z[,,z] * parm$e$z[,z]
  }


  parm

}



fn.Estimate_mgPS <- function(data)
{

  parm = fn.RF_ZFirst(data, svdCutoff=.6)

  eps = 1e-2

  tmp.ar = parm$e$sz
  tmp.ar[which(tmp.ar>(1-eps),arr.ind = TRUE)] = (1-eps)
  tmp.ar[which(tmp.ar<eps,arr.ind = TRUE)] = eps

  for (ii in 1:data$N)
    {tmp.ar[ii,,] = tmp.ar[ii,,]/sum(tmp.ar[ii,,])
    }

  parm$tmp.ar = tmp.ar

  parm$tmp1.ar = -log(tmp.ar)

  parm

}




fn.IC <- function(theta_z, gamma_s, data, parm)
{

  tmp1.ar = parm$tmp1.ar

  constant.mt = array(0, c(data$K,data$J))
  constant.mt = constant.mt + log(theta_z)
  constant.mt = t(t(constant.mt) + log(gamma_s))

  wt.v = array(,data$N)

  for (s in 1:data$J)
    for (z in 1:data$K)
    {indx.sz = data$indx.sz.lst[[s,z]]
    wt.v[indx.sz] = tmp1.ar[indx.sz,z,s] + constant.mt[z,s]
    }

  maxx = max(wt.v)

  wt.v = wt.v - maxx

  wt.v = exp(wt.v)

  wt.v = wt.v/mean(wt.v)

  parm$wt.v = wt.v

  parm$ESS = data$N/mean(parm$wt.v^2)

  parm$gamma_s = gamma_s
  parm$theta_z = theta_z

  parm
}



fn.ESS_FLEXOR <- function(sub_gamma_s, theta_z, parm, data)
{
  gamma_s = c(sub_gamma_s, 1-sum(sub_gamma_s))

  tmp1.ar = parm$tmp1.ar

  constant.mt = array(0, c(data$K,data$J))
  constant.mt = constant.mt + log(theta_z)
  constant.mt = t(t(constant.mt) + log(gamma_s))
  constant.v = as.vector(constant.mt)

  tmp2.mt = t(apply(tmp1.ar, 1, "+", 2*constant.mt) )

  max.v = apply(tmp2.mt, 1, max)

  tmp4.mt = tmp2.mt - max.v

  tmp5.mt = exp(tmp4.mt)

  tmp6.mt = tmp5.mt/ rowSums(tmp5.mt)

  wt.mt = t(t(tmp6.mt) / exp(constant.v))

  wt.v = array(,data$N)

  for (s in 1:data$J)
    for (z in 1:data$K)
    {indx.sz = data$indx.sz.lst[[s,z]]
    indx2.sz = (s-1)*data$K +  z
    wt.v[indx.sz] = wt.mt[indx.sz, indx2.sz]
    }

  wt.v = wt.v/mean(wt.v)

  parm$wt.v = wt.v

  parm$ESS = data$N/mean(parm$wt.v^2)

  -parm$ESS/data$N
}






fn.Omnibus <- function(theta_z, gamma_s, data, parm)
{

  tmp1.ar = parm$tmp1.ar

  constant.mt = array(0, c(data$K,data$J))
  constant.mt = constant.mt + log(theta_z)
  constant.mt = t(t(constant.mt) + log(as.vector(gamma_s)))
  constant.v = as.vector(constant.mt)

  tmp2.mt = t(apply(tmp1.ar, 1, "+", 2*constant.mt) )

  max.v = apply(tmp2.mt, 1, max)

  tmp4.mt = tmp2.mt - max.v

  tmp5.mt = exp(tmp4.mt)

  tmp6.mt = tmp5.mt/ rowSums(tmp5.mt)

  wt.mt = t(t(tmp6.mt) / exp(constant.v))

  wt.v = array(,data$N)

  for (s in 1:data$J)
    for (z in 1:data$K)
    {indx.sz = data$indx.sz.lst[[s,z]]
    indx2.sz = (s-1)*data$K +  z
    wt.v[indx.sz] = wt.mt[indx.sz, indx2.sz]
    }

  wt.v = wt.v/mean(wt.v)

  parm$wt.v = wt.v

  parm$ESS = data$N/mean(parm$wt.v^2)

  parm$gamma_s = gamma_s
  parm$theta_z = theta_z

  parm
}



fn.loop <- function(start_theta_z, start_gamma_s, data, parm, FUN, FUN2, tol=1e-3, optim.maxit=2, loop.maxit=50)
{

  change = Inf
  count = 0

  parm$start_gamma_s = start_gamma_s

  parm$ESS = 0

  while ((change > tol) & (count < loop.maxit))
  {
    old.parm = parm
    sub_gamma_s = head(old.parm$start_gamma_s, -1)
    tmp = optim(sub_gamma_s, FUN, gr=NULL, start_theta_z, old.parm, data, control= list(maxit=optim.maxit), method = "CG")

    sub_gamma_s = as.numeric(tmp$par)
    parm$gamma_s = c(sub_gamma_s, 1-sum(sub_gamma_s))
    parm$theta_z = start_theta_z

    parm$ESS = -tmp$value*data$N

    change = abs(parm$ESS-old.parm$ESS)/old.parm$ESS
    count = count + 1
  }

  parm$count = count
  parm$change = change

  parm = FUN2(parm$theta_z, parm$gamma_s, data, parm)

  parm
}



fn.FLEXOR <- function(start_theta_z, start_gamma_s, data, parm)
{

  parm = fn.loop(start_theta_z, start_gamma_s, data, parm, FUN=fn.ESS_FLEXOR,
                 FUN2=fn.Omnibus)


  parm
}

write_res = function(estimates, CI){
  lower_bound <- CI[, 1]  # 2.5% confidence bound
  upper_bound <- CI[, 2]  # 97.5% confidence bound

  sapply(1:length(estimates), function(i) {
    paste0(round(estimates[i],2), " (", round(lower_bound[i], 2), ",", round(upper_bound[i], 2), ")")
  })
}

write_sigma_ratio = function(output){
  B = length(output$collatedESS)
  num_outcomes = length(output$otherFeatures.v)
  sigma_ratio_est = rep(NA, num_outcomes)
  for (i in 1:num_outcomes) {
    sigma_ratio_est[i] = output$moments.ar[,,i][2,1]/output$moments.ar[,,i][2,2]
  }

  CI_sigma_ratio = matrix(0,nrow = num_outcomes, ncol = 2)
  for (i in 1:num_outcomes) {
    sigma_ratio = rep(NA,B)
    for(b in 1:B){
      sigma_ratio[b] = (output$collatedMoments.ar[,,i,b][2,1]/output$collatedMoments.ar[,,i,b][2,2])
    }
    CI_sigma_ratio[i,] = quantile(sigma_ratio, probs = c(0.025,0.975))
  }
  write_res(sigma_ratio_est,CI_sigma_ratio)
}
