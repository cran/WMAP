

fn.Stage2 <- function(local_data, local_wt.v)
{

  featureNames = c("mean","sd", "median")
  moments.ar = array(,c(length(featureNames),local_data$K, ncol(local_data$Y)))
  dimnames(moments.ar) = list(featureNames, paste("group",1:local_data$K), NULL)


  moreFeatureNames = "m1-m2"
  moreFeatures.v = array(,ncol(local_data$Y))
  names(moreFeatures.v) = moreFeatureNames

  for (z in 1:local_data$K)
  {
    flag.z = local_data$Z==z

    #\mu^{z}
    moments.ar[1,z,] = colMeans(local_data$Y *flag.z *as.vector(local_wt.v))/mean(flag.z *local_wt.v)

    moments.ar[2,z,] = colMeans(local_data$Y^2 *flag.z *as.vector(local_wt.v))/mean(flag.z *local_wt.v)

    #\sigma^{z}
    moments.ar[2,z,] = sqrt(abs(moments.ar[2,z,] - moments.ar[1,z,]^2))

    #median^{z}
    moments.ar[3,z,] = apply(local_data$Y[flag.z,], 2, fn.PMF_quantile, q=.5, mass=local_wt.v[flag.z])

  }

  #\mu^{1} - \mu^{2}
  moreFeatures.v = moments.ar[1,1,] - moments.ar[1,2,]

  list(moments.ar, moreFeatures.v)
}



checkSample <- function(data.in, wt_in.v)
{

  data.local = data.in

  okayFlag = FALSE

  while(!okayFlag)
  {
    indx.b = sample(1:data.local$N, size=data.local$N, replace=TRUE, prob=wt_in.v)

    data.local$X = data.in$X[indx.b,]
    data.local$S = data.in$S[indx.b]
    data.local$Z = data.in$Z[indx.b]
    data.local$Y = data.in$Y[indx.b,]

    data.local = preStage1(data.local)

    okayFlag =sum(data.local$N.zs==0) == 0

  }

  data.local
}



bootStrap <- function(data, moments.ar, moreFeatures.v, B, method, method_wt.v, naturalGroupProp, num.random.b, gammaMin, gammaMax, seed, verbose)
{

  collatedMoments.ar = array(,c(dim(moments.ar), B))
  dimnames(collatedMoments.ar) = c(dimnames(moments.ar), list(NULL))
  #
  collatedMoreFeatures.mt = array(,c(length(moreFeatures.v), B))
  dimnames(collatedMoreFeatures.mt) = c(names(moreFeatures.v), list(NULL))

  ## ESS
  collatedESS = rep(NA,B)


  for (b in 1:B)
  {
    if(verbose){
      print("*************************")
      print(paste("Bootstrap::",b))
      print("*************************")
    }

    data.b = checkSample(data.in=data, wt_in.v=rep(1,data$N))

    output1.b = balancing.weights(S=data.b$S, Z=data.b$Z, X=data.b$X, method, naturalGroupProp, num.random.b, gammaMin, gammaMax, seed, verbose)

    data.b = checkSample(data.in=data, wt_in.v=output1.b$wt.v)

    c(moments_b.ar, moreFeatures_b.v) %<-% fn.Stage2(local_data=data.b, local_wt.v=rep(1,data.b$N))

    collatedMoments.ar[,,,b] = moments_b.ar
    collatedMoreFeatures.mt[,b] = moreFeatures_b.v

    ## botstrap ESS
    collatedESS[b] = output1.b$percentESS

  }

  list(collatedMoments.ar, collatedMoreFeatures.mt,collatedESS)
}


