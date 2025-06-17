`%||%` <- function(a, b) if (!is.null(a)) a else b



fn.Stage2 <- function(local_data, local_wt.v)
{

  if (is.null(dim(local_data$Y))) {
    local_data$Y <- matrix(local_data$Y, ncol = 1)
  }

  featureNames = c("mean","sd", "median")
  moments.ar = array(,c(length(featureNames),local_data$K, ncol(local_data$Y)))
  dimnames(moments.ar) = list(
    featureNames,
    paste("group", 1:local_data$K),
    colnames(local_data$Y) %||% paste("Y", seq_len(ncol(local_data$Y)))
  )


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
    if (sum(flag.z) == 0)
      stop(paste("No samples found for group", z, "in Stage 2."))

    moments.ar[3,z,] = apply(local_data$Y[flag.z, , drop = FALSE], 2, fn.PMF_quantile, q = 0.5, mass = local_wt.v[flag.z])

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


bootStrap <- function(data, moments.ar, moreFeatures.v, B, method, method_wt.v, naturalGroupProp, num.random.b, gammaMin, gammaMax, seed, verbose) {

  local_data <- data
  featureNames <- c("mean", "sd", "median")

  collatedMoments.ar = array(NA, c(dim(moments.ar), B))
  dimnames(collatedMoments.ar) = list(
    featureNames,
    paste("group", 1:local_data$K),
    colnames(local_data$Y) %||% paste("Y", seq_len(ncol(local_data$Y))),
    paste0("b", 1:B)
  )

  collatedMoreFeatures.mt = array(NA, c(length(moreFeatures.v), B))
  dimnames(collatedMoreFeatures.mt) = list(names(moreFeatures.v), paste0("b", 1:B))

  ## ESS
  collatedESS = rep(NA, B)

  for (b in seq_len(B)) {
    if (verbose) {
      cat("*************************\n")
      cat(paste("Bootstrap::", b), "\n")
      cat("*************************\n")
    }

    data.b = checkSample(data.in = data, wt_in.v = rep(1, data$N))
    output1.b = balancing.weights(S = data.b$S, Z = data.b$Z, X = data.b$X,
                                  method, naturalGroupProp, num.random.b,
                                  gammaMin, gammaMax, seed, verbose)
    data.b = checkSample(data.in = data, wt_in.v = output1.b$wt.v)

    # Safe fn.Stage2 call
    tryCatch({
      c(moments_b.ar, moreFeatures_b.v) %<-% fn.Stage2(local_data = data.b, local_wt.v = rep(1, data.b$N))
      collatedMoments.ar[,,,b] = moments_b.ar
      collatedMoreFeatures.mt[,b] = moreFeatures_b.v
      collatedESS[b] = output1.b$percentESS
    }, error = function(e) {
      if (verbose) {
        message(sprintf("Skipping bootstrap sample %d due to error in fn.Stage2: %s", b, e$message))
      }
    })
  }

  list(collatedMoments.ar, collatedMoreFeatures.mt, collatedESS)
}


