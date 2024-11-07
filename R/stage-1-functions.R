

Stage1 <- function(data, parm, method, flexorFlag, num.random, naturalGroupProp, gammaMin, gammaMax, verbose)
{

  if (!flexorFlag){
    startLabels = method

    if (method %in% c("IC","ic")){
      return.parm = fn.IC(theta_z=rep(1,data$K)/data$K, gamma_s=rep(1,data$J)/data$J, data, parm)
    }

    if (method %in% c("IGO","igo")){
      return.parm = fn.Omnibus(theta_z=rep(1,data$K)/data$K, gamma_s=rep(1,data$J)/data$J, data, parm)
    }

  } # end   if (!flexorFlag)



  if (flexorFlag){

    startLabels = c("IGO_FLEXOR","Data-driven_FLEXOR")

    num.fixed = length(startLabels)

    if (num.random > 0){
      startLabels = c(startLabels, paste("randomGamma_FLEXOR_",1:num.random,sep=""))
    }

    num.starts = num.fixed + num.random

    percentESS.v = array(,num.starts)
    names(percentESS.v) = startLabels

    gammaStart.mt = gammaEnd.mt = array(,c(num.starts, data$J))
    thetaStart.mt = thetaEnd.mt =  array(,c(num.starts, data$K))
    dimnames(gammaStart.mt) = dimnames(gammaEnd.mt) = list(startLabels, NULL)
    dimnames(thetaStart.mt) = dimnames(thetaEnd.mt) = list(startLabels, NULL)

    # FLEXOR w/ IGO initialization
    parm = fn.FLEXOR(start_theta_z=naturalGroupProp, start_gamma_s=rep(1,data$J)/data$J, data, parm)
    gammaStart.mt[1,] = rep(1,data$J)/data$J
    thetaStart.mt[1,] = naturalGroupProp
    gammaEnd.mt[1,] = parm$gamma_s
    thetaEnd.mt[1,] = parm$theta_z
    percentESS.v[1] = parm$ESS/data$N*100

    # FLEXOR w/ data-driven initialization
    parm = fn.FLEXOR(start_theta_z=naturalGroupProp, start_gamma_s=data$N.s/data$N, data, parm)
    gammaStart.mt[2,] = data$N.s/data$N
    thetaStart.mt[2,] = naturalGroupProp
    gammaEnd.mt[2,] = parm$gamma_s
    thetaEnd.mt[2,] = parm$theta_z
    percentESS.v[2] = parm$ESS/data$N*100

    if (num.random>0){

      for (cc in (num.fixed+1):num.starts){
        mass = 50
        gammaStart.mt[cc,] = fn.clip(fn.rdirichlet(alpha=mass*data$N.s/data$N), gammaMin, gammaMax)

        thetaStart.mt[cc,] = naturalGroupProp
        tmp = try(fn.FLEXOR(start_theta_z=naturalGroupProp, start_gamma_s=gammaStart.mt[cc,], data, parm), silent = TRUE)

        if (length(tmp)>1){
          parm = tmp
          gammaEnd.mt[cc,] = parm$gamma_s
          thetaEnd.mt[cc,] = parm$theta_z
          percentESS.v[cc] = parm$ESS/data$N*100
        }

        if (verbose & cc %% 10 == 0){
          print(paste("FLEXOR... estimate", cc))
        }
      } # for (cc in (num.fixed+1):num.starts)
    } # if (num.random>0)

    admissibleGamma = rowSums((gammaEnd.mt > gammaMin) & (gammaEnd.mt < gammaMax)) == data$J
    admissible = which(admissibleGamma)

    cc.max = admissible[which.max(percentESS.v[admissible])]

    return.parm = fn.FLEXOR(start_theta_z=naturalGroupProp, start_gamma_s=gammaStart.mt[cc.max,], data, parm)
  } # if (flexorFlag)

  percentESS = return.parm$ESS/data$N*100

  list(return.parm, percentESS)
}




fn.Weights <- function(data, method, flexorFlag, num.random, naturalGroupProp, gammaMin, gammaMax, verbose)
    {

    parm =  fn.Estimate_mgPS(data)

    c(return.parm, percentESS) %<-% Stage1(data, parm, method, flexorFlag, num.random, naturalGroupProp, gammaMin, gammaMax, verbose)

    list(return.parm, percentESS)

    }


