Stage1 <- function(data, parm, method, flexorFlag, num.random, naturalGroupProp, gammaMin, gammaMax, verbose) {

  if (!flexorFlag) {
    startLabels = method
    if (method %in% c("IC","ic")) {
      return.parm = fn.IC(theta_z=rep(1, data$K)/data$K, gamma_s=rep(1, data$J)/data$J, data, parm)
    }
    if (method %in% c("IGO","igo")) {
      return.parm = fn.Omnibus(theta_z=rep(1, data$K)/data$K, gamma_s=rep(1, data$J)/data$J, data, parm)
    }

  } else { # flexorFlag is TRUE

    startLabels = c("IGO_FLEXOR", "Data-driven_FLEXOR")
    num.fixed = length(startLabels)
    if (num.random > 0) {
      startLabels = c(startLabels, paste("randomGamma_FLEXOR_", seq_len(num.random), sep = ""))
    }

    num.starts = length(startLabels)

    percentESS.v = array(, num.starts)
    names(percentESS.v) = startLabels

    gammaStart.mt = gammaEnd.mt = array(, c(num.starts, data$J))
    thetaStart.mt = thetaEnd.mt = array(, c(num.starts, data$K))
    dimnames(gammaStart.mt) = dimnames(gammaEnd.mt) = list(startLabels, NULL)
    dimnames(thetaStart.mt) = dimnames(thetaEnd.mt) = list(startLabels, NULL)

    # FLEXOR with IGO initialization
    parm = fn.FLEXOR(start_theta_z = naturalGroupProp, start_gamma_s = rep(1, data$J)/data$J, data, parm)
    if (!is.null(parm)) {
      gammaStart.mt[1,] = rep(1, data$J)/data$J
      thetaStart.mt[1,] = naturalGroupProp
      gammaEnd.mt[1,] = parm$gamma_s
      thetaEnd.mt[1,] = parm$theta_z
      percentESS.v[1] = parm$ESS / data$N * 100
    }

    # FLEXOR with data-driven initialization
    parm = fn.FLEXOR(start_theta_z = naturalGroupProp, start_gamma_s = data$N.s / data$N, data, parm)
    if (!is.null(parm)) {
      gammaStart.mt[2,] = data$N.s / data$N
      thetaStart.mt[2,] = naturalGroupProp
      gammaEnd.mt[2,] = parm$gamma_s
      thetaEnd.mt[2,] = parm$theta_z
      percentESS.v[2] = parm$ESS / data$N * 100
    }

    # Random starts
    if (num.random > 0) {
      for (cc in (num.fixed + 1):num.starts) {
        mass = 50
        gammaStart.mt[cc,] = fn.clip(fn.rdirichlet(alpha = mass * data$N.s / data$N), gammaMin, gammaMax)
        thetaStart.mt[cc,] = naturalGroupProp

        tmp = try(fn.FLEXOR(start_theta_z = naturalGroupProp, start_gamma_s = gammaStart.mt[cc,], data, parm), silent = TRUE)

        valid_result = (
          is.list(tmp) &&
            !is.null(tmp$gamma_s) && length(tmp$gamma_s) == data$J &&
            !any(is.na(tmp$gamma_s)) && all(is.finite(tmp$gamma_s)) &&
            !is.null(tmp$theta_z) && length(tmp$theta_z) == data$K &&
            !is.null(tmp$ESS) && is.finite(tmp$ESS)
        )

        if (valid_result) {
          parm = tmp
          gammaEnd.mt[cc,] = parm$gamma_s
          thetaEnd.mt[cc,] = parm$theta_z
          percentESS.v[cc] = parm$ESS / data$N * 100
        } else if (verbose) {
          message(paste("FLEXOR estimate", cc, "failed or returned invalid output."))
        }

        if (verbose && cc %% 10 == 0) {
          print(paste("FLEXOR... estimate", cc))
        }
      }
    }

    admissibleGamma <- apply(gammaEnd.mt, 1, function(row) {
      if (any(is.na(row)) || any(!is.finite(row))) return(FALSE)
      all(row > gammaMin & row < gammaMax)
    })

    admissible <- which(admissibleGamma)
    if (length(admissible) == 0) {
      stop("No admissible gamma vectors found. All FLEXOR estimates failed.")
    }

    percentESS.v[is.na(percentESS.v)] <- -Inf
    cc.max <- admissible[which.max(percentESS.v[admissible])]

    return.parm <- fn.FLEXOR(start_theta_z = naturalGroupProp, start_gamma_s = gammaStart.mt[cc.max,], data, parm)
    if (is.null(return.parm) || is.null(return.parm$ESS) || !is.finite(return.parm$ESS)) {
      stop("Final FLEXOR call failed or returned invalid output.")
    }
  }

  percentESS = return.parm$ESS / data$N * 100
  list(return.parm, percentESS)
}



fn.Weights <- function(data, method, flexorFlag, num.random, naturalGroupProp, gammaMin, gammaMax, verbose) {

  parm <- fn.Estimate_mgPS(data)

  success <- FALSE
  attempts <- 3
  i <- 1

  while (!success && i <= attempts) {
    tryCatch({
      c(return.parm, percentESS) %<-% Stage1(data, parm, method, flexorFlag,
                                             num.random, naturalGroupProp,
                                             gammaMin, gammaMax, verbose)
      success <- TRUE
    }, error = function(e) {
      message(sprintf("Attempt %d failed: %s", i, e$message))
    })
    i <- i + 1
  }

  if (!success) {
    warning("All Stage1 attempts failed after 3 tries. Returning NULL.")
    return(NULL)
  }

  return(list(return.parm, percentESS))
}
