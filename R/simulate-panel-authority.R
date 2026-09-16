#' Simulate a small multivariate panel data set
#'
#' Generates synthetic multivariate panel trajectories with cluster structure
#' and optional observation, cell, time-point, segment, multivariate-segment,
#' and global sensor-time contamination. This is the simulation generator used
#' for package examples; large manuscript simulations should not be run during
#' package checks.
#'
#' @param n Number of units.
#' @param tt Number of time points.
#' @param p Number of variables.
#' @param k Number of latent clusters.
#' @param rho AR(1) correlation parameter for noise.
#' @param noise_sd Innovation standard deviation.
#' @param observation_rate Whole-observation contamination rate.
#' @param cell_rate Cell contamination rate.
#' @param time_rate Unit-time contamination rate.
#' @param segment_rate Unit-variable segment contamination rate.
#' @param multivariate_segment_rate Unit-time segment contamination rate.
#' @param separation Cluster separation multiplier.
#' @param anomaly_scale Contamination shift multiplier.
#' @param cluster_probs Optional cluster probabilities.
#' @param global_fault Whether to include one global sensor-time fault.
#' @param controlled_overlap Whether to add structured overlapping faults.
#' @param seed Integer random seed.
#' @return A list with `x`, true cluster labels, true means, and contamination
#'   truth components.
#' @export
simulate_mtr_panel <- function(n = 90L, tt = 30L, p = 8L, k = 3L,
                               rho = 0.5, noise_sd = 0.35,
                               observation_rate = 0.06,
                               cell_rate = 0.004,
                               time_rate = 0.03,
                               segment_rate = 0.04,
                               multivariate_segment_rate = 0.04,
                               separation = 1,
                               anomaly_scale = 1,
                               cluster_probs = NULL,
                               global_fault = TRUE,
                               controlled_overlap = TRUE,
                               seed = 20260827L) {
  set.seed(seed)
  if(is.null(cluster_probs)) z <- rep(seq_len(k), length.out = n) else {
    stopifnot(length(cluster_probs)==k,all(cluster_probs>0))
    cluster_probs<-cluster_probs/sum(cluster_probs);counts<-floor(n*cluster_probs)
    if(sum(counts)<n)counts[seq_len(n-sum(counts))]<-counts[seq_len(n-sum(counts))]+1L
    z<-rep(seq_len(k),counts)
  }
  z <- sample(z)
  grid <- seq(0, 1, length.out = tt)
  means <- array(0, c(k, tt, p))
  for (cl in seq_len(k)) for (j in seq_len(p)) {
    amp <- 0.7 + 0.12 * j
    means[cl, , j] <- separation*switch(as.character(cl),
      "1" = amp * (2 * grid - 1) + 0.15 * sin(2*pi*grid),
      "2" = -amp * (2 * grid - 1) + 0.15 * cos(2*pi*grid),
      0.9 * amp * exp(-((grid - 0.55)/0.18)^2) - 0.35)
  }
  x <- array(0, c(n, tt, p))
  for (i in seq_len(n)) for (j in seq_len(p)) {
    e <- numeric(tt); e[1] <- rnorm(1, 0, noise_sd / sqrt(1-rho^2))
    for (t in 2:tt) e[t] <- rho * e[t-1] + rnorm(1, 0, noise_sd)
    x[i, , j] <- means[z[i], , j] + e
  }
  truth <- list(observation=integer(), cell=matrix(integer(), 0, 3),
                time=matrix(integer(), 0, 2), sensor_segments=data.frame(),
                multivariate_segments=data.frame(), global=matrix(integer(), 0, 2))
  nobs <- max(1L, round(observation_rate*n)); obs <- sample(seq_len(n), nobs)
  x[obs, , ] <- x[obs, , ] + 5*anomaly_scale
  truth$observation <- obs
  nc <- max(1L, round(cell_rate*n*tt*p))
  cells <- cbind(sample(n, nc, TRUE), sample(tt, nc, TRUE), sample(p, nc, TRUE))
  for (r in seq_len(nrow(cells))) x[cells[r,1], cells[r,2], cells[r,3]] <- x[cells[r,1], cells[r,2], cells[r,3]] + sample(c(-8,8),1)*anomaly_scale
  truth$cell <- cells
  ntm <- max(1L, round(time_rate*n)); tu <- sample(setdiff(seq_len(n), obs), ntm)
  times <- sample(3:(tt-2), ntm, TRUE)
  for (r in seq_along(tu)) x[tu[r], times[r], ] <- x[tu[r], times[r], ] + 6*anomaly_scale
  truth$time <- cbind(unit=tu, time=times)
  nseg <- max(1L, round(segment_rate*n)); su <- sample(setdiff(seq_len(n), c(obs,tu)), nseg)
  for (r in seq_along(su)) {
    v <- sample(seq_len(p),1); st <- sample(3:(tt-6),1); en <- st+4
    x[su[r], st:en, v] <- x[su[r], st:en, v] + 6*anomaly_scale
    truth$sensor_segments <- rbind(truth$sensor_segments, data.frame(unit=su[r],variable=v,start=st,end=en))
  }
  nms <- max(1L, round(multivariate_segment_rate*n))
  mu <- sample(setdiff(seq_len(n), c(obs,su,tu)), nms)
  for (r in seq_along(mu)) {
    if (controlled_overlap && r==1L && nrow(truth$time)) {
      mu[r] <- truth$time[1,"unit"]
      st <- max(2L, min(tt-5L, truth$time[1,"time"]-2L))
    } else st <- sample(3:(tt-6),1)
    en <- st+4L
    x[mu[r],st:en,] <- x[mu[r],st:en,] + 6*anomaly_scale
    truth$multivariate_segments <- rbind(truth$multivariate_segments,
      data.frame(unit=mu[r],start=st,end=en))
  }
  if (global_fault) {
    gt <- round(tt*0.7); gp <- min(2L,p); x[,gt,gp] <- x[,gt,gp]+5*anomaly_scale
    truth$global <- matrix(c(gt,gp),1,2,dimnames=list(NULL,c("time","variable")))
  }
  overlap_manifest <- data.frame(type=character(),detail=character())
  if (controlled_overlap) {
    # A cell inside a sensor segment: the segment must take precedence.
    q <- truth$sensor_segments[1,]
    oc <- c(q$unit, q$start+2L, q$variable)
    x[oc[1],oc[2],oc[3]] <- x[oc[1],oc[2],oc[3]] + 8*anomaly_scale
    truth$cell <- rbind(truth$cell,oc)
    overlap_manifest <- rbind(overlap_manifest,data.frame(type="cell_in_sensor_segment",detail=paste(oc,collapse=":")))
    # A sensor segment inside an observation outlier: observation takes precedence.
    ou <- obs[1]; ov <- sample(seq_len(p),1); os <- sample(3:(tt-6),1); oe <- os+4L
    x[ou,os:oe,ov] <- x[ou,os:oe,ov] + 6*anomaly_scale
    truth$sensor_segments <- rbind(truth$sensor_segments,data.frame(unit=ou,variable=ov,start=os,end=oe))
    overlap_manifest <- rbind(overlap_manifest,data.frame(type="sensor_segment_in_observation",detail=paste(ou,ov,os,oe,sep=":")))
    # A cell at the global sensor-time coordinate: global fault takes precedence.
    gu <- setdiff(seq_len(n),obs)[1]; gc <- c(gu,truth$global[1,1],truth$global[1,2])
    x[gc[1],gc[2],gc[3]] <- x[gc[1],gc[2],gc[3]] - 8*anomaly_scale
    truth$cell <- rbind(truth$cell,gc)
    overlap_manifest <- rbind(overlap_manifest,data.frame(type="cell_in_global_fault",detail=paste(gc,collapse=":")))
    overlap_manifest <- rbind(overlap_manifest,data.frame(type="temporal_point_in_multivariate_segment",
      detail=paste(truth$time[1,],collapse=":")))
  }
  truth$overlap_manifest <- overlap_manifest
  list(x=x, true_cluster=z, true_means=means, truth=truth)
}
