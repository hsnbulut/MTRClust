mtr_scale <- function(x, eps=1e-8) {
  stopifnot(length(dim(x))==3L); p <- dim(x)[3]; center <- scale <- numeric(p); y <- x
  for (j in seq_len(p)) {
    v <- as.vector(x[,,j]); center[j] <- median(v,na.rm=TRUE)
    scale[j] <- 1.4826*median(abs(v-center[j]),na.rm=TRUE)
    if (!is.finite(scale[j]) || scale[j]<eps) scale[j] <- IQR(v,na.rm=TRUE)/1.349
    if (!is.finite(scale[j]) || scale[j]<eps) scale[j] <- 1
    y[,,j] <- (x[,,j]-center[j])/scale[j]
  }
  list(y=y,center=center,scale=scale)
}

tukey_rho <- function(r,c=4.685) {
  u <- r/c; out <- rep(c^2/6,length(r)); inside <- is.finite(u)&abs(u)<=1
  out[inside] <- (c^2/6)*(1-(1-u[inside]^2)^3); out[!is.finite(r)] <- c^2/6
  array(out,dim=dim(r))
}

tukey_weight <- function(r,c=4.685) {
  u <- r/c; out <- numeric(length(r)); inside <- is.finite(u)&abs(u)<=1
  out[inside] <- (1-u[inside]^2)^2; array(out,dim=dim(r))
}

ari_index <- function(a,b) {
  ok <- is.finite(a)&is.finite(b); tab <- table(factor(a[ok]),factor(b[ok])); c2 <- function(z) z*(z-1)/2
  nij <- sum(c2(tab)); ai <- sum(c2(rowSums(tab))); bj <- sum(c2(colSums(tab))); n2 <- c2(sum(tab))
  if (n2==0) return(NA_real_); ex <- ai*bj/n2; den <- .5*(ai+bj)-ex
  if (abs(den)<.Machine$double.eps) return(1); (nij-ex)/den
}

robust_profile_outlyingness <- function(y,k_nn) {
  n <- dim(y)[1]; d <- matrix(0,n,n)
  for (i in seq_len(n-1L)) for (j in (i+1L):n) d[i,j] <- d[j,i] <- median(abs(y[i,,]-y[j,,]),na.rm=TRUE)
  k_nn <- max(1L,min(as.integer(k_nn),n-1L))
  score <- vapply(seq_len(n),function(i) sort(d[i,-i],partial=k_nn)[k_nn],numeric(1))
  list(score=score,distance=d)
}

smooth_low_rank_centers <- function(y,z,keep,total_w,k,latent_rank=2L,spar=.55,old_centers=NULL) {
  dd <- dim(y); tt <- dd[2]; pp <- dd[3]; centers <- array(0,c(k,tt,pp)); grid <- seq_len(tt)
  for (cl in seq_len(k)) {
    ids <- which(z==cl & keep)
    if (!length(ids)) { if (!is.null(old_centers)) centers[cl,,] <- old_centers[cl,,]; next }
    for (t in seq_len(tt)) for (p in seq_len(pp)) {
      ww <- total_w[ids,t,p]; yy <- y[ids,t,p]; good <- is.finite(yy)&is.finite(ww)&ww>0
      centers[cl,t,p] <- if (any(good)) sum(ww[good]*yy[good])/sum(ww[good]) else median(yy,na.rm=TRUE)
    }
    for (p in seq_len(pp)) centers[cl,,p] <- predict(smooth.spline(grid,centers[cl,,p],spar=spar),grid)$y
  }
  stacked <- do.call(rbind,lapply(seq_len(k),function(cl) centers[cl,,])); rr <- min(latent_rank,nrow(stacked),pp)
  sv <- svd(stacked,nu=rr,nv=rr); rec <- sv$u[,seq_len(rr),drop=FALSE] %*% diag(sv$d[seq_len(rr)],rr) %*% t(sv$v[,seq_len(rr),drop=FALSE])
  for (cl in seq_len(k)) centers[cl,,] <- rec[((cl-1L)*tt+1L):(cl*tt),,drop=FALSE]
  list(centers=centers,loadings=sv$v[,seq_len(rr),drop=FALSE])
}

profile_trimmed_loss <- function(loss,local_trim=.10) {
  vapply(seq_len(dim(loss)[1]),function(i) { v <- sort(as.vector(loss[i,,]),na.last=NA); mean(v[seq_len(max(1L,floor((1-local_trim)*length(v))))]) },numeric(1))
}

calibrate_residuals_v7 <- function(y,centers,z,keep,k,adaptive=FALSE,
  cutoff_prob=.99,cutoff_bounds=c(2.5,3.5),scale_floor=.10) {
  pp<-dim(y)[3]; scales<-matrix(1,k,pp); cutoffs<-matrix(NA_real_,k,pp)
  for(cl in seq_len(k)) for(p in seq_len(pp)) {
    ids<-which(z==cl&keep); if(!length(ids)) next
    v<-as.vector(sweep(y[ids,,p,drop=FALSE],2,centers[cl,,p],"-")); med<-median(v,na.rm=TRUE)
    sc<-1.4826*median(abs(v-med),na.rm=TRUE)
    if(!is.finite(sc)||sc<scale_floor) sc<-IQR(v,na.rm=TRUE)/1.349
    if(!is.finite(sc)||sc<scale_floor) sc<-1
    scales[cl,p]<-sc
    if(adaptive) {
      q<-unname(quantile(abs(v/sc),cutoff_prob,na.rm=TRUE,type=8))
      cutoffs[cl,p]<-min(max(q,cutoff_bounds[1]),cutoff_bounds[2])
    }
  }
  list(scales=scales,cutoffs=cutoffs)
}

calc_components_v7 <- function(y,centers,scales=NULL,cutoffs=NULL,c_cell=4.685,cell_threshold=.20,time_fraction_cut=.50,local_trim=.10) {
  n <- dim(y)[1]; k <- dim(centers)[1]; dist <- matrix(0,n,k)
  pp<-dim(y)[3]; if(is.null(scales)) scales<-matrix(1,k,pp); if(is.null(cutoffs)) cutoffs<-matrix(NA_real_,k,pp)
  cbk <- fbk <- tbk <- tsbk <- rbk <- ubk <- vector("list",k)
  for (cl in seq_len(k)) {
    r <- sweep(y,c(2,3),centers[cl,,],"-"); u<-sweep(r,3,scales[cl,],"/"); cw <- tukey_weight(u,c_cell)
    flag<-cw<cell_threshold
    if(any(is.finite(cutoffs[cl,]))) flag<-abs(u)>array(cutoffs[cl,],dim=dim(u))
    ts <- apply(flag,c(1,2),mean)
    tw <- 1-ts/time_fraction_cut
    tw[tw<0] <- 0
    dist[,cl] <- profile_trimmed_loss(tukey_rho(u,c_cell),local_trim)
    cbk[[cl]]<-cw; fbk[[cl]]<-flag; tbk[[cl]]<-tw; tsbk[[cl]]<-ts; rbk[[cl]]<-r; ubk[[cl]]<-u
  }
  list(distance=dist,cell_by_k=cbk,flag_by_k=fbk,time_by_k=tbk,time_score_by_k=tsbk,residual_by_k=rbk,standardized_residual_by_k=ubk)
}

repair_min_clusters_v12 <- function(z,keep,distance,min_cluster,k) {
  counts<-tabulate(z[keep],nbins=k);moves<-0L;moved_ids<-integer()
  while(any(counts<min_cluster)) {
    target<-which.min(counts);donors<-which(counts>min_cluster)
    if(!length(donors))return(list(ok=FALSE,z=z,counts=counts,moves=moves,moved_ids=unique(moved_ids)))
    ids<-which(keep&z%in%donors);delta<-distance[cbind(ids,rep(target,length(ids)))]-distance[cbind(ids,z[ids])]
    good<-is.finite(delta);if(!any(good))return(list(ok=FALSE,z=z,counts=counts,moves=moves,moved_ids=unique(moved_ids)))
    move_id<-ids[which.min(replace(delta,!good,Inf))];old<-z[move_id];z[move_id]<-target
    counts[old]<-counts[old]-1L;counts[target]<-counts[target]+1L;moves<-moves+1L;moved_ids<-c(moved_ids,move_id)
  }
  list(ok=TRUE,z=z,counts=counts,moves=moves,moved_ids=unique(moved_ids))
}

mtrclust_v12 <- function(x,k,alpha=.08,latent_rank=2L,spar=.55,c_cell=4.685,
  cell_threshold=.20,time_fraction_cut=.50,global_fraction_cut=.50,local_trim=.10,
  nstart=8L,max_iter=100L,tol_center=1e-5,min_cluster=NULL,min_segment=3L,
  method=c("v6_global","cluster_scale","cluster_adaptive"),cutoff_prob=.99,
  cutoff_bounds=c(2.5,3.5),scale_floor=.10,repair=TRUE,
  rollback_nonconverged=TRUE,trace_history=FALSE,
  selection_rule="objective",cycle_detection=TRUE,
  cycle_max_period=6L,cycle_tol=1e-7,seed=20260827L,verbose=FALSE,
  return_candidates=FALSE) {
  method<-match.arg(method)
  selection_rule<-match.arg(selection_rule,"objective")
  stopifnot(length(dim(x))==3L,k>=2L,alpha>=0,alpha<.5)
  scaled <- mtr_scale(x); y <- scaled$y; n <- dim(y)[1]; tt <- dim(y)[2]; pp <- dim(y)[3]; h <- n-floor(alpha*n)
  if (is.null(min_cluster)) min_cluster <- max(5L,ceiling(.35*h/k)); min_cluster <- as.integer(min_cluster)
  out_init <- robust_profile_outlyingness(y,max(floor(alpha*n)+2L,min_cluster))
  init_keep <- rank(out_init$score,ties.method="first")<=h; flat <- matrix(y,nrow=n); candidates <- list()
  rejected_starts<-0L;valid_starts<-0L;set.seed(seed)
  for (st in seq_len(nstart)) {
    if(st==1L) {
      iz<-cutree(hclust(as.dist(out_init$distance[init_keep,init_keep]),method="average"),k=k)
      init_centers<-do.call(rbind,lapply(seq_len(k),function(cl)colMeans(flat[init_keep,,drop=FALSE][iz==cl,,drop=FALSE])))
      init_type<-"robust_hclust"
    } else {
      km<-kmeans(flat[init_keep,,drop=FALSE],centers=k,nstart=1,iter.max=100)
      iz<-km$cluster;init_centers<-km$centers;init_type<-"kmeans"
    }
    z <- integer(n); z[init_keep] <- iz
    if (any(!init_keep)) {
      ids <- which(!init_keep); dd <- matrix(0,length(ids),k)
      for (cl in seq_len(k)) dd[,cl] <- rowSums((flat[ids,,drop=FALSE]-matrix(init_centers[cl,],length(ids),ncol(flat),byrow=TRUE))^2)
      z[ids] <- max.col(-dd,ties.method="random")
    }
    initial_z<-z;initial_keep_start<-init_keep
    keep <- init_keep; total_w <- array(rep(as.numeric(keep),tt*pp),dim(y))
    upd <- smooth_low_rank_centers(y,z,keep,total_w,k,latent_rank,spar); centers <- upd$centers
    objective <- center_hist <- numeric(); member_hist <- keep_hist <- repair_hist <- integer()
    converged <- FALSE; valid <- TRUE; cycle_detected<-FALSE;cycle_period<-NA_integer_;stable_run<-0L
    repair_moves_total<-0L;repair_units_all<-integer();best_iter_state<-NULL
    for (iter in seq_len(max_iter)) {
      old_z<-z; old_keep<-keep; old_centers<-centers;old_loadings<-upd$loadings
      cal<-if(method=="v6_global") list(scales=matrix(1,k,pp),cutoffs=matrix(NA_real_,k,pp)) else
        calibrate_residuals_v7(y,centers,z,keep,k,method=="cluster_adaptive",cutoff_prob,cutoff_bounds,scale_floor)
      comp <- calc_components_v7(y,centers,cal$scales,cal$cutoffs,c_cell,cell_threshold,time_fraction_cut,local_trim)
      z <- max.col(-comp$distance,ties.method="random"); dmin <- comp$distance[cbind(seq_len(n),z)]
      keep <- rank(dmin,ties.method="first")<=h;counts<-tabulate(z[keep],nbins=k);iter_moves<-0L
      if(repair){
        repa<-repair_min_clusters_v12(z,keep,comp$distance,min_cluster,k)
        if(!repa$ok){valid<-FALSE;rejected_starts<-rejected_starts+1L;break}
        z<-repa$z;counts<-repa$counts;iter_moves<-repa$moves
        repair_moves_total<-repair_moves_total+iter_moves;repair_units_all<-union(repair_units_all,repa$moved_ids)
      } else if(any(counts<min_cluster)) {
        valid<-FALSE;rejected_starts<-rejected_starts+1L;break
      }
      dmin[keep]<-comp$distance[cbind(which(keep),z[keep])]
      cell_w <- array(0,dim(y)); cell_flag<-array(FALSE,dim(y)); time_w <- time_score <- matrix(0,n,tt); residual <- array(0,dim(y))
      for (i in seq_len(n)) { cl<-z[i]; cell_w[i,,]<-comp$cell_by_k[[cl]][i,,]; cell_flag[i,,]<-comp$flag_by_k[[cl]][i,,]; time_w[i,]<-comp$time_by_k[[cl]][i,]; time_score[i,]<-comp$time_score_by_k[[cl]][i,]; residual[i,,]<-comp$residual_by_k[[cl]][i,,] }
      raw_flag <- cell_flag; gf <- apply(raw_flag[keep,,,drop=FALSE],c(2,3),mean); gw <- ifelse(gf>=global_fraction_cut,0,1)
      total_w <- cell_w; for(i in seq_len(n)) total_w[i,,] <- cell_w[i,,]*time_w[i,]*gw*as.numeric(keep[i])
      upd <- smooth_low_rank_centers(y,z,keep,total_w,k,latent_rank,spar,centers); centers <- upd$centers
      cc <- sqrt(sum((centers-old_centers)^2))/(sqrt(sum(old_centers^2))+1e-12); mc<-sum(z!=old_z); kc<-sum(keep!=old_keep)
      objective<-c(objective,sum(dmin[keep])); center_hist<-c(center_hist,cc); member_hist<-c(member_hist,mc); keep_hist<-c(keep_hist,kc);repair_hist<-c(repair_hist,iter_moves)
      if(is.null(best_iter_state)||tail(objective,1)<best_iter_state$objective)
        best_iter_state<-list(objective=tail(objective,1),centers=old_centers,loadings=old_loadings,z=z,keep=keep,iteration=iter)
      if(verbose) message(sprintf("start=%d iter=%d objective=%.6f center=%.3e membership=%d keep=%d sizes=%s",st,iter,tail(objective,1),cc,mc,kc,paste(counts,collapse=",")))
      if(mc==0L&&kc==0L&&cc<tol_center) { converged<-TRUE; break }
      stable_run<-if(mc==0L&&kc==0L)stable_run+1L else 0L
      if(cycle_detection&&stable_run>=4L&&length(objective)>=4L) {
        max_lag<-min(as.integer(cycle_max_period),floor(length(objective)/2L),floor(stable_run/2L))
        if(max_lag>=2L)for(lag in 2L:max_lag) {
          current_block<-tail(objective,lag);previous_block<-objective[(length(objective)-2L*lag+1L):(length(objective)-lag)]
          repeated<-all(abs(current_block-previous_block)<=cycle_tol*(1+abs(current_block)))
          oscillating<-(max(current_block)-min(current_block))>cycle_tol*(1+abs(mean(current_block)))
          if(repeated&&oscillating) {cycle_detected<-TRUE;cycle_period<-lag;break}
        }
        if(cycle_detected)break
      }
    }
    if(!valid||!length(objective)) next
    rolled_back<-FALSE
    if(!converged&&rollback_nonconverged&&!is.null(best_iter_state)) {
      centers<-best_iter_state$centers;upd$loadings<-best_iter_state$loadings;z<-best_iter_state$z;keep<-best_iter_state$keep;rolled_back<-TRUE
    }
    cal<-if(method=="v6_global") list(scales=matrix(1,k,pp),cutoffs=matrix(NA_real_,k,pp)) else
      calibrate_residuals_v7(y,centers,z,keep,k,method=="cluster_adaptive",cutoff_prob,cutoff_bounds,scale_floor)
    comp <- calc_components_v7(y,centers,cal$scales,cal$cutoffs,c_cell,cell_threshold,time_fraction_cut,local_trim)
    z <- max.col(-comp$distance,ties.method="first"); dmin<-comp$distance[cbind(seq_len(n),z)]; keep<-rank(dmin,ties.method="first")<=h
    counts<-tabulate(z[keep],nbins=k);final_repair_moves<-0L;final_repair_units<-integer()
    if(repair){
      repa<-repair_min_clusters_v12(z,keep,comp$distance,min_cluster,k);if(!repa$ok){rejected_starts<-rejected_starts+1L;next}
      z<-repa$z;counts<-repa$counts;final_repair_moves<-repa$moves;final_repair_units<-repa$moved_ids
      repair_moves_total<-repair_moves_total+final_repair_moves;repair_units_all<-union(repair_units_all,final_repair_units)
    } else if(any(counts<min_cluster)){rejected_starts<-rejected_starts+1L;next}
    dmin[keep]<-comp$distance[cbind(which(keep),z[keep])]
    valid_starts<-valid_starts+1L
    cell_w<-array(0,dim(y)); cell_flag<-array(FALSE,dim(y)); time_w<-time_score<-matrix(0,n,tt); residual<-array(0,dim(y))
    for(i in seq_len(n)) { cl<-z[i]; cell_w[i,,]<-comp$cell_by_k[[cl]][i,,]; cell_flag[i,,]<-comp$flag_by_k[[cl]][i,,]; time_w[i,]<-comp$time_by_k[[cl]][i,]; time_score[i,]<-comp$time_score_by_k[[cl]][i,]; residual[i,,]<-comp$residual_by_k[[cl]][i,,] }
    raw_flag<-cell_flag; gf<-apply(raw_flag[keep,,,drop=FALSE],c(2,3),mean); gflag<-gf>=global_fraction_cut
    cand<-list(objective=sum(dmin[keep]),z=z,keep=keep,centers=centers,loadings=upd$loadings,distance=comp$distance,
      cell_weights=cell_w,cell_flags=cell_flag,time_weights=time_w,time_scores=time_score,global_weights=ifelse(gflag,0,1),global_fraction=gf,residuals=residual,
      residual_scales=cal$scales,adaptive_cutoffs=cal$cutoffs,
      objective_history=objective,center_change_history=center_hist,membership_change_history=member_hist,keep_change_history=keep_hist,
      iterations=iter,converged=converged,rolled_back=rolled_back,
      cycle_detected=cycle_detected,cycle_period=cycle_period,
      termination=if(converged)"converged"else if(cycle_detected)"cycle_detected"else"max_iter",
      best_iteration=if(rolled_back)best_iter_state$iteration else iter,
      retained_cluster_sizes=counts,initialization_used=init_type,repair_enabled=repair,
      initial_z=initial_z,initial_keep=initial_keep_start,
      repair_moves=repair_moves_total,repair_units_unique=length(repair_units_all),
      final_repair_moves=final_repair_moves,final_repair_units=length(unique(final_repair_units)),
      repair_moves_history=repair_hist)
    if(trace_history)cand$iteration_trace<-data.frame(iteration=seq_along(objective),objective=objective,
      center_change=center_hist,membership_change=member_hist,keep_change=keep_hist,repair_moves=repair_hist)
    candidates[[length(candidates)+1L]]<-cand
  }
  if(!length(candidates)) stop(sprintf("No valid solution satisfied the minimum cluster-size constraint after %d starts.",nstart))
  objectives<-vapply(candidates,`[[`,numeric(1),"objective");objective_best<-which.min(objectives)
  if(length(candidates)==1L){stability<-1;selected<-1L}else{
    agreement<-diag(1,length(candidates))
    for(a in seq_len(length(candidates)-1L))for(b in (a+1L):length(candidates))
      agreement[a,b]<-agreement[b,a]<-ari_index(candidates[[a]]$z,candidates[[b]]$z)
    stability<-rowMeans(agreement)
    stable_best<-which(stability>=max(stability,na.rm=TRUE)-1e-10)
    selected<-stable_best[which.min(objectives[stable_best])]
  }
  if(selection_rule=="objective")selected<-objective_best
  best<-candidates[[selected]]
  best$selection_rule<-selection_rule;best$selected_candidate<-selected;best$objective_best_candidate<-objective_best
  best$selection_changed<-selected!=objective_best;best$selected_stability<-stability[selected]
  best$candidate_summary<-data.frame(candidate=seq_along(candidates),objective=objectives,stability=stability,
    converged=vapply(candidates,`[[`,logical(1),"converged"),cycle_detected=vapply(candidates,`[[`,logical(1),"cycle_detected"),
    termination=vapply(candidates,`[[`,character(1),"termination"),
    rolled_back=vapply(candidates,`[[`,logical(1),"rolled_back"),
    iterations=vapply(candidates,`[[`,integer(1),"iterations"),
    best_iteration=vapply(candidates,`[[`,integer(1),"best_iteration"),
    initialization=vapply(candidates,`[[`,character(1),"initialization_used"),
    repair_moves=vapply(candidates,`[[`,integer(1),"repair_moves"),
    repair_units_unique=vapply(candidates,`[[`,integer(1),"repair_units_unique"),
    selected=seq_along(candidates)==selected)
  best$nstarts_attempted<-as.integer(nstart);best$valid_starts<-valid_starts;best$rejected_starts<-rejected_starts
  raw_cell_flag<-best$cell_flags
  raw_time_flag<-best$time_scores>=time_fraction_cut
  get_runs <- function(flag,level) {
    ans<-list(); id<-0L
    if(level=="cell") for(i in seq_len(n)) for(p in seq_len(pp)) { rr<-rle(flag[i,,p]); ends<-cumsum(rr$lengths); starts<-ends-rr$lengths+1L; hit<-which(rr$values&rr$lengths>=min_segment); for(j in hit){id<-id+1L;ans[[id]]<-data.frame(unit=i,variable=p,start=starts[j],end=ends[j],length=rr$lengths[j])} }
    else for(i in seq_len(n)) { rr<-rle(flag[i,]); ends<-cumsum(rr$lengths); starts<-ends-rr$lengths+1L; hit<-which(rr$values&rr$lengths>=min_segment); for(j in hit){id<-id+1L;ans[[id]]<-data.frame(unit=i,start=starts[j],end=ends[j],length=rr$lengths[j])} }
    if(!length(ans)) data.frame() else do.call(rbind,ans)
  }
  runs_to_mask <- function(runs,level) {
    if(level=="time") {
      mask<-matrix(FALSE,n,tt)
      if(nrow(runs)) for(j in seq_len(nrow(runs))) mask[runs$unit[j],runs$start[j]:runs$end[j]]<-TRUE
    } else {
      mask<-array(FALSE,c(n,tt,pp))
      if(nrow(runs)) for(j in seq_len(nrow(runs))) mask[runs$unit[j],runs$start[j]:runs$end[j],runs$variable[j]]<-TRUE
    }
    mask
  }

  # V7 exclusive hierarchy. Raw flags remain available for diagnostics.
  observation_outliers<-which(!best$keep)
  global_outliers<-which(best$global_weights==0,arr.ind=TRUE)

  eligible_time<-raw_time_flag
  if(length(observation_outliers)) eligible_time[observation_outliers,]<-FALSE
  multivariate_segments<-get_runs(eligible_time,"time")
  multivariate_segment_mask<-runs_to_mask(multivariate_segments,"time")
  temporal_point_flag<-eligible_time & !multivariate_segment_mask

  eligible_cell<-raw_cell_flag
  if(length(observation_outliers)) eligible_cell[observation_outliers,,]<-FALSE
  if(nrow(global_outliers)) for(j in seq_len(nrow(global_outliers)))
    eligible_cell[,global_outliers[j,1],global_outliers[j,2]]<-FALSE
  time_idx<-which(eligible_time,arr.ind=TRUE)
  if(nrow(time_idx)) for(j in seq_len(nrow(time_idx)))
    eligible_cell[time_idx[j,1],time_idx[j,2],]<-FALSE
  sensor_segments<-get_runs(eligible_cell,"cell")
  sensor_segment_mask<-runs_to_mask(sensor_segments,"cell")
  isolated_cell_flag<-eligible_cell & !sensor_segment_mask

  best$cluster<-ifelse(best$keep,best$z,0L); best$nearest_cluster<-best$z
  best$observation_outliers<-observation_outliers
  best$global_sensor_time_outliers<-global_outliers
  best$time_point_outliers<-which(temporal_point_flag,arr.ind=TRUE)
  best$multivariate_segments<-multivariate_segments
  best$sensor_segments<-sensor_segments
  best$isolated_cell_outliers<-which(isolated_cell_flag,arr.ind=TRUE)
  best$cell_outliers<-best$isolated_cell_outliers
  best$raw_flags<-list(
    cell_outliers=which(raw_cell_flag,arr.ind=TRUE),
    time_point_outliers=which(raw_time_flag,arr.ind=TRUE),
    sensor_segments=get_runs(raw_cell_flag,"cell"),
    multivariate_segments=get_runs(raw_time_flag,"time"))
  best$exclusive_masks<-list(temporal_point=temporal_point_flag,
    multivariate_segment=multivariate_segment_mask,sensor_segment=sensor_segment_mask,
    isolated_cell=isolated_cell_flag)
  best$initial_outliers<-which(!init_keep)
  best$initial_outlyingness<-out_init$score; best$standardization<-scaled[c("center","scale")]
  best$settings<-list(k=k,alpha=alpha,latent_rank=latent_rank,spar=spar,c_cell=c_cell,cell_threshold=cell_threshold,
    time_fraction_cut=time_fraction_cut,global_fraction_cut=global_fraction_cut,local_trim=local_trim,min_cluster=min_cluster,
    method=method,cutoff_prob=cutoff_prob,cutoff_bounds=cutoff_bounds,scale_floor=scale_floor,
    repair=repair,rollback_nonconverged=rollback_nonconverged,trace_history=trace_history,
    selection_rule=selection_rule,cycle_detection=cycle_detection,cycle_max_period=cycle_max_period,cycle_tol=cycle_tol,seed=seed)
  if(return_candidates) best$candidates<-candidates
  class(best)<-"mtrclust_v12"; best
}

print.mtrclust_v12 <- function(x,...) {
  cat("MTRClust prototype v12 [",x$settings$method,"]\nRetained cluster sizes:",paste(x$retained_cluster_sizes,collapse=", "),"\n",sep="")
  cat("Observation outliers:",length(x$observation_outliers),"\nIsolated cell outliers:",nrow(x$isolated_cell_outliers),"\n")
  cat("Temporal point outliers:",nrow(x$time_point_outliers),"\nGlobal sensor-time outliers:",nrow(x$global_sensor_time_outliers),"\n")
  cat("Iterations:",x$iterations," Converged:",x$converged,"\n"); invisible(x)
}
