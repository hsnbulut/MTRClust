robust_residual_scale_v12 <- function(v,floor=.05) {
  med<-median(v,na.rm=TRUE);s<-1.4826*median(abs(v-med),na.rm=TRUE)
  if(!is.finite(s)||s<floor)s<-IQR(v,na.rm=TRUE)/1.349
  if(!is.finite(s)||s<floor)s<-1
  s
}

apply_diagnostics_v12 <- function(fit,method=c("v6_fixed","pooled_all","hierarchical_dual"),min_segment=3L) {
  method<-match.arg(method);out<-fit;out$diagnostic_method<-method
  if(method=="v6_fixed") {
    out$diagnostic_scales<-list(pooled=rep(1,dim(fit$residuals)[3]),applied=NULL)
    class(out)<-c("mtrclust_v12_diagnostics",class(out));return(out)
  }
  r<-fit$residuals;n<-dim(r)[1];tt<-dim(r)[2];pp<-dim(r)[3];k<-dim(fit$centers)[1]
  keep<-fit$keep;z<-fit$nearest_cluster;c_cell<-fit$settings$c_cell;cell_threshold<-fit$settings$cell_threshold
  time_fraction_cut<-fit$settings$time_fraction_cut;global_fraction_cut<-fit$settings$global_fraction_cut
  pooled<-vapply(seq_len(pp),function(p)robust_residual_scale_v12(as.vector(r[keep,,p,drop=FALSE])),numeric(1))
  applied<-matrix(rep(pooled,each=k),k,pp)
  u<-array(0,dim(r));for(i in seq_len(n))u[i,,]<-sweep(r[i,,],2,applied[z[i],],"/")
  cutoff<-c_cell*sqrt(1-sqrt(cell_threshold));raw_cell_flag<-abs(u)>cutoff
  v6_cell_flag<-array(FALSE,dim(r));v6_cell_flag[fit$raw_flags$cell_outliers]<-TRUE
  raw_time_score<-apply(raw_cell_flag,c(1,2),mean);raw_time_flag<-raw_time_score>=time_fraction_cut
  global_fraction<-apply(raw_cell_flag[keep,,,drop=FALSE],c(2,3),mean);global_flag<-global_fraction>=global_fraction_cut

  get_runs<-function(flag,level){ans<-list();id<-0L
    if(level=="cell")for(i in seq_len(n))for(p in seq_len(pp)){rr<-rle(flag[i,,p]);ends<-cumsum(rr$lengths);starts<-ends-rr$lengths+1L;hit<-which(rr$values&rr$lengths>=min_segment);for(j in hit){id<-id+1L;ans[[id]]<-data.frame(unit=i,variable=p,start=starts[j],end=ends[j],length=rr$lengths[j])}}
    else for(i in seq_len(n)){rr<-rle(flag[i,]);ends<-cumsum(rr$lengths);starts<-ends-rr$lengths+1L;hit<-which(rr$values&rr$lengths>=min_segment);for(j in hit){id<-id+1L;ans[[id]]<-data.frame(unit=i,start=starts[j],end=ends[j],length=rr$lengths[j])}}
    if(!length(ans))data.frame()else do.call(rbind,ans)}
  runs_to_mask<-function(runs,level){if(level=="time"){mask<-matrix(FALSE,n,tt);if(nrow(runs))for(j in seq_len(nrow(runs)))mask[runs$unit[j],runs$start[j]:runs$end[j]]<-TRUE
    }else{mask<-array(FALSE,c(n,tt,pp));if(nrow(runs))for(j in seq_len(nrow(runs)))mask[runs$unit[j],runs$start[j]:runs$end[j],runs$variable[j]]<-TRUE};mask}
  observation_outliers<-which(!keep);global_outliers<-which(global_flag,arr.ind=TRUE)
  eligible_time<-raw_time_flag;if(length(observation_outliers))eligible_time[observation_outliers,]<-FALSE
  multivariate_segments<-get_runs(eligible_time,"time");multivariate_segment_mask<-runs_to_mask(multivariate_segments,"time")
  temporal_point_flag<-eligible_time&!multivariate_segment_mask
  eligible_group_cell<-raw_cell_flag;if(length(observation_outliers))eligible_group_cell[observation_outliers,,]<-FALSE
  if(nrow(global_outliers))for(j in seq_len(nrow(global_outliers)))eligible_group_cell[,global_outliers[j,1],global_outliers[j,2]]<-FALSE
  time_idx<-which(eligible_time,arr.ind=TRUE);if(nrow(time_idx))for(j in seq_len(nrow(time_idx)))eligible_group_cell[time_idx[j,1],time_idx[j,2],]<-FALSE
  sensor_segments<-get_runs(eligible_group_cell,"cell");sensor_segment_mask<-runs_to_mask(sensor_segments,"cell")
  isolated_base<-if(method=="hierarchical_dual")v6_cell_flag else raw_cell_flag
  eligible_isolated<-isolated_base;if(length(observation_outliers))eligible_isolated[observation_outliers,,]<-FALSE
  if(nrow(global_outliers))for(j in seq_len(nrow(global_outliers)))eligible_isolated[,global_outliers[j,1],global_outliers[j,2]]<-FALSE
  if(nrow(time_idx))for(j in seq_len(nrow(time_idx)))eligible_isolated[time_idx[j,1],time_idx[j,2],]<-FALSE
  isolated_cell_flag<-eligible_isolated&!sensor_segment_mask
  out$cell_weights<-tukey_weight(u,c_cell);out$cell_flags<-raw_cell_flag;out$time_scores<-raw_time_score
  tw<-1-raw_time_score/time_fraction_cut;tw[tw<0]<-0;out$time_weights<-tw
  out$global_fraction<-global_fraction;out$global_weights<-ifelse(global_flag,0,1)
  out$global_sensor_time_outliers<-global_outliers;out$time_point_outliers<-which(temporal_point_flag,arr.ind=TRUE)
  out$multivariate_segments<-multivariate_segments;out$sensor_segments<-sensor_segments
  out$isolated_cell_outliers<-which(isolated_cell_flag,arr.ind=TRUE);out$cell_outliers<-out$isolated_cell_outliers
  out$raw_flags<-list(cell_outliers=which(raw_cell_flag,arr.ind=TRUE),time_point_outliers=which(raw_time_flag,arr.ind=TRUE),
    sensor_segments=get_runs(raw_cell_flag,"cell"),multivariate_segments=get_runs(raw_time_flag,"time"))
  out$exclusive_masks<-list(temporal_point=temporal_point_flag,multivariate_segment=multivariate_segment_mask,
    sensor_segment=sensor_segment_mask,isolated_cell=isolated_cell_flag)
  out$diagnostic_scales<-list(pooled=pooled,applied=applied,cutoff=cutoff)
  class(out)<-c("mtrclust_v12_diagnostics",class(out));out
}

print.mtrclust_v12_diagnostics<-function(x,...){cat("MTRClust prototype v12 diagnostics [",x$diagnostic_method,"]\n",sep="");NextMethod("print")}
