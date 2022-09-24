#a function call in code #3
adjust_dws<-function(acause,effect_size_draws=NULL,effect_size=NULL,effect_size_sd=NULL, 
                     optimal_e=NULL, optimal_se=NULL,return_draws=F,collapse_idividuals=T){
  # set seed
  set.seed(555)
  # read in draw files
  dw_draws_dir<-paste0("PATH")
  message(paste("Reading in draws from", dw_draws_dir))
  dw_dt<-fread(paste0(dw_draws_dir,"t",acause,".csv"))
  
  drop<-c(paste0("DW_diff_pred_",0:999),paste0("DW_pred_",0:999),paste0("DW_counter_",0:999))
  dw_dt<-dw_dt[,c(names(dw_dt)[!(names(dw_dt)%in%drop)]),with=F]
  
  dw_dt[,row_num:=.I]
  
  # melt draws long
  id_vars<-names(dw_dt)[grepl("\\d",names(dw_dt))==F]
  dw_dt<-melt(dw_dt,id.vars =id_vars,value.name = "dw")
  
  # Pull out draw_num from variable col
  dw_dt[,draw:= paste0("draw_",as.numeric(gsub("\\D","",variable)))]
  dw_dt[,variable:=gsub("_\\d*","",variable)]
  dw_dt[,dw_t := dw]
  
  #readin individual level
  ind <- fread("FILEPATH/adjust_beta_all_panels_long.csv")[cond_name==paste0("t",acause), -c("physio", "physio_psych","surg","opioids","psych","epidural_steroids_inject",
                                                                                                                                            "pharmaNonOp","physio_pharmaNonOp","physio_psych_pharmaNonOp","beta","pro_beta_n")]
  ind<- unique(ind)
  dw_dt <- merge(dw_dt, ind, by=c('id'))
  
  
  # Read in SF-12 to DW map 
  message("Reading in SF-12 to DW map from MR BERT object")
  remap_dw <- function(x){
    draw_specific_lbp_dws <- dw_dt[draw == x, ]
    s<-sort(dw_map[draw == x,dw_t])
    draw_specific_lbp_dws[,closest_dw:=s[findInterval(dw_t,s,all.inside = T)]]
    draw_specific_lbp_dws<-merge(draw_specific_lbp_dws,dw_map[draw == x,.(sf = sf12, dw_t)],by.x="closest_dw",by.y="dw_t")
  }
  
  lbp_dws_list <- lapply(unique(dw_dt$draw), remap_dw)
  lbp_dws <- Reduce(function(x, y){rbind(x, y)}, lbp_dws_list)
  rm(lbp_dws_list) # to save memory
  
  # calculate the SD of sf for each draw
  lbp_dws[,sf_sd:=sd(sf),by=c("draw")]
  
  if (is.null(effect_size_draws)){
    # create draws of intervention effect
    es_draws<-rnorm(1000,effect_size,effect_size_sd)
    es_dt<-data.table(draw_num=0:999,es=es_draws)
    
    # create draws of optimal effect
    os_draws<-rnorm(1000,optimal_e,optimal_se)
    os_dt<-data.table(draw_num=0:999,os=os_draws) 
    
    # merge 
    es_dt <- merge(es_dt, os_dt, by="draw_num")
  }else{
    es_dt<-copy(effect_size_draws)
  }
  names(lbp_dws)[names(lbp_dws) == "draw"] <- "draw_num"
  # apply the effect size
  dw_dt<-merge(lbp_dws[,-c('cond_name')],es_dt,by=c("draw_num"))
  
  #no treatment
  dw_dt[sum==0,sf_adj:= sf]
  dw_dt[sum>=1,sf_adj:= sf + sf_sd*adjust_es]
  dw_dt[,sf_os_adj := sf_adj - sf_sd*os]
  
  # Map back from SF to DW
  remap_sf <- function(x){
    draw_specific_lbp_dws <- dw_dt[draw_num == x, ]
    s<-sort(dw_map[draw == x,sf12])
    ## for adj_dw
    draw_specific_lbp_dws[,closest_sf:=s[findInterval(sf_adj,s,all.inside = T)]]
    draw_specific_lbp_dws<-merge(draw_specific_lbp_dws,dw_map[draw == x,.(sf = sf12, dw_adj = dw_t)],by.x="closest_sf",by.y="sf")
    ## for FUOT
    draw_specific_lbp_dws[,closest_sf:=s[findInterval(sf_os_adj,s,all.inside = T)]]
    draw_specific_lbp_dws<-merge(draw_specific_lbp_dws,dw_map[draw == x,.(sf = sf12, dw_adj_os = dw_t)],by.x="closest_sf",by.y="sf")
  }
  
  lbp_dws_list <- lapply(unique(lbp_dws$draw_num), remap_sf)
  message("map back from SF to DW")
  lbp_dws <- Reduce(function(x, y){rbind(x, y)}, lbp_dws_list)
  rm(lbp_dws_list) # to save memory
  
  
  lbp_dws[sf_adj > 120, dw_adj := 0]
  lbp_dws[sf_os_adj > 120, dw_adj_os:= 0]
  
  # Clean up
  dw_dt<-lbp_dws[,.(id,draw_num,dw,dw_adj, dw_adj_os)]
  
  # Melt pre and post adjustment DWs 
  dw_dt<-melt(dw_dt,id.vars = c("id","draw_num"))
  
  # Cap values below 0 or above 1
  dw_dt[value>1,value:=1]
  dw_dt[value<0,value:=0]
  
  if (collapse_idividuals){
    dw_dt<-dw_dt[,mean_dw:=mean(value),by=c("draw_num","variable")]
    
    if (!return_draws){
      # calculate mean of draws 
      dw_dt<-dw_dt[,.(avg_dw=mean(mean_dw),lower_avg_dw=quantile(mean_dw,0.025),upper_avg_dw=quantile(mean_dw,0.975)),
                   by="variable"]
    }
  }
  
  return(dw_dt)
}
