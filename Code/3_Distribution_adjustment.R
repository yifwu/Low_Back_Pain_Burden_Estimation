#####################################
## Docstring ----------------------------------------------------------- ----
## Project: GBD LBP Severity Distribution varying by HAQI Paper
## Details: Step 3:Create Intermediate file/calculation, distribution adjustment,
## --------------------------------------------------------------------- ----

## Environment Prep ---------------------------------------------------- ----

#Align with the write up
# clear workspace
rm(list=ls())

version_id <- 3
pacman::p_load(msm, data.table, readstata13, ggplot2, openxlsx, scales, gridExtra)

library(mrbrt001, lib.loc = '/ihme/code/mscm/R/packages/')
# source adjustment function
# it is sourcing an DW adjustment function
source("FILEPATH")
source(file.path(central_lib,"get_location_metadata.R"))
source(file.path(central_lib,"get_draws.R"))
source(file.path(central_lib,"get_population.R"))
source(file.path(central_lib,"get_covariate_estimates.R"))
source("FILEPATH","mr_brt_functions.R")
rlogit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}

set.seed(1640)

data_dir<-paste0("FILEPATH")
input_dir<-paste0("FILEPATH")

map_dir<-paste0(input_dir,"maps/")
output_dir<-paste0("FILEPATH")
plot_dir <- paste0("FILEPATH")

remake <-FALSE
# weight effect size by coverage ------------------------------------------

# pull in coverage estimates
if (remake==TRUE){
  cov_ests<-fread(paste0(output_dir,"utilization_by_condition_with_sd.csv"))
  cov_ests <- cov_ests[!(cov_ests$variable=="analgesic" | cov_ests$variable=="nsaid"),]
  # add scenario var
  cov_ests[,scenario:="current"]
  
  # cast wide by intervention
  cov_ests[,t:=1]
  cov_ests<-dcast(cov_ests,icd9codx+cond_name+scenario+beta+n+value+value_sd+lower+upper~variable,value.var = "t")
  cols = c("epidural_steroids_inject","opioids","pharmaNonOp","physio","physio_pharmaNonOp", "physio_psych",
           "physio_psych_pharmaNonOp","psych","surg" )
  cov_ests[ , (cols) := lapply(.SD, nafill, fill=0), .SDcols = cols]
  
  # pull in opitmal coverage df
  cov_ests_op<-fread(paste0(output_dir,"optimals_coverage_ests_test.csv"))
  
  # bind 
  cov_ests<-rbindlist(list(cov_ests,cov_ests_op),use.names = T,fill=T)
  
  # create draws of covariate estimates
  coverage_draws<-list()
  for (i in 1:nrow(cov_ests)){
    cov_ids<-cov_ests[i,.(scenario,epidural_steroids_inject,surg,opioids,pharmaNonOp,physio,physio_pharmaNonOp,
                          physio_psych,physio_psych_pharmaNonOp,psych,
                          icd9codx,cond_name,beta,n,value,value_sd)]
    temp<-data.table(cov_ids,draws=rnorm(n=1000,mean=cov_ests[i,value],sd = cov_ests[i,value_sd]))
    temp[,draw_num:=.I-1]
    temp[,draw_num:=paste0("draw_",draw_num)]
    coverage_draws[[i]]<-dcast(temp,scenario+epidural_steroids_inject+surg+opioids+pharmaNonOp+physio+physio_pharmaNonOp+
                               physio_psych+physio_psych_pharmaNonOp+psych+icd9codx+cond_name+beta+n+value+value_sd~draw_num,value.var = "draws")
    merge_vars<-names(cov_ids)
  }
  coverage_draws<-rbindlist(coverage_draws,use.names = T)
  coverage_draws<-melt(coverage_draws,id.vars =merge_vars, value.name = "utilization",variable.name = "draw")
  
  # merge coverage draws onto cov_ests
  cov_ests<-merge(cov_ests,coverage_draws,by=merge_vars)
  # weight overall es by coverage
  #get the coefficient from our Mr BERT analysis 
  
  # calculate weighted effect size
  #made a change here if optimal, the utilization for optimal treatment should be 1
  #0.42 came from for those lbp with leg pain, 0.42 are eligible for surgery
  #cov_ests[(scenario=="Optimal" & cond_name=="tmsk_pain_lowback_wleg"), utilization := 0.42]
  #cov_ests[(scenario=="Optimal" & cond_name!="tmsk_pain_lowback_wleg"), utilization := 1]
  #then multiple them together 
  #ended up changing the input file optimal file
  cov_ests[,weighted_eff:=utilization*beta]
  
  # sum by condition and scenario
  cov_ests<-cov_ests[,.(weighted_eff=sum(weighted_eff)),by=c("cond_name","scenario","draw")]
  
  #all of sudden at this point, no uncertainty, all the uncertainty disappeared 
  # get the mean,sd,lower,and,upper for reporting purposes 
  cov_ests_collapsed<-copy(cov_ests[,.(weighted_eff=mean(weighted_eff),weighted_eff_sd=sd(weighted_eff),
                                       weighted_eff_lower=quantile(weighted_eff,0.025),weighted_eff_upper=quantile(weighted_eff,0.975)),
                                    by=c("cond_name","scenario")])
  cov_ests_collapsed<-dcast(cov_ests_collapsed,cond_name~scenario,
                            value.var = c("weighted_eff","weighted_eff_sd","weighted_eff_lower","weighted_eff_upper"))
  write.csv(cov_ests_collapsed,paste0(output_dir,"pooled_effect_estimate_summaries_v3.csv"),row.names=F)
  
  # make draw numeric
  #cov_ests[,draw:=as.numeric(gsub(pattern = "\\D*","",draw))]
  
  # cast wide by scenario
  cov_ests<-dcast(cov_ests,cond_name+draw~scenario,
                  value.var = c("weighted_eff"))
  setnames(cov_ests,c("current","Optimal"),c("es","os"))
  setnames(cov_ests,"draw","draw_num")
  #   -----------------------------------------------------------------------
  
  # adust LBP DWs -----------------------------------------------------------
  
  # generate map of SF12 to DW values - takes a long time, only run this if the function
  # to map needs to be updated
  
  ## conduct standard SF-12 to DW model and estimate an sf-12 score for each person attributable to anxiety disorders
  message("Reading in SF-12 to DW map from MR BERT object")
  dw_map_model <- py_load_object(filename = "/home/j/temp/Holly/damians/publications/severity_by_haqi/sf12_to_dw.pkl", pickle = "dill")
  
  dw_map <- data.table(intercept = 1, sf12_c = seq(40, 120, 0.01)-120)
  predict_data <- MRData()
  predict_data$load_df(data = dw_map, col_covs=as.list(dw_map_model$cov_names))
  
  beta_samples <- mrbrt001::core$other_sampling$sample_simple_lme_beta(1000L, dw_map_model)
  gamma_outer_samples <- matrix(rep(dw_map_model$gamma_soln, each = 1000L), nrow = 1000L)
  draws <- dw_map_model$create_draws(predict_data, beta_samples = beta_samples, gamma_samples = gamma_outer_samples, random_study = T)
  
  dw_map$dw_logit_mean <- apply(draws, 1, function(x) mean(x))
  dw_map$dw_logit_lo <- apply(draws, 1, function(x) quantile(x, 0.025))
  dw_map$dw_logit_hi <- apply(draws, 1, function(x) quantile(x, 0.975))
  dw_map[, `:=` (sf12 = sf12_c + 120, dw_logit_se = (dw_logit_hi - dw_logit_lo)/3.92)]
  
  draws <- data.table(draws)
  names(draws) <- paste0("draw_", as.numeric(gsub("V", "", names(draws)))-1)
  dw_map <- cbind(dw_map, draws)
  dw_map <- melt.data.table(dw_map, id.vars = names(dw_map)[!(names(dw_map) %like% "draw")], value.name="dw_logit", variable.name="draw")
  dw_map[, dw_t := rlogit(dw_logit)]
  
  # adjust dws for w/ and w/o leg pain separately
  # this step takes forever to run
  # getting sf-12 scores and convert the disability weight back to sf-12 scores. Shift the 
  lbp_noleg<-adjust_dws(acause="msk_pain_lowback_noleg",
                        effect_size_draws = cov_ests[cond_name=="tmsk_pain_lowback_noleg"],
                        collapse_idividuals = F)
  
  lbp_wleg<-adjust_dws(acause="msk_pain_lowback_wleg",
                       effect_size_draws = cov_ests[cond_name=="tmsk_pain_lowback_wleg"],
                       collapse_idividuals = F)
  
  # combine datasets
  lbp<-rbindlist(list(lbp_noleg[,condition:="LBP - no leg involvement"],lbp_wleg[,condition:="LBP - leg involvement"]),use.names=T)
  write.csv(lbp, paste0(output_dir,"lbp_intermediate_",version_id,".csv"), row.names = F)
} else {
  lbp<-fread(paste0(output_dir,"lbp_intermediate_",2,".csv"))
 }

#test<-lbp[id==916290108]
#lbp <- lbp_noleg[,condition:="LBP - no leg involvement"]
  
# plot 
plot_dt<-copy(lbp)
lbp_int <- copy(lbp)

# average each individuals dw across the draws
plot_dt<-plot_dt[,.(mean_dw=mean(value)),by=c("id","variable","condition")]
# average across individuals
plot_dt[,avg_dw:=mean(mean_dw),by=c("variable","condition")]
# label the variable
plot_dt[,variable:=factor(variable,levels=c("dw_adj","dw","dw_adj_os"),labels = c("With no treatment","With current treatment",
                                                                                  "Full utilization of optimal treatment"))]

#I tried the overlapping histagram too, it does not look good. It is still better to separated them out.
#ggnoleg<-ggnoleg+coord_cartesian(xlim=c(0,.25))
#gg<-grid.arrange(ggleg,ggnoleg,ncol=2)
#ggleg<-ggplot(data=plot_dt[condition=="LBP - leg involvement"])+
#  geom_density(aes(x=mean_dw, colour=variable, fill=variable),alpha=0.5)+
#  scale_x_continuous(expand = c(0,0))+
#  scale_y_continuous(expand = c(0,0))
#  geom_vline(aes(xintercept=avg_dw),size=1)
#ggleg<-ggleg+coord_cartesian(xlim=c(0,.25))


# generate density plots
#partially displayed problem
#scaleFUN <- function(x) sprintf("%.1f", x)
ggleg<-ggplot(data=plot_dt[condition=="LBP - leg involvement"])+
  geom_density(aes(x=mean_dw),alpha=0.4,fill="#9999FF")+
  geom_vline(aes(xintercept=avg_dw),size=1)+
  facet_wrap(~variable,scales="free_y",ncol = 1,nrow = 3)+
  theme_test()+
  xlab("LBP Disability Weight")+
  ylab("Density")+
  ggtitle("LBP - leg involvement")+
  #scale_x_continuous(expand = c(0,0),limits = c(0,1),labels=scaleFUN)+
  scale_x_continuous(expand = c(0,0),limits = c(0,1))+
  scale_y_continuous(expand = c(0,0))+
  theme(strip.background = element_blank(),strip.text.x = element_text(size=10,face = "italic"),
        plot.title = element_text(hjust=0.5,face="bold"))
ggnoleg<-ggplot(data=plot_dt[condition=="LBP - no leg involvement"])+
  geom_density(aes(x=mean_dw),alpha=0.4,fill="#CC9900")+
  geom_vline(aes(xintercept=avg_dw),size=1)+
  facet_wrap(~variable,scales="free_y",ncol = 1,nrow = 3)+
  theme_test()+
  xlab("LBP Disability Weight")+
  ylab("Density")+
  ggtitle("LBP - no leg involvement")+
  scale_x_continuous(expand = c(0,0))+
  #scale_x_continuous(expand = c(0,0), labels=scaleFUN)+
  scale_y_continuous(expand = c(0,0))+
  theme(strip.background = element_blank(),strip.text.x = element_text(size=10,face = "italic"),
        plot.title = element_text(hjust=0.5,face="bold"))
gg<-grid.arrange(ggleg,ggnoleg,ncol=2)
ggsave(filename = paste0(plot_dir,"avg_dw_density_plot_zoom_out_",version_id,".png"),plot = gg,
       height=6,width=10)


#scaleFUN <- function(x) sprintf("%.2f", x)
ggleg<-ggleg+coord_cartesian(xlim=c(0,0.3))
ggnoleg<-ggnoleg+coord_cartesian(xlim=c(0,0.3))
gg<-grid.arrange(ggleg,ggnoleg,ncol=2)
ggsave(filename = paste0(plot_dir,"avg_dw_density_plot_",version_id,".png"),plot = gg,
       height=6,width=10)
#ggsave(filename = paste0(plot_dir,"avg_dw_density_plot_setting_opioid_zero.png"),plot = gg,
#       height=6,width=10)

# Average across everyone - keep one dt with avg DW by condition for reporting purposes
lbp_by_cond<-copy(lbp[,.(mean_dw=mean(value)),by=c("draw_num","variable","condition")])
lbp<-lbp[,.(mean_dw=mean(value)),by=c("draw_num","variable")]
#lbp[,draw_num:=paste0("draw_",draw_num)]

# calculate mean of draws by condition (for reporting) and across whole sample (for regression)
lbp_by_cond<-lbp_by_cond[,.(avg_dw=mean(mean_dw),lower_avg_dw=quantile(mean_dw,0.025),upper_avg_dw=quantile(mean_dw,0.975)),
                         by=c("variable","condition")]
write.csv(lbp_by_cond,paste0(output_dir,"LBP_burden_meps_per_sample_dw_bycond_",version_id,".csv"),row.names=F)
#write.csv(lbp_by_cond,paste0(output_dir,"LBP_burden_meps_per_sample_dw_bycond_setting_opoid_zero.csv"),row.names=F)


lbp_collapsed<-copy(lbp[,.(avg_dw=mean(mean_dw),lower_avg_dw=quantile(mean_dw,0.025),upper_avg_dw=quantile(mean_dw,0.975)),
                        by=c("variable")])
write.csv(lbp_collapsed,paste0(output_dir,"LBP_burden_meps_per_sample_dw_all_",version_id,".csv"),row.names=F)
#write.csv(lbp_collapsed,paste0(output_dir,"LBP_burden_meps_per_sample_dw_all_setting_opioid_zero.csv"),row.names=F)

####ADAPTING THE severity distribution###########
### Conduct severity analysis ###
cause <- c("msk_pain_lowback_noleg","msk_pain_lowback_wleg")

# load severity cutoffs
data <- fread("/home/j/Project/GBD/Systematic Reviews/ANALYSES/MEPS/2018_update/output/3b_meps_severity_cutoffs.csv")

## make 1000 distribution variables
data <- data[yld_cause%in%cause,]
data <- melt.data.table(data, id.vars = names(data)[!(names(data) %like% "MID")], value.name = "MID", variable.name = "draw")
data <- data[draw !="MID_calc",]
data[, `:=` (MID_type = draw)]
data[, `:=` (draw = substring(gsub("MID", "draw_", draw), 1, nchar(gsub("MID", "draw_", draw))-1))]

data <- merge(data[grepl("a", MID_type),.(yld_cause, healthstate_id, severity, draw, MID_lower = MID)],
              data[grepl("b", MID_type),.(yld_cause, healthstate_id, severity, draw, MID_upper = MID)],
              by = c("yld_cause", "healthstate_id", "severity", "draw"))

plot_dt$variable <- factor(plot_dt$variable, levels = c("With current treatment", "With no treatment","Full utilization of optimal treatment"))
plot_dt$variable <- relevel(plot_dt$variable, ref ="With no treatment")
plot_lines <- data.table(severity = c("Asymptomatic", "Mild", "Moderate","Severe"), cutoff = c(-0.02, data[severity == 1, mean(MID_upper)], data[severity == 2, mean(MID_upper)], data[severity==3,mean(MID_upper)]))

plot <- ggplot(plot_dt[condition=="LBP - leg involvement", .(variable, mean_dw = ifelse(mean_dw == 0, -0.02, mean_dw))], aes(x=mean_dw, fill=variable, color =variable)) +
  geom_histogram(position="identity", alpha=0.3)+
  geom_vline(data=plot_lines, aes(xintercept=cutoff), linetype="dashed") + ## If want severity splits
  #scale_x_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-0.02, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0")) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1100, 100))+
  annotate("text", x = plot_lines[severity == "Asymptomatic", cutoff-0.02], y = -12, label = "Asymp.", size = 3) + ## If want severity splits
  annotate("text", x = mean(c(plot_lines[severity == "Asymptomatic", cutoff], plot_lines[severity == "Mild", cutoff])), y = -12, label = "Mild", size = 3) + ## If want severity splits
  annotate("text", x = mean(c(plot_lines[severity == "Mild", cutoff], plot_lines[severity == "Moderate", cutoff])), y = -12, label = "Moderate", size = 3) + ## If want severity splits
  annotate("text", x = mean(c(plot_lines[severity == "Moderate", cutoff], plot_lines[severity == "Severe", cutoff])), y = -12, label = "Severe", size = 3) + ## If want severity splits
  annotate("text", x = mean(c(plot_lines[severity == "Severe", cutoff], 1)), y = -12, label = "Most Severe", size = 3) + ## If want severity splits
  coord_cartesian(ylim = c(0, 1100), clip = "off")+
  ylab("Frequency") +
  xlab("\nDisability Weight") +
  theme(legend.position="top")
plot
ggsave(filename = paste0(plot_dir,"figure_lbp_dw_shift_with_leg_pain_",version_id,".pdf"),plot = plot,
       height=6,width=10)



plot <- ggplot(plot_dt[condition=="LBP - no leg involvement", .(variable, mean_dw = ifelse(mean_dw == 0, -0.02, mean_dw))], aes(x=mean_dw, fill=variable, color =variable)) +
  geom_histogram(position="identity", alpha=0.3)+
  geom_vline(data=plot_lines, aes(xintercept=cutoff), linetype="dashed") + ## If want severity splits
  #scale_x_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-0.02, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0")) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 6600, 1100))+
  annotate("text", x = plot_lines[severity == "Asymptomatic", cutoff-0.02], y = -12, label = "Asymp.", size = 3) + ## If want severity splits
  annotate("text", x = mean(c(plot_lines[severity == "Asymptomatic", cutoff], plot_lines[severity == "Mild", cutoff])), y = -12, label = "Mild", size = 3) + ## If want severity splits
  annotate("text", x = mean(c(plot_lines[severity == "Mild", cutoff], plot_lines[severity == "Moderate", cutoff])), y = -12, label = "Moderate", size = 3) + ## If want severity splits
  annotate("text", x = mean(c(plot_lines[severity == "Moderate", cutoff], plot_lines[severity == "Severe", cutoff])), y = -12, label = "Severe", size = 3) + ## If want severity splits
  annotate("text", x = mean(c(plot_lines[severity == "Severe", cutoff], 1)), y = -12, label = "Most Severe", size = 3) + ## If want severity splits
  coord_cartesian(ylim = c(0, 6600), clip = "off")+
  ylab("Frequency") +
  xlab("\nDisability Weight") +
  theme(legend.position="top")
plot
ggsave(filename = paste0(plot_dir,"figure_lbp_dw_shift_wo_leg_pain_",version_id,".pdf"),plot = plot,
       height=6,width=10)

#rename some columns for merging
#names(data)[names(data) == "yld_cause"] <- "condition"
names(data)[names(data) == "draw"] <- "draw_num"
lbp_int[condition == "LBP - no leg involvement", yld_cause:= "msk_pain_lowback_noleg"]
lbp_int[condition == "LBP - leg involvement", yld_cause:="msk_pain_lowback_wleg"]

#adjustment_factor_level <- "Binned" # Choose Person vs Binned
#found an issue with dcast and change fun.aggregate to mean, it should clearly not be sum
lbp_dws <- merge(lbp_int, data,by=c("draw_num","yld_cause"), allow.cartesian=TRUE)
lbp_dws <-dcast(lbp_dws[,.(draw_num, yld_cause, id, variable, value, condition, healthstate_id, severity, MID_lower,
                        MID_upper)], ...~variable, value.var="value", fun.aggregate=mean)

lbp_dws[severity == 0 , case_severity := ifelse(dw <= 0, 1, 0)]
lbp_dws[severity == 0, case_severity_adj := ifelse(dw_adj <= 0, 1, 0)]
lbp_dws[severity == 0, case_severity_fcot := ifelse(dw_adj_os <= 0, 1, 0)]

lbp_dws[severity != 0, case_severity := ifelse(dw > MID_lower & dw <= MID_upper, 1, 0)]
lbp_dws[severity != 0, case_severity_adj := ifelse(dw_adj > MID_lower & dw_adj <= MID_upper, 1, 0)]
lbp_dws[severity != 0, case_severity_fcot := ifelse(dw_adj_os > MID_lower & dw_adj_os <= MID_upper, 1, 0)]
#test<-lbp_dws[severity==0 &variable=='dw_adj',]
#table(lbp_dws$case_severity, lbp_dws$varaible)

#total cases 69453(54798, 14655)
lbp_dws[, total_cases := max(seq_len(.N)), by = c("draw_num", "severity",'yld_cause')]
lbp_dws[, sev_cases := sum(case_severity), by = c("draw_num", "severity",'yld_cause')]
lbp_dws[, sev_cases_adj := sum(case_severity_adj), by = c("draw_num", "severity",'yld_cause')]
lbp_dws[, sev_cases_fcot := sum(case_severity_fcot), by = c("draw_num", "severity",'yld_cause')]
lbp_dws[, sev_prop := sev_cases / total_cases]
lbp_dws[, sev_prop_adj := sev_cases_adj / total_cases]
lbp_dws[, sev_prop_fcot := sev_cases_fcot / total_cases]
lbp_dws <- unique(lbp_dws[,.(yld_cause, healthstate_id, severity, sev_prop, sev_prop_adj, sev_prop_fcot, draw_num)])

#read in the weights
weights <- fread("/home/j/WORK/04_epi/03_outputs/01_code/02_dw/02_standard/dw_full.csv")
weights <- melt.data.table(weights, id.vars = c("hhseqid", "healthstate_id", "healthstate"), variable.name = "draw", value.name = "dw")
weights[, draw_num := gsub("draw", "draw_", draw)]


lbp_dws <- merge(lbp_dws, weights, by = c("healthstate_id", "draw_num"))
unique(lbp_dws$healthstate)
lbp_dws[severity == 0, dw := 0]
lbp_dws[, `:=` (dw_mean = sev_prop*dw, dw_mean_adj = sev_prop_adj * dw, dw_mean_fcot = sev_prop_fcot * dw)]
lbp_dws[, `:=` (dw_mean = sum(dw_mean), dw_mean_adj = sum(dw_mean_adj), dw_mean_fcot = sum(dw_mean_fcot)), by =c("draw_num",'yld_cause')]

severity_table <- rbind(lbp_dws[,.(Scenario = "Observed", mean= mean(sev_prop), lower = quantile(sev_prop, 0.025), upper = quantile(sev_prop, 0.975)), by = c('severity','yld_cause')],
                        lbp_dws[,.(Scenario = "No treatment", mean= mean(sev_prop_adj), lower = quantile(sev_prop_adj, 0.025), upper = quantile(sev_prop_adj, 0.975)), by = c('severity','yld_cause')],
                        lbp_dws[,.(Scenario = "FUOT", mean= mean(sev_prop_fcot), lower = quantile(sev_prop_fcot, 0.025), upper = quantile(sev_prop_fcot, 0.975)), by = c('severity','yld_cause')])

severity_table[, Sequela := ifelse(severity == 0, "Asymptomatic", ifelse(severity == 1, "Mild", ifelse(severity == 2, "Moderate", ifelse(severity == 3,"Severe", "Most Severe"))))]
severity_table[, `:=` (Proportion = paste0(round(mean*100, 1), " (", round(lower*100, 1), " to ", round(upper*100, 1), ")"))]
severity_table
severity_table <- dcast(severity_table[,.(Scenario, Sequela, yld_cause, Proportion)], ...~Scenario, value.var="Proportion")

write.csv(severity_table, paste0(output_dir,"severity_table_USA_",version_id,".csv", row.names = F))
write.csv(lbp_dws, paste0(output_dir,"lbp_dws_intermediate_",version_id,".csv"),row.names = F )
