## Docstring ----------------------------------------------------------- ----
## Project: GBD LBP Severity Distribution varying by HAQI Paper
## Script: 
## Description: Prepare for Mr Bert Run
## --------------------------------------------------------------------- ----

## Environment Prep ---------------------------------------------------- ----
rm(list=ls(all.names=T))

#set input
treatment_dir<-paste0("FILEPATH")
table_dir<-paste0(treatment_dir,"FILEPATH")

library(crosswalk, lib.loc = "FILEPATH")
pacman::p_load(dplyr, gtools, readxl, msm, data.table,metafor,ggplot2 )

meta <- read_excel("FILEPATH of the INPUT DATA")
meta <- as.data.table(meta)


#remove those lines that contains information about SF-12 and SF-36. 
#those should not be counted in the disability measurement
meta <- meta[meta$`Effect Size Type` %in% c('SMD','MD'),]

# make it consistent so that negative effect is improvement
effect_cols<-c("Intervention Mean","Comparison Mean")
meta[`Improvement Direction`==1,(effect_cols):=lapply(.SD,"*",-1),.SDcols=effect_cols]


dt_meta <- as.data.table(
  escalc(measure = "SMD",
         m1i = `Intervention Mean`,
         sd1i = `Intervention SD`,
         n1i = `Intervention N`,
         m2i = `Comparison Mean`,
         sd2i = `Comparison SD`,
         n2i = `Comparison N`,
         data=meta)
)

#the output yi:vector to specify the observed effect size or outcomes.
#the output vi: vector to specify the corresponding sampling variances.
#there are also rows that contains calculated or extracted SF-36 and SF-12, need to exclude them.
#only keeping MD and SMD
meta$logit_diff<- dt_meta$yi
meta$logit_diff_se <- sqrt(dt_meta$vi)

dt <- meta[!is.na(meta$logit_diff),]
dt$Intervention_class_flag <- tolower(dt$Intervention_class_flag)
dt$Comparison_flag <- tolower(dt$Comparison_flag)

colnames(dt)[colnames(dt) == 'Intervention_class_flag'] <- 'dorm_alt'
colnames(dt)[colnames(dt) == 'Comparison_flag'] <- 'dorm_ref'

dt_1 <- dt[,c('dorm_alt','dorm_ref','logit_diff','logit_diff_se','Study Name',"Effect Size Mean","follow_up_weeks",'Instrument/Scale')]

dt_3<-dt_1[!is.na(dorm_ref)& !is.na(dorm_alt)& !is.na(logit_diff),]

dt_3$study <- dt_3$`Study Name`
#physio and psychological are the same categories
dt_3$dorm_alt <-ifelse(dt_3$dorm_alt %in% c('active physio + passive physio','active + education', 
                                            'active + education + psychological', "active physio", 
                                            "passive physio","active + passive + education","passive + education"), 'physio',
                       ifelse(dt_3$dorm_alt %in% c("active + education + pharma_non_op","active + passive + education +  pharma_non_op",
                                                   "passive physio + pharma_non_op","pharma_non_op + active physio"),'physio_pharmaNonOp',
                              ifelse(dt_3$dorm_alt %in% c('active + passive + education + psychological','active + passive + psychological',
                                                          'passive physio + psychological','active + psychological'),"physio_psych",
                                     ifelse(dt_3$dorm_alt %in% c('active + education + psychological + pharma_non_op',
                                                                 'active + passive + education + psychological + pharma_non_op'),'physio_psych_pharmaNonOp',
                                            ifelse(dt_3$dorm_alt %in% c('psychological + education'),'psychological',dt_3$dorm_alt)))))


dt_3$dorm_ref <-ifelse(dt_3$dorm_ref %in% c('active physio + passive physio','active + education', 
                                            'active + education + psychological', "active physio", 
                                            "passive physio","active + passive + education","passive + education", 'placebo + active physio','passive physio + placebo','care as usual + passive physio'), 'physio',
                       ifelse(dt_3$dorm_ref %in% c("active + education + pharma_non_op","active + passive + education +  pharma_non_op",
                                                   "passive physio + pharma_non_op","pharma_non_op + active physio"),'physio_pharmaNonOp',
                              ifelse(dt_3$dorm_ref %in% c('active + passive + education + psychological','active + passive + psychological',
                                                          'passive physio + psychological','active + psychological','active physio + passive physio + psychological'),"physio_psych",
                                     ifelse(dt_3$dorm_ref %in% c('active + education + psychological + pharma_non_op',
                                                                 'active + passive + education + psychological + pharma_non_op'),'physio_psych_pharmaNonOp',
                                            ifelse(dt_3$dorm_ref %in% c('psychological + education'),'psychological',
                                                   ifelse(dt_3$dorm_ref %in% c('care as usual','placebo','care as usual + education','wait list control'),'control',dt_3$dorm_ref))))))

#what about those groups with placebo
dt_check<-dt_3[(grepl('+',dt_3$dorm_alt, fixed = TRUE) | grepl('+',dt_3$dorm_ref, fixed=TRUE))=='FALSE', ]

#n=118 had the within group comparison
dt_final <- dt_check[dt_check$dorm_alt!=dt_check$dorm_ref,]
#n=188, 188 studies were kept
dt_final$dorm_ref <- ifelse(dt_final$dorm_ref=='wait list control', 'placebo', dt_final$dorm_ref)


#only limit to 3 months to 12 months
#dropped 31 observations
dt_final_1 <- dt_final[dt_final$follow_up_weeks <=52 | is.na(dt_final$follow_up_weeks), ]

#MR_BRT Model
dat4 <- CWData(
  df = dt_final_1,    # dataset for metaregression
  obs = "logit_diff",  # column name for the observation mean
  obs_se = "logit_diff_se", # column name for the observation standard error
  alt_dorms = "dorm_alt", # variable with alternative def/method
  ref_dorms = "dorm_ref", # variable with reference def/method
  study_id = "study", # name of the column indicating group membership
  add_intercept = TRUE,
  dorm_separator = " + "
)

fit1 <- CWModel(
  cwdata = dat4,           # result of CWData() function call
  obs_type = "diff_logit", # must be "diff_logit" or "diff_log"
  cov_models = list( # specify covariate details
    CovModel("intercept")),
  gold_dorm = "control",  # level of 'ref_dorms' that's the gold standard
  max_iter=200L,
  inlier_pct = 0.9
)

#calculate adjusted sd and adjusted confidence interval 
#the adjusted sd is clculated by using the sqrt(old sd^2+ gamma)
#after adjusting the values, only two treatments are significant
summaries <- data.table(fit1$create_result_df()[,1:5])
summaries[,`:=`(lower=beta-(qnorm(0.975)*beta_sd), upper=beta+(qnorm(0.975)*beta_sd), p=(1-pnorm(abs(beta/beta_sd)))*2)]
summaries$adjusted_sd <- sqrt(summaries$beta_sd^2 + 0.02087869)
summaries[,`:=`(adjusted_lower=beta-(qnorm(0.975)*adjusted_sd), adjusted_upper=beta+(qnorm(0.975)*adjusted_sd), adjusted_p=(1-pnorm(abs(beta/adjusted_sd)))*2)]


treatment_dir<-paste0("FILEPATH")
table_dir<-paste0("FILEPATH")

if (comparison_tables){
  
  # write descriptive files
  icc<-dt_final_1[,list(n_studies=.N),by=c("dorm_alt","dorm_ref")]
  setnames(icc,names(icc),c("Intervention class","Comparison class","Number of trials"))
  setorderv(icc,c("Number of trials"),order = -1)
  setorderv(icc,"Intervention class")
  write.csv(icc,paste0("FILEPATH"),row.names=F)
  
  isc<-dt_final_1[,list(n_studies=.N),by=c("Instrument/Scale")]
  setnames(isc,names(isc),c("Instrument","Number of trials"))
  write.csv(isc,paste0("FILEPATH"),row.names=F)
}

one_vs_one<-as.data.frame(fit1$fixed_vars)
write.csv(summaries,"FILEPATH", row.names = FALSE)

# prep for plotting 
data <- as.data.table(dt_final_1)
data[, group := paste0(dorm_alt, " - ", dorm_ref)]
data[, `:=` (lower = logit_diff - 1.96*logit_diff_se, upper = logit_diff + 1.96*logit_diff_se)]
data[, trimmed := fit1$w]
data<-data[!duplicated(data$`Study Name`),]
data<-data[data$group %like% "control", ]

#rename columns
colnames(summaries)[colnames(summaries) == 'upper'] <- 'beta_upper'
colnames(summaries)[colnames(summaries) == 'lower'] <- 'beta_lower'
colnames(summaries)[colnames(summaries) == 'dorms'] <- 'dorm_alt'
predicts<-merge(data, summaries[,c('dorm_alt','beta','beta_upper','beta_lower')])

table(predicts$dorm_alt)
predicts <- mutate(predicts, plots_name = ifelse(dorm_alt %in% 'bed_rest', "Best Rest",
                                              ifelse(dorm_alt %in% "pharma_non_op", "Non-opioid analgesics",
                                                     ifelse(dorm_alt %in% 'physio_pharmaNonOp', "Psychological therapies and Non-opioid analgesics",
                                                            ifelse(dorm_alt %in% 'psychological', "Psychological",
                                                                   ifelse(dorm_alt %in% 'education', "Education",
                                                                          ifelse(dorm_alt %in% 'pharma_op', "Opioid analgesics",
                                                                                 ifelse(dorm_alt %in% 'physio_psych', "Psychological and Physical therapies",
                                                                                        ifelse(dorm_alt %in% 'surgery', "Surgery", 
                                                                                               ifelse(dorm_alt %in% 'epidural injection steroids', "Epidural injection steroids",
                                                                                                      ifelse(dorm_alt %in% 'physio', "Physical therapies",
                                                                                                             ifelse(dorm_alt %in% 'physio_psych_pharmaNonOp', "Psychological, Non-opioid analgesics and Physical therapies", dorm_alt))))))))))))



plot_dir <- paste0("FILEPATH")
pdf(paste0("FILEPATH"), width = 12)

for (j in unique(predicts$group)){
  graph_dt<-predicts[predicts$group %in% j, ]
gg<- ggplot() +
    geom_point(data = graph_dt, aes(x = logit_diff, y = `Study Name`, color = as.factor(trimmed)))+
    geom_errorbarh(data = graph_dt, aes(y = `Study Name`, xmin = lower, xmax = upper, color = as.factor(trimmed)),height=0)+
    geom_vline(xintercept = graph_dt$beta, linetype = "dashed", color = "darkorchid")+
    facet_grid(plots_name~.,space="free",scales="free",labeller = labeller(facet_category = label_wrap_gen(width=25)))+
    geom_rect(data = graph_dt[1,], xmin = graph_dt$beta_lower[1], xmax = graph_dt$beta_upper[1], ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "darkorchid")+
    geom_vline(xintercept = 0)+
    labs(x = "Effect Size", y = "") +
    xlim(-3, 2) +
    scale_color_manual(name = "", values = c("1" = "midnightblue", "0" = "red"), labels = c("0" = "Trimmed", "1" = "Included")) +
    ggtitle(paste0(unique(graph_dt$plots_name)," vs. Placebo or Usual care")) +
    theme_bw()+
    theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          axis.text.y = element_text(size = 6,hjust = 0))
print(gg)
}
dev.off()
