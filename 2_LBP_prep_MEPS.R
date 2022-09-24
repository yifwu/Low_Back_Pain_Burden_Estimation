## Docstring ----------------------------------------------------------- ----
## Project: GBD LBP Severity Distribution varying by HAQI Paper
## Details: Step 2: Process MEPS file for utilization estimation of treatment in LBP, Prepare for Mr Bert Run
## Script: 
## --------------------------------------------------------------------- ----

## Environment Prep ---------------------------------------------------- ----

rm(list=ls())
pacman::p_load(foreign, data.table, readxl, stringi, ggplot2, haven, dplyr)
setDTthreads(5)

#read in the raw extraction sheet
meta <- read_excel("FILEPATH of the INPUT DT")

data_dir<-paste0("FILEPATH")
input_dir<-paste0("FILEPATH")
map_dir<-paste0("FILEPATH")
output_dir<-paste0("FILEPATH")
retab<-F

# only panel 4-17 were included, the panels after 18 removed those key variables, physical therapy,
# psycho therapy, and rehab extra. We attemped to include more data but the related variables were all removed
# from the new panel
# MEPS data can be downloaded by the MEPS website
for (MEPS_panel in c(4:17)){   
  
    message(paste("looping through panel "), MEPS_panel)
  if (MEPS_panel %in% c(4:17)){
    # Pull in all medical provider component files
    mc<-as.data.table(read.dta((paste0(data_dir,"/PANEL_",MEPS_panel,"/IND_MEDICAL_CONDITIONS.dta"))))
    dr<-as.data.table(read.dta((paste0(data_dir,"/PANEL_",MEPS_panel,"/PRESCRIBED_MEDICATIONS.dta"))))
    e1<-as.data.table(read.dta((paste0(data_dir,"/PANEL_",MEPS_panel,"/OFFICE_VISITS.dta"))))
    e2<-as.data.table(read.dta((paste0(data_dir,"/PANEL_",MEPS_panel,"/OUTPT_VISITS.dta"))))
    e3<-as.data.table(read.dta((paste0(data_dir,"/PANEL_",MEPS_panel,"/ER_VISITS.dta"))))
    e4<-as.data.table(read.dta((paste0(data_dir,"/PANEL_",MEPS_panel,"/HOSP_INPT_STAYS.dta"))))
    e5<-as.data.table(read.dta((paste0(data_dir,"/PANEL_",MEPS_panel,"/HOME_HEALTH.dta"))))
    link<-as.data.table(read.dta((paste0(data_dir,"/PANEL_",MEPS_panel,"/LINK_FILES.dta"))))
   
    # Panels 4-17 also have an adendum to the prescribed medications file that needs to be pulled in 
    # and merged on
    
      # Read in addendum file 
      ml<-as.data.table(read.dta((paste0(data_dir,"/PANEL_",MEPS_panel,"/MULTINUM_LEXICON.dta"))))
      
      # If both sets of data have therapeutic classes already mapped, use the classification in the 
      # addendum file, not the classification in the original panel data
      merge_vars<-c("dupersid","rxrecidx","datayear")
      drop_vars<-intersect(names(ml),names(dr))[!(intersect(names(ml),names(dr))%in%merge_vars)]
      if(length(drop_vars)>0) dr[,(drop_vars):=NULL]
      
      # Merge on extra info from addendum
      dr<-merge(dr,ml,by=merge_vars,all.x=T)
      
      thera_classes<-grep("tc",names(dr),value=T)
      dr<-dr[,c("dupersid","rxrecidx","linkidx","rxname","rxdrgnam",
                "purchrd",thera_classes,"datayear"),with=F]
      e1<-e1[,.(dupersid,evntidx,eventrn,obdateyr,obdatemm,seetlkpv,
                seedoc,medptype,docatloc,vstctgry,vstrelcn,physth,
                occupth,psychoth,sonogram,anesth,surgproc,medpresc,
                obicd1x,obicd2x,obicd3x,obicd4x,obpro1x,datayear)] 
      
      #in the outpatient, the op is outpatient, procedure is oppro
      e2<-e2[,.(dupersid,evntidx,eventrn,opdateyr,opdatemm,seetlkpv,
                seedoc,medptype,vstctgry,vstrelcn,physth,occupth,
                psychoth,sonogram,anesth,surgproc,medpresc,opicd1x,
                opicd2x,opicd3x,opicd4x,oppro1x)] 
      #emergency room erpro
      e3<-e3[,.(dupersid,evntidx,eventrn,erhevidx,erdateyr,erdatemm,
                seedoc,vstctgry,vstrelcn,sonogram,anesth,surgproc,
                medpresc,ericd1x,ericd2x,ericd3x,erpro1x,datayear)] 
      e4<-e4[,.(dupersid,evntidx,eventrn,ipbegyr,ipbegmm,ipendyr,
                ipendmm,ipenddd,rsninhos,anyoper,ipicd1x,ipicd2x,
                ipicd3x,ipicd4x,ippro1x,ippro2x,dschpmed,datayear)] 
      e5<-e5[,.(dupersid,evntidx,eventrn,hhdateyr,hhdatemm,hhtype,
                cna,compann,dieticn,hhaide,hospice,hmemaker,ivthp,
                medldoc,nurpract,nuraide,occupthp,personal,physlthp,
                respthp,socialw,speecthp,othrhcw,nonskill,skillwos,
                othcw,treatmt,medequip,dailyact,company,othsvce,
                othsvcos,datayear)]
  }
  
    file_list<-list(mc,dr,e1,e2,e3,e4,e5,link)
    
    #invisible is a function that prevents from printing. 
    invisible(lapply(file_list,function(x) setnames(x,names(x),tolower(names(x)))))
    mc[,panel:=as.numeric(gsub("^.*\\D","",panel))]
    icd_map <- as.data.table(
      
    read_excel("FILEPATH"))[,.(icd = `ICD-9`,
                                                              cond_name = paste0("t", `non-fatal cause`),
                                                                   Acute)]
    icd_map <- icd_map[!grepl("N", icd),]
    icd_map <- icd_map[!grepl("V", icd),]
    setnames(icd_map,"icd","icd9codx")
 
    # Read in procedure map
    #INTERVIEW/CONSULT/EXAM 89 (ICD-9)
    #PT, REHAB & RELATED PROC 93 (rehab) (ICD-9)
    #PSYCHE RELATED PROCEDURE (psychological) 94
    proc_map<-as.data.table(read_excel("FILEPATH"))
    proc_map[,proc_label:=tolower(proc_label)]   
    setnames(proc_map, "proc",'proc_description')
    
    # Read in drug code map
    #epidural not listed, manually went back to the map and added, drug class 98, 109 are in 
    #epidural steroids injection 
    drug_map<-as.data.table(read_excel("FILEPATH"))
    drug_map[,drug_class:=tolower(drug_class)]    

    # Merge ICD map on medical conditions file
    #only contains two code, with LBP with leg pain and without legpain
    mc<-merge(mc,icd_map,by="icd9codx")
    
    
    # Only keep relevant variables for mc 
    mc<-mc[,.(condidx,dupersid,panel,condrn,icd9codx,cond_name)]
    
    # Subset medical conditions file to just LBP cases
    mc<-mc[cond_name%in%c("tmsk_pain_lowback_wleg","tmsk_pain_lowback_noleg")]
    
    # Grab condition ids corresponding to LBP
    lbp_condids<-unique(mc[,condidx])
    
    # Grab event ids correspond to LBP condition ids
    lbp_event_ids<-unique(link[condidx%in%lbp_condids,evntidx])

    # Grab rxids that correspond to LBP events
    #getting the medciation that were related to LBP events
    lbp_rx_ids<-unique(link[evntidx%in%lbp_event_ids,linkidx])    
    
    # Grab dupersids that correspond ppl with LBP
    #DUPERSID looks like a unique identifier for a pat across different panels
    lbp_dupersids<-unique(mc[,dupersid])
    
    # Process drug file
    # Subset to rxids that correspond to LBP events
    # Set to dupersids that correpond to respondants with LBP
    dr<-dr[dupersid%in%lbp_dupersids]
    # Melt therapeutic classes
    dr_long<-melt(dr,id.vars = names(dr)[!(names(dr)%in%thera_classes)])
    # Drop rows where value is NA
    dr_long<-dr_long[!is.na(value)]
    # Pull out drug code
    dr_long[,drug_code:=as.numeric(stri_split_fixed(value," ",2,simplify = T)[,1])]
    # Drop rows where drug code is -1 (inapplicable) or -9 (not ascertained)
    dr_long<-dr_long[!(drug_code%in%c(-1,-9))]
    # Merge on mapping file 
    dr_long<-merge(dr_long,drug_map,by=c("drug_code"),all.x=T)
    # Only keep most detailed classifications
    dr_long<-dr_long[most_detailed==1]
    
    #0 rows got epidural injection in 
    dr_long[dr_long$drug_code==109]

    # Fix the miscellaneous analgesic category
    dr_long[grep("TRAM",rxdrgnam),analgesic:=0]
    dr_long[grep("TRAM",rxdrgnam),opioids:=1]
    
    #get the row for epidural
    dr_long[drug_code %in% (109), epidural_steroids_inject:=1]
    
    # Only keep drugs of interest
    dr_long<-dr_long[opioids==1|nsaid==1|analgesic==1|epidural_steroids_inject==1]
    # Make unique by individual
    dr_long<-dr_long[,.(opioids=max(opioids),nsaid=max(nsaid),analgesic=max(analgesic),epidural_steroids_inject=max(epidural_steroids_inject)),
                     by=c("dupersid")]
    dr_long$epidural_steroids_inject <- replace(dr_long$epidural_steroids_inject, is.na(dr_long$epidural_steroids_inject), 0)
  
    
    # Process each event file
    # Office visits
    # Subset to visits for lbp events
    e1<-e1[evntidx%in%lbp_event_ids]
    # Only count surgical procedures where patient also got anesthesia
    message(paste("Office visits: Surgical procedures + anesthesia combos for panel "),MEPS_panel)
    print(e1[,.N,by=c("surgproc","anesth")])
    
    #It looks like there is surgery only but does not have a lot about the typpe of surgery
    e1[surgproc=="1 YES" & anesth=="1 YES", surg:=1]
    e1[is.na(surg),surg:=0]
    
    # Count anyone receiving physth, occupth, or psychoth as msk intervention
    #physth: physcial therapy
    #occupth: OCCUPATIONAL THERAPY
    #psychoth: DID P HAVE PSYCHOTHERAPY/COUNSELING
    e1[,.N,by=c("physth","occupth","psychoth")]
    e1[physth=="1 YES" | occupth == "1 YES" | obpro1x %in% (93), physio:=1]
    #adding psych category and physio category
    e1[psychoth == "1 YES" | obpro1x %in% (94), psych:=1]
    e1[is.na(physio),physio:=0]
    
    e1[is.na(psych),psych:=0]
    
    #also adding combined categories in there
    e1[psych ==1 & physio==1, physio_psych:=1]
    e1[is.na(physio_psych),physio_psych:=0]
    
    # Make unique by individual
    e1<-e1[,.(physio=max(physio),surg=max(surg), psych=max(psych), physio_psych=max(physio_psych)),by=c("dupersid")]
    
    # Outpatient visits
    # Subset to visits for lbp events
    e2<-e2[evntidx%in%lbp_event_ids]
    # Only count surgical procedures where patient also got anesthesia
    message(paste("Outpatient visits: Surgical procedures + anesthesia combos for panel "),MEPS_panel)
    print(e2[,.N,by=c("surgproc","anesth")])
    e2[surgproc=="1 YES" & anesth=="1 YES", surg:=1]
    e2[is.na(surg),surg:=0]
    
    # Count anyone receiving physth, occupth, or psychoth as msk intervention
    e2[,.N,by=c("physth","occupth","psychoth")]
    e2[physth=="1 YES" | occupth == "1 YES" | oppro1x %in% (93), physio:=1]
    #adding psych category and physio category
    e2[psychoth == "1 YES" | oppro1x %in% (94), psych:=1]
    e2[is.na(physio),physio:=0]
    
    e2[is.na(psych),psych:=0]
    
    #also adding combined categories in there
    e2[psych ==1 & physio==1, physio_psych:=1]
    e2[is.na(physio_psych),physio_psych:=0]
    
    # Make unique by individual
    e2<-e2[,.(physio=max(physio),surg=max(surg), psych=max(psych), physio_psych=max(physio_psych)),by=c("dupersid")]
    
    # ER visits
    # Subset to visits for lbp events
    e3<-e3[evntidx%in%lbp_event_ids]
    # Only count surgical procedures where patient also got anesthesia
    message(paste("ER visits: Surgical procedures + anesthesia combos for panel "),MEPS_panel)
    print(e3[,.N,by=c("surgproc","anesth")])
    e3[surgproc=="1 YES" & anesth=="1 YES", surg:=1]
    e3[is.na(surg),surg:=0]
    # Assuming you can't get msk interventions at the ER
    e3[,physio:=0]
    e3[,psych:=0]
    e3[,physio_psych:=0]
    # Make unique by individual
    e3<-e3[,.(physio=max(physio),surg=max(surg), psych=max(psych), physio_psych=max(physio_psych)),by=c("dupersid")]
    
    # Hospital inpatient stays
    # Subset to visits for lbp events
    e4<-e4[evntidx%in%lbp_event_ids]
    # Count any operation as surgery
    e4[anyoper=="1 YES", surg:=1]
    e4[is.na(surg),surg:=0]
    # Assuming you can't get msk interventions as an inpt
    e4[,physio:=0]
    e4[,psych:=0]
    e4[,physio_psych:=0]
    # Make unique by individual
    e4<-e4[,.(physio=max(physio),surg=max(surg), psych=max(psych), physio_psych=max(physio_psych)),by=c("dupersid")]
    
    # Home health
    # Subset to visits for lbp events
    e5<-e5[evntidx%in%lbp_event_ids]
    # Count anyone getting visited by a physical or occupational therapist as
    # getting an msk intervention 
    message(paste("Home health: PT and OT combos for panel "),MEPS_panel)
    print(e5[,.N,by=c("physlthp","occupthp")])
    e5[physlthp =="1 YES" | occupthp == "1 YES", physio:=1]
    e5[is.na(physio),physio:=0]
    # Assuming you can't get surg interventions at home
    e5[,surg:=0]
    #could not find any psychological related variables, treating it as 0
    e5[,psych:=0]
    e5[,physio_psych:=0]
    # Make unique by individual
    e5<-e5[,.(physio=max(physio),surg=max(surg), psych=max(psych), physio_psych=max(physio_psych)),by=c("dupersid")]

    # Combine all event files
    evnts<-rbindlist(list(e1,e2,e3,e4,e5),use.names=T) 

    # Merge on drug information
    full_data<-merge(evnts,dr_long,by=c("dupersid"),all=T)    
    
    #creating a few new categories
    #combine nsaid and analgesic together
    full_data[(nsaid==1 | analgesic==1) & opioids!=1, pharmaNonOp:=1 ]
    full_data[physio==1 & pharmaNonOp==1, physio_pharmaNonOp:=1]
    full_data[physio==1 & psych==1 & pharmaNonOp==1, physio_psych_pharmaNonOp:=1]
    
    # Fill in NAs with 0s
    trt_cols<-c("surg","opioids","physio","psych","physio_psych","nsaid","analgesic","epidural_steroids_inject","pharmaNonOp","physio_pharmaNonOp","physio_psych_pharmaNonOp")
    full_data[,(trt_cols):=lapply(.SD,function(x) ifelse(is.na(x),0,x)),
              .SDcols=trt_cols]

    # Make unique by individual
    full_data<-full_data[,.(physio=max(physio),
                            physio_psych= max(physio_psych),
                            surg=max(surg),
                            opioids=max(opioids),
                            nsaid=max(nsaid),
                            analgesic=max(analgesic),
                            psych=max(psych),
                            epidural_steroids_inject=max(epidural_steroids_inject),
                            pharmaNonOp=max(pharmaNonOp),
                            physio_pharmaNonOp=max(physio_pharmaNonOp),
                            physio_psych_pharmaNonOp=max(physio_psych_pharmaNonOp)),
                         by=c("dupersid")]
    
    # Merge utilization data onto lbp sample
    mc_full<-merge(mc,full_data,by=c("dupersid"),all=T)  
    mc_full[, id:= paste0(dupersid, MEPS_panel)]

    # Fill in NAs with 0s
    mc_full[,(trt_cols):=lapply(.SD,function(x) ifelse(is.na(x),0,x)),
            .SDcols=trt_cols]
    
    # Make unique by individual and condition
    mc_full<-mc_full[,.(physio=max(physio),
                        physio_psych= max(physio_psych),
                        surg=max(surg),
                        opioids=max(opioids),
                        nsaid=max(nsaid),
                        analgesic=max(analgesic),
                        psych=max(psych),
                        epidural_steroids_inject=max(epidural_steroids_inject),
                        pharmaNonOp=max(pharmaNonOp),
                        physio_pharmaNonOp=max(physio_pharmaNonOp),
                        physio_psych_pharmaNonOp=max(physio_psych_pharmaNonOp)),
                     by=c("dupersid","cond_name","icd9codx","panel")]
    mc_full[, id:= paste0(dupersid, MEPS_panel)]
    
    assign(paste0("mc_full_",MEPS_panel),copy(mc_full))
}
full<-rbindlist(lapply(paste0("mc_full_",4:17),get),use.names=T)

#read in coefficient dataset
coefficient <- read.csv("FILEPATH")

# Bind all utilization datasets together 
# the other lines are not created yet 
output <- full[,-c('nsaid','analgesic')]
output[, i := .I]
#calculate the sum of treatment they received? 
output[,sum:= sum(physio, surg, opioids, psych, epidural_steroids_inject, pharmaNonOp), by=i]

#change wide to long dataset
output_long<-melt(output,
                    id.vars = c("dupersid","cond_name","icd9codx","panel","id",'sum','i'),
                    variable.name = "dorms",
                    value.name = "n")
beta_output<- merge(output_long, coefficient[, c('dorms','beta')], by='dorms')
sort<-beta_output[order(id, beta, cond_name)]
sort[, pro_beta_n:= 1-abs(beta)*n]

#here it is doing the 1-(1-beta*n)(1-beta*n)
out<-sort %>%
group_by(id, cond_name) %>%
  summarize(Product = prod(pro_beta_n, na.rm = TRUE)) %>%
  inner_join(sort, ., by = c('id', 'cond_name'))

data_wide <- as.data.table(dcast(out, ...~dorms, value.var="n", fun.aggregate = sum))
data_wide[, adjust_es:= -(1-Product)]

write.csv(data_wide,paste0(output_dir,"adjust_beta_all_panels_long.csv"),row.names=F )

full<-as.data.table(read.csv(paste0(output_dir,"utilization_all_panels_long.csv")))

# Tabulate utilization - first check by panel and condition 
full_tab<-full[,.(physio=mean(physio),
                  physio_psych= mean(physio_psych),
                  surg=mean(surg),
                  opioids=mean(opioids),
                  nsaid=mean(nsaid),
                  analgesic=mean(analgesic),
                  psych=mean(psych),
                  epidural_steroids_inject=mean(epidural_steroids_inject),
                  pharmaNonOp=mean(pharmaNonOp),
                  physio_pharmaNonOp=mean(physio_pharmaNonOp),
                  physio_psych_pharmaNonOp=mean(physio_psych_pharmaNonOp),
                  n=.N),
               by=c("panel","icd9codx","cond_name")]
setorderv(full_tab,c("cond_name","panel"))

# Plot 
pdat<-melt(full_tab,id.vars=c("panel","icd9codx","cond_name","n"))
pdat[variable=="physio",intervention:="Physical therapies"]
pdat[variable=="surg",intervention:="Surgical interventions"]
pdat[variable=="opioids",intervention:="Opioid analgesics"]
pdat[variable=="pharmaNonOp",intervention:="Non-opioid, NSAID analgesics"]
pdat[variable=="physio_psych",intervention:="CBT + Physical therapies"]
pdat[variable=="psych",intervention:="Cognitive behavioral therapies"]
pdat[variable=="epidural_steroids_inject",intervention:="Epidural steroids injection"]
pdat[variable=="physio_psych_pharmaNonOp",intervention:="CBT + Physical therapies + N
     on-opioid, NSAID analgesics"]
pdat[variable=="physio_pharmaNonOp",intervention:="Physical therapies + Non-opioid, NSAID analgesics"]
pdat <- pdat[complete.cases(pdat$intervention),]

pdat[cond_name=="tmsk_pain_lowback_noleg",condition:="LBP - no leg involvement"]
pdat[cond_name=="tmsk_pain_lowback_wleg",condition:="LBP - leg involvement"]
pdat[,value_sd:=sqrt((value*(1-value))/n)]
pdat[,lower:=(value-1.96*value_sd)*100]
pdat[,upper:=(value+1.96*value_sd)*100]
pdat[,value:=value*100]

pdat_individual <- pdat[pdat$intervention %in% c("Cognitive behavioral therapies","Epidural steroids injection",
                                                 "Non-opioid, NSAID analgesics","Opioid analgesics",
                                                 "Physical therapies","Surgical interventions"),]

pdat_combined <- pdat[pdat$intervention %in% c("CBT + Physical therapies + N\n     on-opioid, NSAID analgesics",
                                               "Physical therapies + Non-opioid, NSAID analgesics",
                                               "CBT + Physical therapies"),]

#individual categories plots
util_gg<-ggplot(data=pdat_individual,aes(x=panel,y=value,color=intervention))+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=intervention),linetype=0, alpha = 0.1)+
  geom_line(size=1)+
  geom_point(size=2)+
  theme_bw()+
  facet_grid(~condition)+
  ylab("Utilization (%)")+
  xlab("MEPS Panel")+
  ggtitle("Intervention-specific utilization among MEPS respondents-Individual treatment")+
  scale_y_continuous(limits=c(0,80),expand=c(0,0))+
  theme(strip.background = element_blank(),strip.text.x = element_text(size=14,hjust = 0),
        legend.text = element_text(size=10), legend.title = element_text(size = 12),
        legend.position = "bottom")
util_gg
ggsave(util_gg,filename = paste0(output_dir,"utilization_plot_Individual treatment.png"),width=15,height=7)

#combined treatment plots
util_gg<-ggplot(data=pdat_combined,aes(x=panel,y=value,color=intervention))+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=intervention),linetype=0, alpha = 0.1)+
  geom_line(size=1)+
  geom_point(size=2)+
  theme_bw()+
  facet_grid(~condition)+
  ylab("Utilization (%)")+
  xlab("MEPS Panel")+
  ggtitle("Intervention-specific utilization among MEPS respondents-combined treatment")+
  scale_y_continuous(limits=c(0,80),expand=c(0,0))+
  theme(strip.background = element_blank(),strip.text.x = element_text(size=14,hjust = 0),
        legend.text = element_text(size=10), legend.title = element_text(size = 12),
        legend.position = "bottom")
util_gg
ggsave(util_gg,filename = paste0(output_dir,"utilization_plot_Combined treatment.png"),width=15,height=7)

# Pool across panel
full_tab<-full[,.(physio=mean(physio),
                   physio_psych= mean(physio_psych),
                   surg=mean(surg),
                   opioids=mean(opioids),
                   nsaid=mean(nsaid),
                   analgesic=mean(analgesic),
                   psych=mean(psych),
                   epidural_steroids_inject=mean(epidural_steroids_inject),
                   pharmaNonOp=mean(pharmaNonOp),
                   physio_pharmaNonOp=mean(physio_pharmaNonOp),
                   physio_psych_pharmaNonOp=mean(physio_psych_pharmaNonOp),
                   n=.N),
               by=c("icd9codx","cond_name")]

full_tab<-melt(full_tab,id.vars = c("icd9codx","cond_name","n"))
full_tab[,value_sd:=sqrt((value*(1-value))/n)]
full_tab[,lower:=(value-1.96*value_sd)]
full_tab[,upper:=(value+1.96*value_sd)]
full_tab[,value:=value]

coefficient <- read.csv("FILEPATH")
coefficient <- coefficient[,c("dorms", "beta")]
names(coefficient)[names(coefficient) == "dorms"] <- "variable"

full_t <- merge(full_tab, coefficient, all.x=TRUE)

# Write this file
write.csv(full_t,paste0(output_dir,"utilization_by_condition_with_sd.csv"),row.names=F)
