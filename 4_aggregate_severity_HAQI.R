## Docstring ----------------------------------------------------------- ----
## Project: GBD LBP Severity Distribution varying by HAQI Paper
## Script: 
## Description: Calculate severity distribution by HAQI
## --------------------------------------------------------------------- ----
rm(list=ls())

version_id <- 3

pacman::p_load(msm, data.table, readstata13, ggplot2, openxlsx, scales,readxl)

set.seed(1640)

data_dir<-paste0("PATH")
input_dir<-paste0("PATH")

map_dir<-paste0("PATH")
output_dir<-paste0("PATH")
plot_dir <- paste0("PATH")

prev_all <- fread("prevalence_agg_2019.csv")

#aggregate the with leg pain and without leg pain together
prev_all <- prev_all[, sum(prev), by = c("metric_id", 'age_group_id', 'location_id',
                                           'sex_id', 'year_id',"version_id", "draw_num")]
names(prev_all)[names(prev_all) == 'V1'] <- 'prev'

lbp_dws<- fread(paste0(output_dir,"lbp_dws_intermediate_",version_id,".csv"))

lbp_dws[, dw_adj_factor := dw_mean/dw_mean_adj]
lbp_dws[, dw_fuot_factor := dw_mean/dw_mean_fcot]

adj_factor <- unique(lbp_dws[,.(draw_num, dw_adj_factor)])
adj_factor[,.(dw_adj_factor = mean(dw_adj_factor), lower = quantile(dw_adj_factor, 0.025), upper = quantile(dw_adj_factor, 0.975))]

raw_dws <- unique(lbp_dws[,.(raw_dw = dw_mean, draw_num)])
fuot_dws <- unique(lbp_dws[,.(fuot_dw = dw_mean_fcot, draw_num)])
adj_dws <- unique(lbp_dws[,.(adj_dw = dw_mean_adj, draw_num)])

##### Predict adjustment factor by HAQI  #####
adj_factor[, haqi := 0]
haqi<- as.data.table(read_excel(paste0(output_dir, 'haqi_2007.xlsx')))
haqi <- data.table(draw_num = paste0("draw_", 0:999), dw_adj_factor = 1, haqi = rnorm(1000, haqi$mean_value, (haqi$upper_value - haqi$lower_value)/(qnorm(0.975, 0, 1)*2)))
haqi[,.(mean(haqi), quantile(haqi, 0.025), quantile(haqi, 0.975))]
adj_factor <- rbind(adj_factor, haqi)

pred_matrix <- data.table(expand.grid(haqi = seq(0, 100, 0.5), draw_num = paste0("draw_", 0:999)))

interpolation_betas <- data.table(draw_num = paste0("draw_", 0:999))
for(d in paste0("draw_", 0:999)){
  adj_mod<-lm(dw_adj_factor~haqi,subset = (draw_num == d), data = adj_factor)
  new <- data.frame(haqi = seq(0, 100, 0.5))
  pred_matrix[draw_num == d, adj_factor := predict(adj_mod, newdata = new)]
  interpolation_betas[draw_num == d, haqi_beta := adj_mod$coefficients["haqi"]]
}
interpolation_betas[,.(mean(haqi_beta), quantile(haqi_beta, 0.025), quantile(haqi_beta, 0.975))]

line_plot <- copy(pred_matrix)
line_plot[, `:=` (mean = mean(adj_factor), lower = quantile(adj_factor, 0.025), upper = quantile(adj_factor, 0.975)), by = "haqi"]
line_plot <- unique(line_plot[,.(haqi, adj_factor = mean, lower, upper)])

plot <- ggplot(data=line_plot, aes(x=haqi, y=adj_factor))+
  geom_line(size=1, color="blue") +
  geom_ribbon(data= line_plot, aes(x=haqi, ymin=lower, ymax=upper),  fill="blue", alpha=.3) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 100, 10), limits = c(0, 101)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0.5, 1.1, 0.1),  limits = c(0.5, 1.1))+
  ylab("Adjustment factor to sequela-weighted disability weight") +
  xlab("Healthcare access quality index") +
  theme(axis.line=element_line(colour="black"),
        legend.title=element_blank())
plot
ggsave(filename = paste0(plot_dir,paste0("adjust_factor.pdf")),plot = plot, width = 8,height = 6)


# pull in LBP prevalence estimates
locs <-as.data.table(read_excel(paste0(output_dir, 'locs.xlsx')))
loc_ids<-locs[level == 3,location_id]

#HAQI by country
haqi_by_country<-as.data.table(read_excel(paste0(output_dir, 'haqi_by_country.xlsx')))

haqi_by_country[, (paste0("draw_", 0:999)) := as.list(rnorm(n = 1000, mean = mean_value, sd = (upper_value - lower_value)/(qnorm(0.975, 0, 1)*2))), by = "location_id"]
haqi_by_country <- melt.data.table(haqi_by_country, id.vars = names(haqi_by_country)[!(names(haqi_by_country) %like% "draw")], value.name="haqi", variable.name="draw_num")

prevalence <- merge(prev_all[,.(location_id, draw_num, prev)], haqi_by_country[,.(location_id, draw_num, haqi)], all.x = T, by = c("location_id", "draw_num"))

prevalence[, haqi := (5*round((haqi*10)/5))/10]

prevalence <- merge(prevalence, pred_matrix, all.x = T, by = c("haqi", "draw_num"))

prevalence <- merge(prevalence, raw_dws, all.x = T, by = c("draw_num"),allow.cartesian=TRUE) # merge with raw DWs

prevalence[, adj_dw := raw_dw / adj_factor] # estimate haqi-specific DWs

prevalence <- merge(prevalence, fuot_dws, all.x = T, by = c("draw_num"),allow.cartesian=TRUE) # merge with fcot DWs

prevalence <- merge(prevalence, adj_dws[,.(draw_num, worst_dw = adj_dw)], all.x = T,  by = "draw_num",allow.cartesian=TRUE) # merge with no-treatment DWs

#get population
population<-get_population(location_id = loc_ids,
                          year_id = 2019,
                          gbd_round_id = 7,
                          sex_id = 3,
                          age_group_id = 22,
                         decomp_step="iterative")
population[,run_id:=NULL]

prevalence <- merge(prevalence, population[,.(location_id, population)], all.x = T, by = "location_id")

prevalence <- merge(prevalence, locs[,.(location_id, super_region_name, region_name)], all.x = T, by = "location_id")

best_haqi <- unique(haqi_by_country[mean_value == max(mean_value),location_id])
haqi_by_country[mean_value == max(mean_value)]
haqi_by_country[mean_value == min(mean_value)]
best_dw <- prevalence[location_id == best_haqi, .(draw_num, best_dw = adj_dw)]

prevalence <- merge(prevalence, best_dw, all.x = T, by = "draw_num",allow.cartesian=TRUE)

prevalence[, cases := population * prev]

prevalence[, `:=` (raw_yld = raw_dw * cases, adj_yld = adj_dw * cases, worst_yld = worst_dw * cases, best_yld = best_dw * cases, fuot_ylds = fuot_dw * cases)]

# Aggregate to region ------------------------------------------------------

prevalence[, `:=` (cases_region = sum(cases), population_region = sum(population), raw_yld_region = sum(raw_yld), worst_yld_region = sum(worst_yld), adj_yld_region = sum(adj_yld), best_yld_region = sum(best_yld), fuot_yld_region = sum(fuot_ylds)), by = c("draw_num", "region_name")]
prevalence[, `:=` (cases_super_region = sum(cases), population_super_region = sum(population), raw_yld_super_region = sum(raw_yld), worst_yld_super_region = sum(worst_yld), adj_yld_super_region = sum(adj_yld), best_yld_super_region = sum(best_yld), fuot_yld_super_region = sum(fuot_ylds)), by = c("draw_num", "super_region_name")]
prevalence[, `:=` (cases_global = sum(cases), population_global = sum(population), raw_yld_global = sum(raw_yld), worst_yld_global = sum(worst_yld), adj_yld_global = sum(adj_yld), best_yld_global = sum(best_yld), fuot_yld_global = sum(fuot_ylds)), by = c("draw_num")]

super_region_data <- unique(prevalence[,.(super_region_name, population = population_super_region, raw_yld = raw_yld_super_region, worst_yld = worst_yld_super_region, adj_yld = adj_yld_super_region, best_yld = best_yld_super_region, fuot_yld = fuot_yld_super_region, draw_num)])
super_region_data[, `:=` (averted_yld = (worst_yld - adj_yld)/worst_yld, brc_yld = (adj_yld - best_yld)/worst_yld, fuot_yld = (best_yld - fuot_yld)/worst_yld)]
super_region_data[, `:=` (unavoidable_yld = 1 - averted_yld - brc_yld - fuot_yld)]

print(super_region_data[,.(averted_yld = mean(averted_yld), averted_lower = quantile(averted_yld, 0.025), averted_upper = quantile(averted_yld, 0.975)), by = "super_region_name"])
print(super_region_data[,.(brc_yld = mean(brc_yld), brc_lower = quantile(brc_yld, 0.025), brc_upper = quantile(brc_yld, 0.975)), by = "super_region_name"])
print(super_region_data[,.(fuot_yld = mean(fuot_yld), fuot_lower = quantile(fuot_yld), fcot_upper = quantile(fuot_yld, 0.975)), by = "super_region_name"])
print(super_region_data[,.(remaining_yld = mean(worst_yld), remaining_lower = quantile(worst_yld, 0.025), remaining_upper = quantile(worst_yld, 0.975)), by = "super_region_name"])

print(super_region_data[,.(total_BRC = mean(brc_yld+averted_yld), lower = quantile(brc_yld+averted_yld, 0.025), upper = quantile(brc_yld+averted_yld, 0.975)), by = "super_region_name"])
print(super_region_data[,.(total_FUOT = mean(fuot_yld+brc_yld+averted_yld), lower = quantile(fuot_yld+brc_yld+averted_yld, 0.025), upper = quantile(fuot_yld+brc_yld+averted_yld, 0.975)), by = "super_region_name"])

print(super_region_data[,.(total_avoidable = mean(brc_yld+fuot_yld), lower = quantile(brc_yld+fuot_yld, 0.025), upper = quantile(brc_yld+fuot_yld, 0.975)), by = "super_region_name"])

super_region_data_summary <- super_region_data[,.(averted_yld = mean(averted_yld), brc_yld = mean(brc_yld), fuot_yld = mean(fuot_yld), unavoidable_yld = mean(unavoidable_yld)), by = "super_region_name"]
super_region_data_summary <- melt.data.table(super_region_data_summary, id.vars = names(super_region_data_summary)[!(names(super_region_data_summary) %like% "yld")], value.name="value", variable.name="YLDs")

super_region_data_summary[, value := value * 100]

super_region_data_summary[, super_region_name := gsub(", ", ",\n", super_region_name)]
super_region_data_summary[, super_region_name := gsub(" and", "\nand", super_region_name)]
super_region_data_summary[, super_region_name := gsub("Saharan", "Saharan\n", super_region_name)]
super_region_data_summary[YLDs == "averted_yld", YLDs := "Averted"]
super_region_data_summary[YLDs == "brc_yld", YLDs := "Avoidable: BRC"]
super_region_data_summary[YLDs == "fuot_yld", YLDs := "Avoidable: FUOT"]
super_region_data_summary[YLDs == "unavoidable_yld", YLDs := "Remaining"]

mean(global_data$averted_yld)
mean(global_data$brc_yld)
#new changes, appending global data into the figure
df2 <- data.frame(super_region_name=c('Global','Global','Global','Global'),
                  YLDs=c("Averted","Avoidable: BRC","Avoidable: FUOT","Remaining"),
                  value=c(15.81639,6.953258,10.213232,67.017118))
df2
super_region_data_summary<-rbind(super_region_data_summary, df2)

gg<-ggplot(data=super_region_data_summary)+
  geom_bar(aes(y=value,x=super_region_name,fill=YLDs),stat="identity")+
  theme_bw()+
  scale_y_continuous(expand=c(0,0))+
  scale_x_discrete(labels=wrap_format(16))+
  ylab("Percent of total burden")+
  xlab("GBD super region")+
  theme(axis.title.x=element_blank(),
        legend.title=element_blank())+
  scale_fill_manual(values=c("#ff6600","#9BC53D","#00A5CE","#006594"))
gg

ggsave(gg, filename=paste0(plot_dir, "ylds_by_super_region_",4,".pdf"), width = 10, height = 6)


global_data <- unique(prevalence[,.(population = population_global, raw_yld = raw_yld_global, worst_yld = worst_yld_global, adj_yld = adj_yld_global, best_yld = best_yld_global, fuot_yld = fuot_yld_global, draw_num)])
global_data[, `:=` (averted_yld = (worst_yld - adj_yld)/worst_yld, brc_yld = (adj_yld - best_yld)/worst_yld, fuot_yld = (best_yld - fuot_yld)/worst_yld)]
global_data[, `:=` (unavoidable_yld = 1 - averted_yld - brc_yld - fuot_yld)]

global_data[,.(averted_yld = mean(averted_yld), averted_lower = quantile(averted_yld, 0.025), averted_upper = quantile(averted_yld, 0.975))]
global_data[,.(brc_yld = mean(brc_yld), brc_lower = quantile(brc_yld, 0.025), brc_upper = quantile(brc_yld, 0.975))]
global_data[,.(fuot_yld = mean(fuot_yld), fuot_lower = quantile(fuot_yld,0.025), fuot_upper = quantile(fuot_yld, 0.975))]
global_data[,.(remaining_yld = mean(worst_yld), remaining_lower = quantile(worst_yld, 0.025), remaining_upper = quantile(worst_yld, 0.975))]

global_data[,.(total_BRC = mean(brc_yld+averted_yld), lower = quantile(brc_yld+averted_yld, 0.025), upper = quantile(brc_yld+averted_yld, 0.975))]
global_data[,.(total_FUOT = mean(fuot_yld+brc_yld+averted_yld), lower = quantile(fuot_yld+brc_yld+averted_yld, 0.025), upper = quantile(fuot_yld+brc_yld+averted_yld, 0.975))]


# Create table with country-level estimates -------------------------------

prevalence[, `:=` (adj_dw_region = adj_yld_region / cases_region)]
prevalence[, `:=` (adj_dw_super_region = adj_yld_super_region / cases_super_region)]
prevalence[, `:=` (adj_dw_global = adj_yld_global / cases_global)]

prevalence[, `:=` (dw_change_percent = 100*(adj_dw - raw_dw)/raw_dw, dw_change_percent_region = 100*(adj_dw_region-raw_dw)/raw_dw, dw_change_percent_super_region = 100*(adj_dw_super_region-raw_dw)/raw_dw, dw_change_percent_global = 100*(adj_dw_global-raw_dw)/raw_dw)]

prevalence[, `:=` (haqi_region = sum(haqi * population)/population_region), by = c("draw_num", "region_name")]
prevalence[, `:=` (haqi_super_region = sum(haqi * population)/population_super_region), by = c("draw_num", "super_region_name")]
prevalence[, `:=` (haqi_global = sum(haqi * population)/population_global), by = "draw_num"]

#taking average
raw_sev_props <- dcast(lbp_dws[,.(draw_num, severity, sev_prop)], ...~severity, value.var="sev_prop", fun.aggregate =mean)
names(raw_sev_props) <- c("draw_num", "raw_asymp_prop", "raw_mild_prop", "raw_mod_prop", "raw_sev_prop", "raw_most_sev_prop")
raw_sev_dws <- dcast(lbp_dws[,.(draw_num, severity, dw)], ...~severity, value.var="dw", fun.aggregate=mean)
names(raw_sev_dws) <- c("draw_num", "asymp_dw", "mild_dw", "mod_dw", "sev_dw", "most_sev_dw")

prevalence <- merge(prevalence, raw_sev_props, by = "draw_num")
prevalence <- merge(prevalence, raw_sev_dws, by = "draw_num")
prevalence[, `:=` (dw_diff = mild_dw * (raw_mild_prop/adj_factor) + mod_dw * (raw_mod_prop/adj_factor) + sev_dw * (raw_sev_prop/adj_factor) + 
                     most_sev_dw *(raw_most_sev_prop/adj_factor) + asymp_dw * (1-(raw_mild_prop+raw_mod_prop+raw_sev_prop+ raw_most_sev_prop)/adj_factor) - raw_dw/adj_factor) ]
prevalence[, `:=` (p_diff = dw_diff / (asymp_dw - most_sev_dw)) ]
prevalence[, `:=` (adj_mild_prop = raw_mild_prop/adj_factor, adj_mod_prop = raw_mod_prop/adj_factor, adj_sev_prop = raw_sev_prop/adj_factor, adj_most_sev_prop= raw_most_sev_prop/adj_factor+ p_diff)]
prevalence[, `:=` (adj_asymp_prop = 1 - adj_mild_prop - adj_mod_prop - adj_sev_prop - adj_most_sev_prop)]

## Correct draws where asymptomatic proportions are negative
prevalence[adj_asymp_prop < 0, dw_diff_2 := mod_dw * adj_mod_prop + sev_dw * adj_sev_prop +most_sev_dw * adj_most_sev_prop  + mild_dw * (1- adj_mod_prop - adj_sev_prop - adj_most_sev_prop) - adj_dw]
prevalence[adj_asymp_prop < 0, p_diff_2 := dw_diff_2 / (mild_dw - most_sev_dw)]
prevalence[adj_asymp_prop < 0, adj_most_sev_prop := adj_most_sev_prop + p_diff_2]
prevalence[adj_asymp_prop < 0, `:=` (adj_asymp_prop = 0, adj_mild_prop = 1-adj_mod_prop-adj_sev_prop -adj_most_sev_prop )]

## Correct draws where mild proportions are negative
prevalence[adj_mild_prop < 0, dw_diff_3 := sev_dw * adj_sev_prop + most_sev_dw * adj_most_sev_prop + mod_dw * (1- adj_sev_prop -adj_most_sev_prop) - adj_dw]
prevalence[adj_mild_prop < 0, p_diff_3 := dw_diff_3 / (mod_dw - most_sev_dw)]
prevalence[adj_mild_prop < 0, adj_most_sev_prop := adj_most_sev_prop + p_diff_3]
prevalence[adj_mild_prop < 0, `:=` (adj_u.mild_prop = 0, adj_mod_prop = 1-adj_sev_prop -adj_most_sev_prop)]

## Correct draws where moderate proportions are negative
prevalence[adj_mild_prop < 0, dw_diff_4 := most_sev_dw * adj_most_sev_prop + sev_dw * (1-adj_most_sev_prop) - adj_dw]
prevalence[adj_mild_prop < 0, p_diff_4 := dw_diff_4 / (sev_dw - most_sev_dw)]
prevalence[adj_mild_prop < 0, adj_most_sev_prop := adj_most_sev_prop + p_diff_4]
prevalence[adj_mild_prop < 0, `:=` (adj_mod_prop = 0, adj_sev_prop = 1-adj_most_sev_prop)]

pres_results_1 <- function(m){return(paste0(sprintf("%.1f", mean(m)), " (", sprintf("%.1f", quantile(m, 0.025)), "–", sprintf("%.1f", quantile(m, 0.975)), ")"))}
pres_results_3 <- function(m){return(paste0(sprintf("%.3f", mean(m)), " (", sprintf("%.3f", quantile(m, 0.025)), "–", sprintf("%.3f", quantile(m, 0.975)), ")"))}

region_table <- unique(prevalence[,.(location_name = region_name, haqi = pres_results_1(haqi_region), raw_dw = pres_results_3(raw_dw), adj_dw = pres_results_3(adj_dw_region), dw_change_percent = pres_results_1(dw_change_percent_region)), by = "region_name"])
super_region_table <- unique(prevalence[,.(location_name = super_region_name, haqi = pres_results_1(haqi_super_region), raw_dw = pres_results_3(raw_dw), adj_dw = pres_results_3(adj_dw_super_region), dw_change_percent = pres_results_1(dw_change_percent_super_region)), by = "super_region_name"])
global_table <- unique(prevalence[,.(location_name = "Global", haqi = pres_results_1(haqi_global), raw_dw = pres_results_3(raw_dw), adj_dw = pres_results_3(adj_dw_global), dw_change_percent = pres_results_1(dw_change_percent_global))])

region_table <- merge(region_table, locs[,.(location_name, lancet_label, sort_order)], by = "location_name")
super_region_table <- merge(super_region_table, locs[,.(location_name, lancet_label, sort_order)], by = "location_name")
global_table <- merge(global_table, locs[,.(location_name, lancet_label, sort_order)], by = "location_name")

country_table <- unique(prevalence[,.(haqi = pres_results_1(haqi), raw_dw = pres_results_3(raw_dw), adj_dw = pres_results_3(adj_dw), dw_change_percent = pres_results_1(dw_change_percent)), by = "location_id"])
country_table <- merge(country_table, locs[,.(location_id, lancet_label, sort_order)], by = "location_id")

final_table <- rbind(country_table, region_table, super_region_table, global_table, fill = T)
final_table <- final_table[order(sort_order), .(Location = lancet_label, HAQI = gsub("\\.", "·", haqi), `Original DW` = gsub("\\.", "·", raw_dw), `Adjusted DW` = gsub("\\.", "·", adj_dw), `DW % change` = gsub("\\.", "·", dw_change_percent))]
final_table[, duplicate := seq_len(.N), by = "Location"]
final_table <- final_table[duplicate == 1,]
final_table[, duplicate := NULL]

write.xlsx(final_table, paste0(output_dir,"appdendix_table_",version_id,".xlsx"))

prevalence[, dw_mean := mean(adj_dw), by = "location_id"]
print(locs[location_id == unique(prevalence[dw_mean == min(dw_mean),location_id]), location_name])

final_table[, num := as.numeric(substring(`Adjusted DW`, 3, 6))]
final_table[num == min(num),]
final_table[num == max(num),]
# Create table with country-level severity splits -------------------------------
prevalence[, `:=` (asymp_region = sum(adj_asymp_prop * population)/population_region,
                   mild_region = sum(adj_mild_prop * population)/population_region,
                   mod_region = sum(adj_mod_prop * population)/population_region,
                   sev_region = sum(adj_sev_prop*population)/population_region,
                   most_sev_region  = sum(adj_most_sev_prop*population)/population_region), by = c("draw_num", "region_name")]

prevalence[, `:=` (asymp_super_region = sum(adj_asymp_prop * population)/population_super_region,
                   mild_super_region = sum(adj_mild_prop * population)/population_super_region,
                   mod_super_region = sum(adj_mod_prop * population)/population_super_region,
                   sev_super_region = sum(adj_sev_prop*population)/population_super_region,
                   most_sev_super_region = sum(adj_most_sev_prop*population)/population_super_region), by = c("draw_num", "super_region_name")]

prevalence[, `:=` (asymp_global = sum(adj_asymp_prop * population)/population_global,
                   mild_global = sum(adj_mild_prop * population)/population_global,
                   mod_global = sum(adj_mod_prop * population)/population_global,
                   sev_global = sum(adj_sev_prop*population)/population_global,
                   most_sev_global =  sum(adj_most_sev_prop*population)/population_global), by = c("draw_num")]

pres_results_prop <- function(m){return(paste0(sprintf("%.1f", mean(m)*100), " (", sprintf("%.1f", quantile(m, 0.025)*100), "–", sprintf("%.1f", quantile(m, 0.975)*100), ")"))}

region_table_sev <- unique(prevalence[,.(location_name = region_name, Asymptomatic = pres_results_prop(asymp_region), Mild = pres_results_prop(mild_region), Moderate = pres_results_prop(mod_region), Severe = pres_results_prop(sev_region), Most_severe = pres_results_prop(most_sev_region)), by = "region_name"])

super_region_table_sev <- unique(prevalence[,.(location_name = super_region_name, Asymptomatic = pres_results_prop(asymp_super_region), Mild = pres_results_prop(mild_super_region), Moderate = pres_results_prop(mod_super_region), Severe = pres_results_prop(sev_super_region),  Most_severe = pres_results_prop(most_sev_super_region)), by = "super_region_name"])
global_table_sev <- unique(prevalence[,.(location_name = "Global", Asymptomatic = pres_results_prop(asymp_global), Mild = pres_results_prop(mild_global), Moderate = pres_results_prop(mod_global), Severe = pres_results_prop(sev_global),Most_severe = pres_results_prop(most_sev_global))])

region_table_sev <- merge(region_table_sev, locs[,.(location_name, lancet_label, sort_order)], by = "location_name")
super_region_table_sev <- merge(super_region_table_sev, locs[,.(location_name, lancet_label, sort_order)], by = "location_name")
global_table_sev <- merge(global_table_sev, locs[,.(location_name, lancet_label, sort_order)], by = "location_name")

country_table_sev <- unique(prevalence[,.(Asymptomatic = pres_results_prop(adj_asymp_prop), Mild = pres_results_prop(adj_mild_prop), Moderate = pres_results_prop(adj_mod_prop), Severe = pres_results_prop(adj_sev_prop), Most_severe = pres_results_prop(adj_most_sev_prop)), by = "location_id"])
country_table_sev <- merge(country_table_sev, locs[,.(location_id, lancet_label, sort_order)], by = "location_id")

final_table_sev <- rbind(country_table_sev, region_table_sev, super_region_table_sev, global_table_sev, fill = T)
final_table_sev <- final_table_sev[order(sort_order), .(Location = lancet_label, Asymptomatic, Mild, Moderate, Severe, Most_severe)]
final_table_sev[, duplicate := seq_len(.N), by = "Location"]
final_table_sev <- final_table_sev[duplicate == 1,]
final_table_sev[, duplicate := NULL]
final_table_sev[, `:=` (Asymptomatic = gsub("\\.", "·", Asymptomatic), Mild = gsub("\\.", "·", Mild), Moderate = gsub("\\.", "·", Moderate), Severe = gsub("\\.", "·", Severe), Most_severe = gsub("\\.", "·", Most_severe))]

write.xlsx(final_table_sev, paste0(output_dir, "appdendix_table_severity_all_",version_id,".xlsx"))
