###############################################################################-#
# giving the core codes of implementing the analysis                           -#
# Author: Wei Wang, Sicheng Li
# 2025-06-01                                                                   -#
###############################################################################-#
path_root <- "G:/crime study/pm2p5-crime"
setwd(path_root)
gcv_logit <- function(model){
  RSS <- 0-logLik(model)
  df <- attr(RSS,"df")
  as.numeric(model$nevent*RSS/(model$nevent-df)^2)
}
load("adata0.Rdata")
adata <- adata0

###### remove the outliers-----
max_lag <- 13
outlier <- quantile(adata$pm2p5_lag0,probs = c(0.01,0.99),na.rm = T)
aa <- adata[,"pm2p5_lag"%+%0:max_lag]
aa[aa < outlier[1] | aa > outlier[2]] <- NA
adata[,"pm2p5_lag"%+%0:max_lag] <- aa

###### DLNM--------
lag_day <- 13
lag_knots <- logknots(lag_day,nk = 1)
cb <- crossbasis(adata[,"pm2p5_lag"%+%0:lag_day],
                 argvar = list(fun = "lin"),
                 arglag = list(fun = "bs",degree = 2,knots = lag_knots)) 
fit_dlm <- clogit(event ~ cb + ns(temperature,3) + ns(wind,3) + ns(sunshine,3) + precipitation_cat + 
                     strata(id),data = adata)
BIC(fit_dlm)
gcv_logit(fit_dlm)

cen_pm2p5 <- 35
at_pm2p5 <- seq(ceiling(min(adata$pm2p5_lag0,na.rm = T)),
                      floor(max(adata$pm2p5_lag0,na.rm = T)),
                      by = 1)
tiff("results/lag effects at per 10 increase.tiff",width = 16,height = 8,units = "cm",pointsize = 11,res = 600)
opar <- par(mar = c(4, 4, 1, 1),mgp = c(2.5,1,0))
pred_dlm <- crosspred(cb,fit_dlm,at = at_pm2p5,bylag = 0.2,cumul = T,cen = cen_pm2p5)
label1 <- paste0("OR at lag0: ",round(pred_dlm$matRRfit["45","lag0"],4)," (95% CI: ",
                 round(pred_dlm$matRRlow["45","lag0"],4),"-",round(pred_dlm$matRRhigh["45","lag0"],4),")") 
label2 <- paste0("OR at lag1: ",round(pred_dlm$matRRfit["45","lag1"],4)," (95% CI: ",
                 round(pred_dlm$matRRlow["45","lag1"],4),"-",round(pred_dlm$matRRhigh["45","lag1"],4),")") 
plot.crosspred(pred_dlm, col = "#1C3D5A",
               ci.arg = list(col = "#D3D3D3"), 
               ylab = expression(paste("OR per ", 10, " ", mu, "g/", m^3," increase")),
               xlab = "Lag days", 
               lwd = 2,cumul = F,var = cen_pm2p5+10)
text(x = c(6,6),y = c(1.004,1.003),labels = c(label1,label2))
par(opar)
dev.off()


tiff("results/3D-effects-dlm.tiff",width = 16,height = 16,units = "cm",pointsize = 12,res = 600)
at_pm2p5 <- seq(0,floor(max(adata$pm2p5_lag0,na.rm = T)),by = 10)
pred_dlm <- crosspred(cb,fit_dlm,at = at_pm2p5,bylag = 1,cumul = T,cen = cen_pm2p5)
plot(pred_dlm, xlab="", main="",theta=120, phi=30,zlab = "",ylab = "Lag days")
text(x = 0.55,y = -0.55,labels = expression(PM[2.5]),srt = 50,cex = 1)  # x轴
text(x = -0.77,y = -0.1,labels = expression("OR referring to 35 "*mu*"g/"*m^3),srt = 95,cex = 1)  # z
dev.off()

###### the moving average exposure: nonlinear --------
type <- "lag01" 
adata$pm2p5 <- (adata$pm2p5_lag0 + adata$pm2p5_lag1)/2
vk_pm2p5 <- quantile(adata$pm2p5,probs = c(0.2,0.5,0.8),na.rm = T)
fit_linear <- clogit(event ~ pm2p5  + ns(temperature,3) + ns(wind,3) + 
                    ns(sunshine,3) + precipitation_cat + strata(id),data = adata)
BIC(fit_linear);gcv_logit(fit_linear)
fit_rcs <- clogit(event ~ rms::rcs(pm2p5,vk_pm2p5)  + ns(temperature,3) + ns(wind,3) + 
                     ns(sunshine,3) + precipitation_cat + strata(id),data = adata)
p_nonlinear <- anova(fit_linear,fit_rcs)[2,4] %>% {ifelse(.<0.001,"< 0.001",paste("=",round(.,3)))}
AIC(fit_rcs); BIC(fit_rcs);gcv_logit(fit_rcs)
at_pm2p5 <- seq(ceiling(min(adata$pm2p5,na.rm = T)),
                      floor(max(adata$pm2p5,na.rm = T)),by = 1)
cen_pm2p5 <- 35

index_aa <- grep("pm2p5",names(coef(fit_rcs)))
coef_aa <- coef(fit_rcs)[index_aa]
vcov_aa <- vcov(fit_rcs)[index_aa,index_aa]
basis_aa <- t(t(rms::rcs(at_pm2p5,vk_pm2p5)) - as.vector(rms::rcs(rep(cen_pm2p5,3),vk_pm2p5)[1,]))
risk <- drop(basis_aa %*% coef_aa) 
risk_se <- apply(basis_aa, 1, function(x) drop(x%*%vcov_aa%*%x)) %>% sqrt
lower_bound <- risk - 1.96 * risk_se
upper_bound <- risk + 1.96 * risk_se
interval <- round((max(upper_bound) - min(lower_bound))/8,2)
# plot
savename <- paste0("fit_rcs_at_",type) 
tiff(filename = paste0("results/",savename,".tiff"),
     width = 15,height = 15,units = "cm",pointsize = 12,res = 300)
g1 <- ggplot(data.frame(pm2p5 = at_pm2p5,risk = risk,lower = lower_bound,upper = upper_bound),
             aes(x = pm2p5, y = risk)) +
  geom_line() +  
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +  
  labs(x = "", y = "lnOR", title = NULL) + 
  theme_minimal() + 
  annotate("text", 
           x = mean(at_pm2p5),          
           y = min(lower_bound) + diff(range(c(lower_bound,upper_bound)))*0.1,           
           label = bquote(italic(P)[nonlinear]~.(p_nonlinear)),  
           size = 4,          
           color = "black") +
  theme(
    panel.grid = element_blank(),    
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),
    plot.margin = margin(t = 5, r = 5, b = -5, l = 5, unit = "pt")) +
  scale_x_continuous(
    breaks = c(seq(0, max(at_pm2p5), by = 25),15,35),  
    labels = c(seq(0, max(at_pm2p5), by = 25),15,35), 
    limits = range(adata$pm2p5,na.rm = T)
  ) +
  scale_y_continuous(
    breaks = seq(round(min(lower_bound)/interval)*interval, max(upper_bound), by = interval), 
    sec.axis = sec_axis(~ exp(.), name = "OR",
                        breaks = round(exp(seq(round(min(lower_bound)/interval)*interval, max(upper_bound), by = interval)),2))  # 右侧轴显示OR (exp(log(OR)))
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +  
  geom_vline(xintercept = c(15,35), linetype = "dashed", color = "grey")   

g2 <- ggplot(adata, aes(x = pm2p5)) +
  geom_histogram(aes(y = ..count..), bins = 30, fill = "gray70", color = "black") +
  labs(x = expression(PM[2.5]), y = "Frequency (K)") +
  scale_y_continuous(labels = function(x) x/1000) + 
  scale_x_continuous(limits = range(adata$pm2p5,na.rm = T)) + 
  theme_minimal()
cowplot::plot_grid(g1, g2, ncol = 1, align = "v", rel_heights = c(2, 1))
dev.off()
###### linear interaction terms --------

adata$pm2p5 <- (adata$pm2p5_lag0 + adata$pm2p5_lag1)/2
adata_merge <- merge(adata,data_crime_with_county[,-9],by = "id")
adata_merge$year <- year(as.Date(adata_merge$incident_date))

apply(adata_merge[,c("pm2p5","wind","sunshine","precipitation")], 2,quantile,probs = (0:10)/10,na.rm = T) 
cor(na.omit(adata_merge[,c("sunshine","precipitation","wind","temperature","pm2p5")]))
violence_keywords <- c("抢夺","抢劫", "故意杀人", "故意伤害", "强奸", "绑架", "虐待", "爆炸", "恐怖活动", 
                   "家暴", "殴打", "拘禁", "组织犯罪", "暴乱", "斗殴",
                   "威胁", "侵犯", "人身自由限制","暴力","寻衅滋事","猥亵")
larceny_keywords <- c("盗窃","偷盗")

adata_merge <- within(adata_merge,{
  temperature_group <- cut(temperature,breaks = c(-50,10,28,Inf))
  wind_group <- cut(wind,breaks = c(-1,1.5,3,Inf))
  sunshine_group <- cut(sunshine,breaks = c(-1,2,8,Inf))
  precipitation_group <- cut(precipitation,breaks = c(-1,0,5,Inf))
  violence_group <- grepl(paste(violence_keywords, collapse = "|"),x = case_type)
  larceny_group <- grepl(paste(larceny_keywords, collapse = "|"),x = case_type)
})
fit_0_group <- clogit(event ~ pm2p5 + ns(temperature,3) + ns(wind,3) + ns(sunshine,3) + 
                        precipitation_cat + strata(id),data = adata_merge)
aa <- summary(fit_0_group)$coefficients %>% .[grep("pm2p5",rownames(.)),] %>% t %>% as.data.frame
res_interaction <- data.frame(group = "All",group_value = "All",n_event = fit_0_group$nevent,aa,p_inter = NA) %>%
  {names(.) <- c("groupname","group_value","n_event","coef","exp_coef","se_coef","z","p_coef","p_inter");.}
for (i in c("temperature","wind","sunshine","precipitation","violence","larceny")) {
  adata_merge$group <- adata_merge[,i%+%"_group"]
  fit_int_group <- clogit(event ~ pm2p5:group - group  + ns(temperature,3) + 
                                    ns(wind,3) + ns(sunshine,3) + precipitation_cat + strata(id),data = adata_merge)
  p_int <- anova(fit_int_group,fit_0_group)[2,4]
  n_group <- table(adata_merge$group[adata_merge$type == 0 & !is.na(adata_merge$pm2p5)])
  coef_group<- summary(fit_int_group)$coefficients %>% .[grep("group",rownames(.)),]
  res_group <- data.frame(group = i,n_event = n_group,coef_group,p_inter = p_int) %>%
    {names(.) <- c("groupname","group_value","n_event","coef","exp_coef","se_coef","z","p_coef","p_inter");.}
  res_interaction <- rbind(res_interaction,res_group)
  cat(i,"\n")
}
write.csv(res_interaction,file = "results/res_interaction.csv")

res_interaction <- read.csv("results/res_interaction.csv")[,-1]
data_forest <- within(res_interaction,{
  group_value_new <- c("All","≤ 10","(10,28]","> 28","[0,1.5]","(1.5,3]","> 3",
                       "[0,2]","(2,8]","> 8","0","(0,5]","> 5","Non-violence",
                       "Violence","Non-larceny","Larceny")
  groupname <- c("All",rep(c("Temperature","Wind speed","Sunshine duration","Precipitation"),each = 3),
                 rep(c("Classification 1","Classification 2"),each = 2))
  p_inter <- ifelse(p_inter < 0.001,"< 0.001",sprintf("%.3f", p_inter))
  p_coef <- ifelse(p_coef < 0.001,"< 0.001",sprintf("%.3f", p_coef))
  OR <- exp(coef*10)
  OR_lower <- exp(coef*10 - 1.96*se_coef*10)
  OR_upper <- exp(coef*10 + 1.96*se_coef*10)
  CI <- paste0(sprintf("%.4f",OR)," (",sprintf("%.4f",OR_lower),"-",sprintf("%.4f",OR_upper),")")
}) 

varnames <- "All"
newdata <- data_forest[1,]
temp_data <- newdata
temp_data[,] <- NA
for (i in 2:nrow(data_forest)) {
  if(!data_forest$groupname[i] %in% varnames){
    temp_data$group_value_new <- data_forest$groupname[i]
    temp_data$p_inter <- data_forest$p_inter[i]
    newdata <- rbind(newdata,temp_data,data_forest[i,])
    varnames <- c(varnames,data_forest$groupname[i])
  } else {
    newdata <- rbind(newdata,data_forest[i,])
  }
}
newdata$p_inter[!is.na(newdata$coef)] <- NA

plotdata <- newdata
plotdata$is.bold <- F
plotdata$is.bold[is.na(plotdata$coef)] <- T
a1 <- as.list(plotdata$group_value_new)
a2 <- as.list(plotdata$n_event)
a3 <- as.list(plotdata$p_inter)
a4 <- as.list(plotdata$CI)
a5 <- as.list(plotdata$p_coef)
labeltext <- list(a1,a2,a3,a4,a5)

fp1 <- plotdata |>
  forestplot(
    labeltext = labeltext,
    mean = OR,
    lower = OR_lower,
    upper = OR_upper,
    is.summary = is.bold,
    col = fpColors(box = "#F8766D", line = "#F8766D"),
    boxsize = 0.2,
    hrzl_lines = list(`1` = gpar(lty = 1, lwd = 2),
                      `2` = gpar(lty = 2, lwd = 1.5),
                      `25` = gpar(lty = 1, lwd = 2)),
    graph.pos  = 4,
    xlab = expression(bold("OR per 10 μg/"*m^3*" increase")),
    grid = c(1.002,1.004,1.006,1.008),
    zero = 1,
    clip = c(0.999,1.008),
    xticks = c(1,1.002,1.004,1.006,1.008),
    align = c("l","c","c","c","c"), # 对齐方式
    xticks.digits = 3,
    colgap = unit(0.2,"cm"),
    txt_gp = fpTxtGp(summary = gpar(cex = 0.9),ticks = gpar(cex = 0.9),xlab = gpar(cex = 0.8)),  # 调整各部分字体
    lineheight = "auto",graphwidth = unit(3,"cm")) |>
  fp_add_header("Item",expression(bold(N[crime])),expression(bolditalic(P)[bold(interaction)]),
                "OR (95%CI)",expression(bolditalic(P)[bold(subgroup)])) 

tiff(filename = "results\\forest_subgroup.tiff",width = 15,height = 15,units = "cm",pointsize = 8,res = 600)
plot(fp1)
dev.off()



###### nonlinear interaction terms-----
covar <- "temperature"
vk_temperature <- quantile(adata_merge$temperature,probs = c(0.2,0.5,0.8),na.rm = T)
fit_linear_int <- clogit(event ~ pm2p5 + pm2p5:temperature + ns(temperature,3) + ns(wind,3) + ns(sunshine,3) + 
                              precipitation_cat + strata(id),data = adata_merge)
fit_nonlinear_int <- clogit(event ~ pm2p5 + pm2p5:rms::rcs(temperature,vk_temperature) + ns(temperature,3) + ns(wind,3) + ns(sunshine,3) + 
                        precipitation_cat + strata(id),data = adata_merge)
p_nonlinear <- anova(fit_nonlinear_int,fit_linear_int)[2,4] %>% {ifelse(.<0.001,"< 0.001",paste("=",round(.,3)))}
linear_int_res <- summary(fit_linear_int)$coefficients["pm2p5:temperature",]
linear_int_coef <- linear_int_res["coef"]*10*5
linear_int_coef_lower <- linear_int_coef - 1.96*linear_int_res["se(coef)"]*10*5
linear_int_coef_upper <- linear_int_coef + 1.96*linear_int_res["se(coef)"]*10*5
aa <- paste0(-round(linear_int_coef,4)," (95%CI:",round(-linear_int_coef_upper,4),"-",round(-linear_int_coef_lower,4),")")
int_coef_text <- bquote(atop(italic(P)[nonlinear]~.(p_nonlinear)~~phantom(),
                             atop("lnOR increase per 5"*degree*"C decrease: "~~phantom(),.(aa)))) 
at_temperature <- seq(ceiling(min(adata_merge$temperature,na.rm = T)),
                      floor(max(adata_merge$temperature,na.rm = T)),by = 1)
index_aa <- grep("pm2p5",names(coef(fit_nonlinear_int)))
coef_aa <- coef(fit_nonlinear_int)[index_aa]
vcov_aa <- vcov(fit_nonlinear_int)[index_aa,index_aa]
basis_aa <- cbind(1,rms::rcs(at_temperature,vk_temperature))
risk <- drop(basis_aa %*% coef_aa)*10
risk_se <- apply(basis_aa, 1, function(x) drop(x%*%vcov_aa%*%x)*100) %>% sqrt
lower_bound <- risk - 1.96 * risk_se
upper_bound <- risk + 1.96 * risk_se
interval <- round((max(upper_bound) - min(lower_bound))/6,3)

g1_temperature <- ggplot(data.frame(temperature = at_temperature,risk = risk,lower = lower_bound,upper = upper_bound),
             aes(x = temperature, y = risk)) +
  geom_line() +  
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +  
  labs(x = "", y = "lnOR", title = NULL) + 
  theme_minimal() + 
  annotate("text",
           x = min(at_temperature)+diff(range(at_temperature))/10*4,          
           y = min(lower_bound)+diff(range(c(lower_bound,upper_bound)))/10,           
           label = int_coef_text,  
           size = 3.5,         
           color = "black") +
  theme(
    panel.grid = element_blank(),    
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 10),
    plot.margin = margin(t = 5, r = 5, b = -5, l = 5, unit = "pt")) +
  scale_x_continuous(
    breaks = c(seq(ceiling(min(at_temperature)/10)*10,max(at_temperature), by = 10)),  
    labels = c(seq(ceiling(min(at_temperature)/10)*10,max(at_temperature), by = 10)), 
    limits = range(adata_merge$temperature,na.rm = T)
  ) +
  scale_y_continuous(
    breaks = seq(round(min(lower_bound)/interval)*interval, max(upper_bound), by = interval),  
    sec.axis = sec_axis(~ exp(.), name = "OR",
                        breaks = round(exp(seq(round(min(lower_bound)/interval)*interval, max(upper_bound), by = interval)),3))  # 右侧轴显示OR (exp(log(OR)))
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") 

g2_temperature <- ggplot(adata_merge, aes(x = temperature)) +
  geom_histogram(aes(y = ..count..), bins = 30, fill = "gray70", color = "black") +
  labs(x = expression("Temperature ("*degree*"C)"), y = "Frequency (K)") +
  scale_y_continuous(labels = function(x) x/1000) + 
  scale_x_continuous(limits = range(adata_merge$temperature,na.rm = T)) +
  theme_minimal() + 
  theme(axis.title = element_text(size = 10))
gg_temperature <- cowplot::plot_grid(g1_temperature, g2_temperature, ncol = 1, align = "v", rel_heights = c(2.5, 1))



covar <- "sunshine"
vk_sunshine <- quantile(adata_merge$sunshine,probs = c(0.2,0.5,0.8),na.rm = T)
fit_linear_int <- clogit(event ~ pm2p5 + pm2p5:sunshine + ns(temperature,3) + ns(wind,3) + ns(sunshine,3) + 
                           precipitation_cat + strata(id),data = adata_merge)
fit_nonlinear_int <- clogit(event ~ pm2p5 + pm2p5:rms::rcs(sunshine,vk_sunshine) + ns(temperature,3) + ns(wind,3) + ns(sunshine,3) + 
                              precipitation_cat + strata(id),data = adata_merge)
p_nonlinear <- anova(fit_nonlinear_int,fit_linear_int)[2,4] %>% {ifelse(.<0.001,"< 0.001",paste("=",round(.,3)))}
int_coef_text <- bquote(italic(P)[nonlinear]~.(p_nonlinear)) 
AIC(fit_nonlinear_int); BIC(fit_nonlinear_int);gcv_logit(fit_nonlinear_int)
at_sunshine <- 0:17
index_aa <- grep("pm2p5",names(coef(fit_nonlinear_int)))
coef_aa <- coef(fit_nonlinear_int)[index_aa]
vcov_aa <- vcov(fit_nonlinear_int)[index_aa,index_aa]
basis_aa <- cbind(1,rms::rcs(at_sunshine,vk_sunshine))
risk <- drop(basis_aa %*% coef_aa)*10
risk_se <- apply(basis_aa, 1, function(x) drop(x%*%vcov_aa%*%x)*100) %>% sqrt
lower_bound <- risk - 1.96 * risk_se
upper_bound <- risk + 1.96 * risk_se
interval <- round((max(upper_bound) - min(lower_bound))/6,3)

g1_sunshine <- ggplot(data.frame(sunshine = at_sunshine,risk = risk,lower = lower_bound,upper = upper_bound),
             aes(x = sunshine, y = risk)) +
  geom_line() +  
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +  
  labs(x = "", y = "lnOR", title = NULL) + 
  theme_minimal() + 
  annotate("text",
           x = min(at_sunshine)+diff(range(at_sunshine))/10*4,         
           y = min(lower_bound)-diff(range(c(lower_bound,upper_bound)))/10*0.6,          
           label = int_coef_text, 
           size = 3.5,         
           color = "black") +
  theme(
    panel.grid = element_blank(),   
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 10),
    plot.margin = margin(t = 5, r = 5, b = -5, l = 5, unit = "pt")) +
  scale_x_continuous(
    breaks = seq(0,17, by = 2),  
    labels = seq(0,17, by = 2),  
    limits = c(-0.5,17)
  ) +
  scale_y_continuous(
    breaks = seq(round(min(lower_bound)/interval)*interval, max(upper_bound), by = interval), 
    sec.axis = sec_axis(~ exp(.), name = "OR",
                        breaks = round(exp(seq(round(min(lower_bound)/interval)*interval, max(upper_bound), by = interval)),3))  # 右侧轴显示OR (exp(log(OR)))
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") 

g2_sunshine <- ggplot(adata_merge, aes(x = sunshine)) +
  geom_histogram(aes(y = ..count..), bins = 30, fill = "gray70", color = "black") +
  labs(x = "Sunshine duration (hours)", y = "Frequency (K)") +
  scale_y_continuous(labels = function(x) x/1000) + # 转换为千单位
  scale_x_continuous(limits = c(-0.5,17)) + # 与上图的x轴范围一致
  theme_minimal() + 
  theme(axis.title = element_text(size = 10))
gg_sunshine <- cowplot::plot_grid(g1_sunshine, g2_sunshine, ncol = 1, align = "v", rel_heights = c(2.5, 1))


tiff(filename = paste0("results/Figure 4.tiff"),
     width = 18,height = 9,units = "cm",pointsize = 8,res = 600)
cowplot::plot_grid(gg_temperature, gg_sunshine, ncol = 2, labels = c("a","b"),vjust = 1)
dev.off()


###### subgroup analysis across year ------
year_width <- 1
effect_years <- NULL
for (i in (2000+year_width):(2019-year_width)) {
  fit_year <- clogit(event ~ pm2p5 + ns(temperature,3) + ns(wind,3) + ns(sunshine,3) + precipitation_cat + strata(id),
                     data = adata_merge,subset = year >= i-year_width & year <= i + year_width)
  aa <- summary(fit_year)$coefficients[1,]
  aa <- c(aa,year = i,N = fit_year$nevent)
  effect_years <- rbind(effect_years,aa)
  cat(i," ")
}
effect_years <- as.data.frame(effect_years)
fit_year0 <- clogit(event ~ pm2p5 + ns(temperature,3) + ns(wind,3) + ns(sunshine,3) + precipitation_cat + strata(id),
                    data = adata_merge)
fit_year1 <- clogit(event ~ pm2p5 + pm2p5:year + ns(temperature,3) + ns(wind,3) + ns(sunshine,3) + 
                      precipitation_cat + strata(id),data = adata_merge)
p_linear <- anova(fit_year0,fit_year1)[2,4] %>% round(3) %>% paste("=",.)
fit_year2 <- clogit(event ~ pm2p5 + pm2p5:rms::rcs(year,3)+ ns(temperature,3) + ns(wind,3) + ns(sunshine,3) + 
                      precipitation_cat + strata(id),data = adata_merge)
p_nonlinear <- anova(fit_year0,fit_year2)[2,4] %>% round(3) %>% paste("=",.)
gg1 <- ggplot(effect_years,aes(x = year, y = exp(coef*10))) +
  geom_point(size = 0.8) + 
  geom_errorbar(aes(ymin = exp(coef*10 - 1.96*`se(coef)`*10), 
                    ymax = exp(coef*10 + 1.96*`se(coef)`*10)), width = 0.1, color = "black") + 
  labs(x = "year", y = expression(paste("OR per ", 10, " ", mu, "g/", m^3," increase")), title = NULL) +
  theme_minimal() + 
  scale_x_continuous(
    breaks = c(2004,2008,2012,2016),  # 设置横轴标签间隔为25
    labels = c(2004,2008,2012,2016)  # 添加35标签
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  theme(panel.grid = element_blank(),    # 去掉网格线
        axis.line = element_line(color = "black"),  # 添加坐标轴线
        axis.ticks = element_line(color = "black"))
tiff(filename = paste0("results/effect_over_year.tiff"),
     width = 15,height = 6.5,units = "cm",pointsize = 10,res = 300)
print(gg1)
dev.off()

###### subgroup analysis across provinces ----------
effect_regions <- NULL
region_names <- unique(adata_merge$province)
for (i in region_names) {
  if(sum(adata_merge$province == i,na.rm = T) < 100) next
  else {
    fit_region <- clogit(event ~ pm2p5 + ns(temperature,3) + ns(wind,3) + ns(sunshine,3) + precipitation_cat + strata(id),
                         data = adata_merge,subset = province == i)
    aa <- summary(fit_region)$coefficients[1,] %>% t %>% data.frame
  }
  aa <- data.frame(aa,region = i,N = fit_region$nevent)
  effect_regions <- rbind(effect_regions,aa)
  cat(i," ")
}
mapdata <- merge(mapdata_province,effect_regions,by.x = "ID",by.y = "region")

## because MCMAR showed an insignificant spatial dependence in the association, meta-analysis was used.
# nb <- poly2nb(mapdata)
# attr(nb,"region.id") <- mapdata$ID
# # 让广东与海南相邻
# nb[[9]] <- as.integer(6)
# nb[[6]] <- sort(c(nb[[6]],9)) %>% as.integer
# W <- nb2mat(nb)
# n <- ncol(W)
# C <- W * apply(W, 1, function(x) sum(x != 0))
# R <- diag(rowSums(C)) - C
# Cmatrix <- diag(n) - R
# ycoef <- cbind(mapdata$coef,mapdata$coef)
# scoef <- as.list(mapdata$se.coef.) %>% lapply(FUN = function(x) matrix(c(x^2,0,0,x^2),2,2))
# fit_mcmar <- mcmar(ycoef,S = scoef,Cmatrix = Cmatrix)
# fit_meta <- mvmeta(ycoef,S = scoef)
# llrtest(fit_mcmar,fit_meta)

fit_meta <- mvmeta(mapdata$coef,S = mapdata$se.coef.^2,method = "ml")
heter <- metafor::rma(mapdata$coef,sei = mapdata$se.coef)
qtest(fit_meta)
blue_meta <- blup.mvmeta(fit_meta,se = T)
mapdata$coef_meta <- blue_meta[,1]*10
mapdata$coef_meta_se <- blue_meta[,2]*10

tiff(filename = paste0("results/effect_over_province_meta.tiff"),
     width = 15,height = 6.5,units = "cm",pointsize = 8,res = 600)
gg <- ggplot(mapdata,aes(x = ID, y = exp(coef_meta))) +
  geom_point(size = 0.8) +  # 绘制预测的非线性曲线
  geom_errorbar(aes(ymin = exp(coef_meta - 1.96*coef_meta_se), 
                    ymax = exp(coef_meta + 1.96*coef_meta_se)), 
                width = 0.2, color = "black") +  # 添加95%置信区间
  labs(x = "Province", y = expression(paste("OR per ", 10, " ", mu, "g/", m^3," increase     ")), title = NULL) +
  theme_minimal() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 8),
        panel.grid = element_blank(),    # 去掉网格线
        axis.line = element_line(color = "black"),  # 添加坐标轴线
        axis.ticks = element_line(color = "black")) +  # 将横轴标签旋转45度
  coord_cartesian(ylim = c(0.997, 1.010)) 
print(gg)
dev.off()



###### attributable burdens ---------
fit_linear <- clogit(event ~ pm2p5 + ns(temperature,3) + ns(wind,3) + ns(sunshine,3) + precipitation_cat + strata(id),
                     data = adata_merge)
aa <- summary(fit_linear)$coefficients[1,]
a1 <- aa["coef"]; a2 <- aa["se(coef)"];
or_point <- a1
or_lower <- a1 - qnorm(0.975)*a2
or_upper <- a1 + qnorm(0.975)*a2

event_data <- within(adata0,{
  pm2p5 <- (pm2p5_dist_day0 + `pm2p5_dist_day-1`)/2
  year <- year(as.Date(incident_date))
}) %>% subset(event == 1,select = c("id","year","pm2p5")) %>%
  merge(data_crime_with_county,by = "id") %>%
  subset(select = c("id","province","city","year","pm2p5")) %>%
  na.omit()
rm(crime_data_pm2p5_60days)

refers <- c(0, 15, 35)
burden_data <- event_data %>% within({
  burden1_point <- (1-1/exp((pm2p5 - refers[1]) * or_point)) %>% pmax(0)
  burden1_lower <- (1-1/exp((pm2p5 - refers[1]) * or_lower)) %>% pmax(0)
  burden1_upper <- (1-1/exp((pm2p5 - refers[1]) * or_upper)) %>% pmax(0)
  burden2_point <- (1-1/exp((pm2p5 - refers[2]) * or_point)) %>% pmax(0)
  burden2_lower <- (1-1/exp((pm2p5 - refers[2]) * or_lower)) %>% pmax(0)
  burden2_upper <- (1-1/exp((pm2p5 - refers[2]) * or_upper)) %>% pmax(0)
  burden3_point <- (1-1/exp((pm2p5 - refers[3]) * or_point)) %>% pmax(0)
  burden3_lower <- (1-1/exp((pm2p5 - refers[3]) * or_lower)) %>% pmax(0)
  burden3_upper <- (1-1/exp((pm2p5 - refers[3]) * or_upper)) %>% pmax(0)
})

# for year-----
burden_data_year <- burden_data %>% group_by(year) %>%
  dplyr::summarize(N = length(burden1_point),
                   burden1_ratio = sum(burden1_point)/length(burden1_point)*100,
                   burden1_ratio_lower = sum(burden1_lower)/length(burden1_point)*100,
                   burden1_ratio_upper = sum(burden1_upper)/length(burden1_point)*100,
                   burden1_point = sum(burden1_point),
                   burden1_lower = sum(burden1_lower),
                   burden1_upper = sum(burden1_upper),
                   burden2_ratio = sum(burden2_point)/length(burden2_point)*100,
                   burden2_ratio_lower = sum(burden2_lower)/length(burden2_point)*100,
                   burden2_ratio_upper = sum(burden2_upper)/length(burden2_point)*100,
                   burden2_point = sum(burden2_point),
                   burden2_lower = sum(burden2_lower),
                   burden2_upper = sum(burden2_upper),
                   burden3_ratio = sum(burden3_point)/length(burden3_point)*100,
                   burden3_ratio_lower = sum(burden3_lower)/length(burden3_point)*100,
                   burden3_ratio_upper = sum(burden3_upper)/length(burden3_point)*100,
                   burden3_point = sum(burden3_point),
                   burden3_lower = sum(burden3_lower),
                   burden3_upper = sum(burden3_upper)) %>% ungroup()
 
cases_point <- apply(burden_data_year[,paste0("burden",1:3,"_point")], 2, sum) %>% round
cases_lower <- apply(burden_data_year[,paste0("burden",1:3,"_lower")], 2, sum) %>% round
cases_upper <- apply(burden_data_year[,paste0("burden",1:3,"_upper")], 2, sum) %>% round
cases_ratio <- (cases_point/sum(burden_data_year$N) * 100)  %>% round(2)
cases_ratio_lower <- (cases_lower/sum(burden_data_year$N) * 100) %>% round(2)
cases_ratio_upper <- (cases_upper/sum(burden_data_year$N) * 100)  %>% round(2)
burden_data_overall <- data.frame(type = rep(c("cases","ratio"),each = 3),
                 point = c(cases_point,cases_ratio),
                 lower = c(cases_lower,cases_ratio_lower),
                 upper = c(cases_upper,cases_ratio_upper))
burden_data_overall$value <- paste0(aa$point," (",aa$lower,"-",aa$upper,")") 
write.csv(burden_data_overall,file = "results/burden_data_overall.csv")

leaf_green <- "#E1CB55"  
leaf_yellow <- "#E07B54" 
leaf_blue <- "#51B1B7"  

jwidth <- 0.25
ratio_expand <- max(burden_data_year$burden1_upper)/max(burden_data_year$burden1_ratio_upper)
gg <- ggplot(burden_data_year) +
  
  geom_bar(aes(y = burden1_point, x = year - jwidth), stat = "identity", fill = leaf_blue, alpha = 0.6, width = 0.25, position = position_dodge(width = 0.3)) +
  geom_bar(aes(y = burden2_point, x = year), stat = "identity", fill = leaf_green, alpha = 0.6, width = 0.25, position = position_dodge(width = 0.3)) +
  geom_bar(aes(y = burden3_point, x = year + jwidth), stat = "identity", fill = leaf_yellow, alpha = 0.6, width = 0.25, position = position_dodge(width = 0.3)) +
  
  geom_line(aes(y = burden1_ratio*ratio_expand, x = year - jwidth), color = leaf_blue, size = 0.7) +
  geom_ribbon(aes(ymin = burden1_ratio_lower*ratio_expand, ymax = burden1_ratio_upper*ratio_expand, x = year - jwidth), 
              fill = leaf_blue, alpha = 0.2) +
  
  geom_line(aes(y = burden2_ratio*ratio_expand, x = year), color = leaf_green, size = 0.7) +
  geom_ribbon(aes(ymin = burden2_ratio_lower*ratio_expand, ymax = burden2_ratio_upper*ratio_expand, x = year), 
              fill = leaf_green, alpha = 0.2) +
  
  geom_line(aes(y = burden3_ratio*ratio_expand, x = year + jwidth), color = leaf_yellow, size = 0.7) +
  geom_ribbon(aes(ymin = burden3_ratio_lower*ratio_expand, ymax = burden3_ratio_upper*ratio_expand, x = year + jwidth),
              fill = leaf_yellow, alpha = 0.2) +
  labs(x = "year") +
  scale_y_continuous(
    name = "Attributable cases", 
    sec.axis = sec_axis(transform = ~ ./ratio_expand, name = "Attributable ratio (%)")  
  ) +
  theme_minimal(base_size = 25) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 16),
    legend.position.inside = c(0.9, 0.9), 
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 12),
    axis.text.x = element_text(angle = 30, hjust = 0.5,color = "black"),
    axis.text.y = element_text(angle = 90, vjust = 0.5,color = "black"),
    axis.title.y.left = element_text(margin = margin(r = 5)),  
    axis.title.y.right = element_text(margin = margin(l = 5)),  
    panel.grid = element_blank(),  
    axis.line.y.left = element_line(color = "black",size = 0.8), 
    axis.line.y.right = element_line(color = "black",size = 0.8),   
    axis.ticks.y.left = element_line(color = "black"), 
    axis.ticks.y.right = element_line(color = "black")  
  )
ggsave(paste0("results/burden_over_year.tiff"), 
       plot = gg, width = 15, height = 7, dpi = 600, device = "tiff")


# for province--------
burden_data_province <- burden_data %>% group_by(province) %>%
  dplyr::summarize(N = length(burden1_point),
                   burden1_ratio = sum(burden1_point)/length(burden1_point)*100,
                   burden1_ratio_lower = sum(burden1_lower)/length(burden1_point)*100,
                   burden1_ratio_upper = sum(burden1_upper)/length(burden1_point)*100,
                   burden1_point = sum(burden1_point),
                   burden1_lower = sum(burden1_lower),
                   burden1_upper = sum(burden1_upper),
                   burden2_ratio = sum(burden2_point)/length(burden2_point)*100,
                   burden2_ratio_lower = sum(burden2_lower)/length(burden2_point)*100,
                   burden2_ratio_upper = sum(burden2_upper)/length(burden2_point)*100,
                   burden2_point = sum(burden2_point),
                   burden2_lower = sum(burden2_lower),
                   burden2_upper = sum(burden2_upper),
                   burden3_ratio = sum(burden3_point)/length(burden3_point)*100,
                   burden3_ratio_lower = sum(burden3_lower)/length(burden3_point)*100,
                   burden3_ratio_upper = sum(burden3_upper)/length(burden3_point)*100,
                   burden3_point = sum(burden3_point),
                   burden3_lower = sum(burden3_lower),
                   burden3_upper = sum(burden3_upper)) %>% ungroup()

mapdata_burden <- merge(mapdata_province,burden_data_province,by.x = "ID",by.y = "province",all = T)
para_list <- list(ref_nums = c(1,1,2,2,3,3),
                  burden_types = c("case","ratio","case","ratio","case","ratio"))
plotlist <- NULL
for (i in 1:length(para_list[[1]])) {
  ref_num <- para_list$ref_nums[i]
  burden_type <- para_list$burden_types[i]
  if(burden_type == "case"){
    valuename <- paste0("burden",ref_num,"_point")
    legendname <- "Attributable cases"
    highcolor <- "#D62728"
  } else{
    valuename <- paste0("burden",ref_num,"_ratio")
    legendname <- "Attributable ratio (%)"
    highcolor <- "blue"
  }
  ref <- refers[ref_num]
  insert_map <- ggplot(mapdata_burden) +
    geom_sf(aes(fill = get(valuename)),color = "white") +
    scale_fill_gradient(low = "lightgray", high = highcolor) + 
    geom_sf(data = mapdata_iland,color = "lightgrey", fill = "white") +
    coord_sf(xlim = c(105, 123), ylim = c(2, 25)) +  
    theme_void() +  
    theme(
      legend.position = "none", 
      plot.background = element_rect(color = "black", size = 1) 
    )
  main_map <- ggplot(mapdata_burden) +
    geom_sf(aes(fill = get(valuename)),color = "white") +
    scale_fill_gradient(low = "lightgray", high = highcolor, 
                        name = bquote(atop(.(legendname), 
                                           paste("[Refer to ",.(ref)," ",mu,"g/",m^3,"]")))) + 
    geom_sf(data = mapdata_iland,color = "lightgrey", fill = "white") +
    geom_text(aes(x = X, y = Y,label = ID), color = "black", size = 3.5) +  
    theme_void(base_size = 14) + 
    labs(title = "") +
    theme(
      legend.position = c(0.91, 0.41),  
      legend.title = element_text(hjust = 0.5), 
      plot.margin = margin(t = -2, r = 1, b = -2, l = 0, unit = "cm"),
      plot.background = element_blank()
    ) + 
    coord_sf(xlim = c(75,NA),ylim = c(18.5, 53))  
  gg <- ggdraw() +
    draw_plot(main_map) +
    draw_plot(insert_map, x = 0.8, y = -0.155, width = 0.1, height = 0.5)  
    plotlist <- c(plotlist,list(gg))
  cat(i," ")
}
tiff(filename = paste0("results/burden_province.tiff"),width = 18*2,height = 19*2,units = "cm",pointsize = 8,res = 600)
plot_grid(plotlist = plotlist,nrow = 3,align = "h",labels = c("c","d","a","b","c","d"),label_size = 22,vjust = 2)
dev.off()

