library(tidyverse)
library(survival)
library(rms)
library(dlnm)
library(MCMAR)
library(spdep)
library(mvmeta)
library(cowplot)
library(patchwork)

`%+%` <- function(x,y) paste0(x,y)
maxlag_day <- 13
gcv_logit <- function(model){
  RSS <- logLik(model)
  df <- attr(RSS,"df")
  as.numeric(model$nevent*RSS/(model$nevent-df)^2)
}

####### 0.data prepare----------------------------------------------
data0 <- load("adata_crime_PM2.5.Rdata") %>% get
# remove outliers
aa <- data0[,grep("pm2p5_lag",colnames(data0))] %>% unlist %>% quantile(na.rm = T,probs = c(0.01,0.99))
bb <- data0[,"pm2p5_lag"%+%0:maxlag_day]
data0[,"pm2p5_lag"%+%0:maxlag_day][bb < aa[1] | bb > aa[2]] <- NA

data0 <- within(data0,{
  temperature <- temperature_lag0
  temperature[is.na(temperature)] <- mean(temperature_lag0,na.rm = T)
  precipitation <- precipitation_lag0
  precipitation[is.na(precipitation)] <- mean(precipitation_lag0,na.rm = T)
  sunshine <- sunshine_lag0
  sunshine[is.na(sunshine)] <- mean(sunshine_lag0,na.rm = T)
  precipitation_cat <- precipitation
  precipitation_cat[precipitation_cat > 0] <- 1
  y <- ifelse(type == 0,TRUE,FALSE)
})

####### 1. using DLM --------------------------------------------
lag_day <- maxlag_day
narow <- apply(data0[,"pm2p5_lag"%+%0:lag_day],1,anyNA)
adata <- subset(data0,!narow) %>%
  group_by(id) %>%
  filter(any(type == 0)) %>%
  ungroup()   
adatax <- adata
vk_temperature <- quantile(adatax$temperature,probs = c(0.1,0.5,0.9))
vk_sunshine <- quantile(adatax$sunshine,probs = c(0.1,0.5,0.9))
vk_precipitation <- quantile(adatax$precipitation[adatax$precipitation != 0],probs = c(0.1,0.5,0.9))

##  selecting knots for the lag dimension 
lagknots_bic <- NULL
for(i in 1:4){
  lag_knots <- logknots(lag_day,nk = i)
  cb <- crossbasis(adatax[,"pm2p5_lag"%+%0:lag_day],argvar = list(fun = "lin"),
                   arglag = list(fun = "bs",degree = 2,knots = lag_knots)) 
  fit_dlm <- clogit(y ~ cb + rcs(temperature,vk_temperature) + 
                      rcs(sunshine,vk_sunshine) + rcs(precipitation,vk_precipitation) + strata(id),data = adatax)
  lagknots_bic <- c(lagknots_bic,BIC(fit_dlm))
  cat(i," ")
}
tiff(filename = "select_number_of_lag.tiff",width = 10,height = 8,
     units = "cm",pointsize = 10,res = 300)
plot(1:4, lagknots_bic, type = "l", xlab = "Number of knots", ylab = "BIC", xaxt = "n")
axis(1, at = 1:4, labels = 1:4)
points(1:4,lagknots_bic,pch = 16,col = c("red","black","black","black"))
dev.off()

# using the optimal knots to fit DLM
lag_knots <- logknots(maxlag_day,nk = 1)
cb <- crossbasis(adatax[,"pm2p5_lag"%+%0:lag_day],argvar = list(fun = "lin"),
                 arglag = list(fun = "bs",degree = 2,knots = lag_knots)) 
fit_dlm <- clogit(y ~ cb + rcs(temperature,vk_temperature) + 
                    rcs(sunshine,vk_sunshine) + rcs(precipitation,vk_precipitation) + strata(id),data = adatax)
BIC(fit_dlm)
concordance(fit_dlm)
gcv_logit(fit_dlm)

pred_dlm <- crosspred(cb,fit_dlm,at = 45,bylag = 0.2,cumul = T,cen = 35)
tiff("lag-effect-13days.tiff",width = 16,height = 8,units = "cm",pointsize = 10,res = 300)
opar <- par(mar = c(4, 4, 1, 1),mgp = c(2.5,1,0))
plot.crosspred(pred_dlm, "slices", var = 45, col = "#1C3D5A",             
               ci.arg = list(col = "#D3D3D3"), 
               ylab = expression(paste("OR per ", 10, " ", mu, "g/", m^3," increase")) ,
               xlab = "Lag days", 
               lwd = 2,cumul = F)  
par(opar)
dev.off()



####### 2. analysis using the MAE----------------------
## 2.0 data processing
lag_day <- 1
narow <- apply(data0[,"pm2p5_lag"%+%0:lag_day],1,anyNA)
adata <- subset(data0,!narow) %>%
  group_by(id) %>%
  filter(any(type == 0)) %>%
  ungroup()   
adatax <- adata
adatax$average_pm2p5 <- apply(adatax[,"pm2p5_lag"%+%0:lag_day],1,mean)
adatax$average_pm2p5_cat <- adatax$average_pm2p5 > 35
vk_temperature <- quantile(adatax$temperature,probs = c(0.1,0.5,0.9))
vk_sunshine <- quantile(adatax$sunshine,probs = c(0.1,0.5,0.9))
vk_precipitation <- quantile(adatax$precipitation[adatax$precipitation != 0],probs = c(0.1,0.5,0.9))

## 2.1 linear ----
type_keywords <- c("抢夺","抢劫", "故意杀人", "故意伤害", "强奸", "绑架", "虐待", "爆炸", "恐怖活动", 
                       "家暴", "殴打", "拘禁", "组织犯罪", "暴乱", "斗殴",
                       "威胁", "侵犯", "人身自由限制","暴力","寻衅滋事","猥亵")
adatax$violence <- grepl(paste(type_keywords, collapse = "|"),x = adata$case_type)
xx1 <- table(subset(adatax,adatax$violence & type == 0,case_type)) %>% as.data.frame
sum(xx1[,2])
xx2 <- table(subset(adatax,!adatax$violence & type == 0,case_type)) %>% as.data.frame

type_keywords <- c("盗窃","偷盗")
adatax$larceny <- grepl(paste(type_keywords, collapse = "|"),x = adata$case_type)
xx1 <- table(subset(adatax,adatax$larceny & type == 0,case_type)) %>% as.data.frame
sum(xx1[,2])
xx2 <- table(subset(adatax,!adatax$larceny & type == 0,case_type)) %>% as.data.frame

linear_effects <- NULL
crime_types <- c("Overall","Violence","Non-violence","Larceny","Non-larceny")
for(i in 1:length(crime_types)){
  crime_type <- crime_types[i]
  if(crime_type == "Overall") data_temp <- adatax
  else if(crime_type == "Violence") data_temp <- subset(adatax,violence == T)
  else if(crime_type == "Non-violence") data_temp <- subset(adatax,violence == F)
  else if(crime_type == "Larceny") data_temp <- subset(adatax,larceny == T)
  else if(crime_type == "Non-larceny") data_temp <- subset(adatax,larceny == F)
  else stop("Incorrect crime type")
  fit_linear0 <- clogit(y ~ average_pm2p5 + strata(id),data = data_temp)
  aa <- summary(fit_linear0)$coefficients[1,]
  a1 <- aa["coef"]*10; a2 <- aa["se(coef)"]*10;
  a3 <- ifelse(aa["Pr(>|z|)"]<0.001,"< 0.001",round(aa["Pr(>|z|)"],3))
  a4 <- qnorm(0.975)
  bb1 <- data.frame(model_type = "model0",OR_type = "continuous",
                    type = crime_type,OR = round(exp(a1),4),
                    lower = round(exp(a1 - a4*a2),4),
                    upper = round(exp(a1 + a4*a2),4),
                    p_value = a3,n_event = fit_linear0$nevent)
  
  fit_linear1 <- clogit(y ~ average_pm2p5 + rcs(temperature,vk_temperature) + 
                          rcs(sunshine,vk_sunshine) + rcs(precipitation,vk_precipitation) + 
                          strata(id),data = data_temp)
  aa <- summary(fit_linear1)$coefficients[1,]
  a1 <- aa["coef"]*10; a2 <- aa["se(coef)"]*10;
  a3 <- ifelse(aa["Pr(>|z|)"]<0.001,"< 0.001",round(aa["Pr(>|z|)"],3))
  a4 <- qnorm(0.975)
  bb2 <- data.frame(model_type = "model1",OR_type = "continuous",
                    type = crime_type,OR = round(exp(a1),4),
                    lower = round(exp(a1 - a4*a2),4),
                    upper = round(exp(a1 + a4*a2),4),
                    p_value = a3,n_event = fit_linear1$nevent)
  
  fit_linear0 <- clogit(y ~ average_pm2p5_cat + strata(id),data = data_temp)
  aa <- summary(fit_linear0)$coefficients[1,]
  a1 <- aa["coef"]; a2 <- aa["se(coef)"];
  a3 <- ifelse(aa["Pr(>|z|)"]<0.001,"<0.001",sprintf(aa["Pr(>|z|)"], fmt = paste0("%#.",3, "f")))
  a4 <- qnorm(0.975)
  bb3 <- data.frame(model_type = "model0",OR_type = "categorical",
                    type = crime_type,OR = round(exp(a1),4),
                    lower = round(exp(a1 - a4*a2),4),
                    upper = round(exp(a1 + a4*a2),4),
                    p_value = a3,n_event = fit_linear0$nevent)
  fit_linear1 <- clogit(y ~ average_pm2p5_cat + rcs(temperature,vk_temperature) + 
                          rcs(sunshine,vk_sunshine) + rcs(precipitation,vk_precipitation) + 
                          strata(id),data = data_temp)
  aa <- summary(fit_linear1)$coefficients[1,]
  a1 <- aa["coef"]; a2 <- aa["se(coef)"];
  a3 <- ifelse(aa["Pr(>|z|)"]<0.001,"< 0.001",round(aa["Pr(>|z|)"],3))
  a4 <- qnorm(0.975)
  bb4 <- data.frame(model_type = "model1",OR_type = "categorical",
                    type = crime_type,OR = round(exp(a1),4),
                    lower = round(exp(a1 - a4*a2),4),
                    upper = round(exp(a1 + a4*a2),4),
                    p_value = a3,n_event = fit_linear1$nevent)
  linear_effects <- rbind(linear_effects,bb1,bb2,bb3,bb4)
  cat(i," ")
}
linear_effects <- linear_effects[order(linear_effects$model_type),]
write.csv(linear_effects,file = "linear_effects.csv")

## 2.2 rcs ---------------
lag_day <- 1
adatax$average_pm2p5 <- apply(adatax[,"pm2p5_lag"%+%0:lag_day],1,mean)
vk_pm2p5 <- quantile(adatax$average_pm2p5,probs = c(0.1,0.5,0.9))
# precipitation_cat has smaller BIC than rcs(precipitation,vk_precipitation)
fit_rcs <- clogit(y ~ rcs(average_pm2p5,vk_pm2p5) + rcs(temperature,vk_temperature) + 
                      rcs(sunshine,vk_sunshine) + precipitation_cat + #rcs(precipitation,vk_precipitation) + 
                      strata(id),data = adatax)
AIC(fit_rcs); BIC(fit_rcs)
ap_seq <- seq(min(adatax$average_pm2p5), max(adatax$average_pm2p5), length.out = 100)

index_aa <- grep("average_pm2p5",names(coef(fit_rcs)))
coef_aa <- coef(fit_rcs)[index_aa]
vcov_aa <- vcov(fit_rcs)[index_aa,index_aa]
basis_aa <- t(t(rcs(ap_seq,vk_pm2p5)) - as.vector(rcs(rep(35,3),vk_pm2p5)[1,]))
risk <- drop(basis_aa %*% coef_aa) 
risk_se <- apply(basis_aa, 1, function(x) drop(x%*%vcov_aa%*%x)) %>% sqrt
lower_bound <- risk - 1.96 * risk_se
upper_bound <- risk + 1.96 * risk_se

tiff(filename = paste0("nonlinear_lag01.tiff"),
     width = 15,height = 8,units = "cm",pointsize = 12,res = 300)
gg <- ggplot(data.frame(average_pm2p5 = ap_seq,risk = risk,lower = lower_bound,upper = upper_bound),
             aes(x = average_pm2p5, y = risk)) +
  geom_line() +  
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +  
  labs(x = expression(paste(PM[2.5]," (",mu,"g/",m^3,")")), y = "log(OR)", title = NULL) + 
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),    
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black")) +
  scale_x_continuous(
    breaks = c(seq(0, max(ap_seq), by = 25),35),  
    labels = c(seq(0, max(ap_seq), by = 25), 35)  
  ) +
  scale_y_continuous(
    breaks = seq(-0.04, max(upper_bound), by = 0.02),  
    sec.axis = sec_axis(~ exp(.), name = "OR",
                        breaks = round(exp(seq(-0.04, max(upper_bound), by = 0.02)),2))   
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +   
  geom_vline(xintercept = 35, linetype = "dashed", color = "black")    
  print(gg)
dev.off()


fit_linear <- clogit(y ~ average_pm2p5 + rcs(temperature,vk_temperature) + 
                    rcs(sunshine,vk_sunshine) + precipitation_cat + #rcs(precipitation,vk_precipitation) + 
                    strata(id),data = adatax)
anova(fit_linear,fit_rcs)


## 2.3 comparison across years ------
lag_day <- 1
adatax$average_pm2p5 <- apply(adatax[,"pm2p5_lag"%+%0:lag_day],1,mean)
year_width <- 1
effect_years <- NULL
for (i in (2000+year_width):(2019-year_width)) {
    fit_year <- clogit(y ~ average_pm2p5 + rcs(temperature,vk_temperature) + 
                       rcs(sunshine,vk_sunshine) + precipitation_cat + #rcs(precipitation,vk_precipitation) + 
                       strata(id),data = adatax,subset = year >= i-year_width & year <= i + year_width)
    aa <- summary(fit_year)$coefficients[1,]
    aa <- c(aa,year = i,N = fit_year$nevent)
    effect_years <- rbind(effect_years,aa)
    cat(i," ")
}
effect_years <- as.data.frame(effect_years)

tiff(filename = paste0("results/effect_over_year.tiff"),
     width = 15,height = 6.5,units = "cm",pointsize = 10,res = 300)
gg <- ggplot(effect_years,aes(x = year, y = exp(coef))) +
  geom_point(size = 0.8) +   
  geom_errorbar(aes(ymin = exp(coef - 1.96*`se(coef)`), 
                    ymax = exp(coef + 1.96*`se(coef)`)), width = 0.1, color = "black") +   
  labs(x = "year", y = "OR", title = NULL) +
  theme_minimal() + 
  scale_x_continuous(
    breaks = c(2004,2008,2012,2016),   
    labels = c(2004,2008,2012,2016)   
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  theme(panel.grid = element_blank(),     
        axis.line = element_line(color = "black"),   
        axis.ticks = element_line(color = "black"))
print(gg)
dev.off()

fit_year0 <- clogit(y ~ average_pm2p5 + rcs(temperature,vk_temperature) + 
                      rcs(sunshine,vk_sunshine) + precipitation_cat + #rcs(precipitation,vk_precipitation) + 
                      strata(id),data = adatax)
fit_year1 <- clogit(y ~ average_pm2p5 + average_pm2p5:year + rcs(temperature,vk_temperature) + 
                      rcs(sunshine,vk_sunshine) + precipitation_cat + #rcs(precipitation,vk_precipitation) + 
                      strata(id),data = adatax)
anova(fit_year0,fit_year1)

fit_year2 <- clogit(y ~ average_pm2p5 + average_pm2p5:rcs(year,3) + rcs(temperature,vk_temperature) + 
                      rcs(sunshine,vk_sunshine) + precipitation_cat + #rcs(precipitation,vk_precipitation) + 
                      strata(id),data = adatax)
anova(fit_year0,fit_year2)
## 2.4 comparison across provinces -----------
lag_day <- 1
adatax$average_pm2p5 <- apply(adatax[,"pm2p5_lag"%+%0:lag_day],1,mean)
effect_regions <- NULL
region_names <- unique(adatax$province)
for (i in region_names) {
  if(sum(adatax$province == i,na.rm = T) < 100) next
  else {
    fit_region <- clogit(y ~ average_pm2p5 + rcs(temperature,vk_temperature) + 
                         rcs(sunshine,vk_sunshine) + precipitation_cat + #rcs(precipitation,vk_precipitation) + 
                         strata(id),data = adatax,subset = province == i)
    aa <- summary(fit_region)$coefficients[1,] %>% t %>% data.frame
  }
  aa <- data.frame(aa,region = i,N = fit_region$nevent)
  effect_regions <- rbind(effect_regions,aa)
  cat(i," ")
}
effect_regions <- within(effect_regions,{
  region[region == "Neimenggu"] <- "NeiMongol"
  region[region == "Shanxi"] <- "Shaanxi"
  region[region == "Shan'xi"] <- "Shanxi"
}) 

mapdata <- merge(mapdata_province,effect_regions,by.x = "ID",by.y = "region")
## because MCMAR showed no spatial dependence, we used the classic meta analysis
# nb <- poly2nb(mapdata)
# attr(nb,"region.id") <- mapdata$ID
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
fit_meta <- mvmeta(mapdata$coef,S = mapdata$se.coef.^2,method = "reml")
heter <- metafor::rma(mapdata$coef,sei = mapdata$se.coef)
qtest(fit_meta)
blue_meta <- blup.mvmeta(fit_meta,se = T)
mapdata$coef_meta <- blue_meta[,1]*10
mapdata$coef_meta_se <- blue_meta[,2]*10

tiff(filename = paste0("effect_over_province_meta.tiff"),
     width = 15,height = 6.5,units = "cm",pointsize = 10,res = 300)
gg <- ggplot(mapdata,aes(x = ID, y = exp(coef_meta))) +
  geom_point(size = 0.8) +  
  geom_errorbar(aes(ymin = exp(coef_meta - 1.96*coef_meta_se), 
                    ymax = exp(coef_meta + 1.96*coef_meta_se)), 
                width = 0.2, color = "black") +  
  labs(x = "Province", y = "OR", title = NULL) +
  theme_minimal() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 8),
        panel.grid = element_blank(),    
        axis.line = element_line(color = "black"),  
        axis.ticks = element_line(color = "black")) +  
  coord_cartesian(ylim = c(0.998, 1.011)) 
print(gg)
dev.off()


## 2.5 comparison across cities -----------
effect_city <- NULL
adatax$city_code <- substr(adatax$county_code,1,4)
region_names <- unique(adatax$city_code)
for (i in region_names) {
  if(sum(adatax$city_code == i,na.rm = T) < 100) next
  else {
    fit_region <- clogit(y ~ average_pm2p5 + rcs(temperature,vk_temperature) + 
                           rcs(sunshine,vk_sunshine) + precipitation_cat + #rcs(precipitation,vk_precipitation) + 
                           strata(id),data = adatax,subset = city_code == i)
    aa <- summary(fit_region)$coefficients[1,] %>% t %>% data.frame
  }
  aa <- data.frame(aa,region = i,N = fit_region$nevent)
  effect_city <- rbind(effect_city,aa)
  cat(i," ")
}
mapdata_city$city_code <- substr(mapdata_city$code,1,4)
mapdata <- merge(mapdata_city,effect_city,by.x = "city_code",by.y = "region")
fit_meta <- mvmeta(mapdata$coef,S = mapdata$se.coef.^2,method = "reml")
qtest(fit_meta)
blue_meta <- blup.mvmeta(fit_meta,se = T)
mapdata$coef_meta <- blue_meta[,1]
mapdata$coef_meta_se <- blue_meta[,2]
mapdata$p_value <- mapdata$Pr...z.. 
metafor::rma(yi = mapdata$coef,sei = mapdata$se.coef.)

######### 3. attributable burdens -----------
fit_linear <- clogit(y ~ average_pm2p5 + rcs(temperature,vk_temperature) + 
                        rcs(sunshine,vk_sunshine) + rcs(precipitation,vk_precipitation) + 
                        strata(id),data = adatax)
aa <- summary(fit_linear)$coefficients[1,]
a1 <- aa["coef"]; a2 <- aa["se(coef)"];
or_point <- a1
or_lower <- a1 - qnorm(0.975)*a2
or_upper <- a1 + qnorm(0.975)*a2

event_data <- load("adata_crime_PM2.5.Rdata") %>% get
rm(adata_crime)
event_data <- subset(event_data,type == 0) %>% 
  within(average_pm2p5 <- (pm2p5_lag0 + pm2p5_lag1)/2) %>%
  subset(select = c("id","province","city","year","average_pm2p5")) %>%
  na.omit()

refers <- c(0, 15, 35)
burden_data <- event_data %>% within({
  burden1_point <- (1-1/exp((average_pm2p5 - refers[1]) * or_point)) %>% pmax(0)
  burden1_lower <- (1-1/exp((average_pm2p5 - refers[1]) * or_lower)) %>% pmax(0)
  burden1_upper <- (1-1/exp((average_pm2p5 - refers[1]) * or_upper)) %>% pmax(0)
  burden2_point <- (1-1/exp((average_pm2p5 - refers[2]) * or_point)) %>% pmax(0)
  burden2_lower <- (1-1/exp((average_pm2p5 - refers[2]) * or_lower)) %>% pmax(0)
  burden2_upper <- (1-1/exp((average_pm2p5 - refers[2]) * or_upper)) %>% pmax(0)
  burden3_point <- (1-1/exp((average_pm2p5 - refers[3]) * or_point)) %>% pmax(0)
  burden3_lower <- (1-1/exp((average_pm2p5 - refers[3]) * or_lower)) %>% pmax(0)
  burden3_upper <- (1-1/exp((average_pm2p5 - refers[3]) * or_upper)) %>% pmax(0)
})

### 3.1 year--------
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
write.csv(burden_data_year,file = "results/burden_data_year.csv")

# 
cases_point <- apply(burden_data_year[,paste0("burden",1:3,"_point")], 2, sum) %>% round
cases_lower <- apply(burden_data_year[,paste0("burden",1:3,"_lower")], 2, sum) %>% round
cases_upper <- apply(burden_data_year[,paste0("burden",1:3,"_upper")], 2, sum) %>% round
cases_ratio <- (cases_point/sum(burden_data_year$N) * 100)  %>% round(2)
cases_ratio_lower <- (cases_lower/sum(burden_data_year$N) * 100) %>% round(2)
cases_ratio_upper <- (cases_upper/sum(burden_data_year$N) * 100)  %>% round(2)
aa <- data.frame(type = rep(c("cases","ratio"),each = 3),
                 point = c(cases_point,cases_ratio),
                 lower = c(cases_lower,cases_ratio_lower),
                 upper = c(cases_upper,cases_ratio_upper))
aa$value <- paste0(aa$point," (",aa$lower,"-",aa$upper,")") 
burden_data_overall <- write.csv(aa,file = "results/burden_data_overall.csv")
# plot 
leaf_green <- "#E1CB55"  # 深绿
leaf_yellow <- "#E07B54" # 柳叶刀黄色
leaf_blue <- "#51B1B7"  # 蓝色

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
ggsave(paste0("results/burden_over_year00.tiff"), 
       plot = gg, width = 15, height = 7, dpi = 300, device = "tiff")


### 3.2 province--------
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
write.csv(burden_data_province,file = "burden_data_procince.csv")

burden_data_province <- within(burden_data_province,{
  province[province == "Neimenggu"] <- "NeiMongol"
  province[province == "Shanxi"] <- "Shaanxi"
  province[province == "Shan'xi"] <- "Shanxi"
}) 

mapdata_burden <- merge(mapdata_province,burden_data_province,by.x = "ID",by.y = "province",all = T)
ref_num <- 1; burden_type <- "ratio" # ratio case
if(burden_type == "case"){
  valuename <- paste0("burden",ref_num,"_point")
  legendname <- "Attributable cases"
} else{
  valuename <- paste0("burden",ref_num,"_ratio")
  legendname <- "Attributable ratio (%)"
}
ref <- refers[ref_num]
insert_map <- ggplot(mapdata_burden) +
  geom_sf(aes(fill = get(valuename)),color = "white") +
  scale_fill_gradient(low = "lightgray", high = "#D62728") + 
  geom_sf(data = mapdata_iland,color = "lightgrey", fill = "white") +
  
  coord_sf(xlim = c(105, 123), ylim = c(2, 25)) +  
  theme_void() +  
  theme(
    legend.position = "none",  
    plot.background = element_rect(color = "black", size = 1)  
  )
main_map <- ggplot(mapdata_burden) +
  geom_sf(aes(fill = get(valuename)),color = "white") +
  scale_fill_gradient(low = "lightgray", high = "#D62728", 
                      name = bquote(atop(.(legendname), 
                                             paste("[Refer to ",.(ref)," ",mu,"g/",m^3,"]")))) + 
  geom_sf(data = mapdata_iland,color = "lightgrey", fill = "white") +
  geom_text(aes(x = X, y = Y,label = ID), color = "black", size = 4) +  
  theme_void(base_size = 14) +  
  labs(title = "") +
  theme(
    legend.position = c(0.88, 0.48),  
    legend.title = element_text(hjust = 0.5), 
    plot.background = element_blank()
  ) + 
  
  coord_sf(xlim = c(75,NA),ylim = c(17, 53))  
map_province <- ggdraw() +
  draw_plot(main_map) +
  draw_plot(insert_map, x = 0.8, y = -0.13, width = 0.12, height = 0.6)  
ggsave(paste0("map_burden_province_",burden_type,ref,".tiff"), 
       plot = map_province, width = 10, height = 8, dpi = 300, device = "tiff")


