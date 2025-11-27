#####################################################################################################################
### Low-Level Airborne Particulate Matter and Risk of Hypertension Hospitalization in Older U.S. Adults
### Main analysis code
### Yechi Zhang, et al. Yale School of Public Health, November 27, 2025
#####################################################################################################################
rm(list=ls())

library(data.table); library(dplyr); library(ggplot2); library(gridExtra); library(mgcv); library(parallel)
library(doParallel); library(xgboost); library(openxlsx)

# Load the data
data <- fread("Medicare_Poisson_FFS_MA.csv")
data <- as.data.frame(data)

covariates <- c( "NDVI_SUMMER","NatWalkInd", "RUCA1",
                 "PCT_AGE65", "PCT_AGE65_POVERTY", "PCT_AGE65_LOWEDU", 
                 "PCT_MINORITY",  "HOUSE_NO_INTERNET", "HOUSE_CROWD", "HOUSE_COST", 
                 "UNEMPLOY", "POP_DENS",  "DRINK",  "SMOKE",  "OBESITY", 
                 "PHYSINACT", "ANNUAL_CHECK", "CHOLSCREEN", "TMEAN_WARM", "RH_WARM", "TMEAN_COLD", "RH_COLD")
data <- data[, c("ZIP5", "YEAR", "AGE_CAT","BENE_SEX_CD","BENE_RACE_CD", "MA",
                           "DUAL", "N_HTN",
                           "PERSON_YEAR_HTN",
                           "PM25_US",
                           "REGION",
                           covariates
)]
# Clean the data
data <- data[data$PERSON_YEAR_HTN > 0,]
data <- as.data.frame(data)
filtered_zips <- data %>%
  group_by(ZIP5) %>%
  summarize(max_pm2.5 = max(PM25_US, na.rm = TRUE)) %>%  
  filter(max_pm2.5 < 9) %>%  
  select(ZIP5)
data <- data %>%
  filter(ZIP5 %in% filtered_zips$ZIP5)

# Construct dataset for double negative control analysis
data_nc <- data %>%
  arrange(ZIP5, YEAR, AGE_CAT, BENE_SEX_CD, BENE_RACE_CD, DUAL, MA) %>%  
  mutate(YEAR_LAG = YEAR - 1) %>%
  left_join(data %>% 
              select(ZIP5, YEAR, AGE_CAT,BENE_SEX_CD, MA,
                     BENE_RACE_CD,DUAL, 
                     N_HTN,
                     PERSON_YEAR_HTN,
                    ) %>% 
              rename(N_HTN_lag = N_HTN, 
                     PERSON_YEAR_HTN_lag = PERSON_YEAR_HTN,
                     YEAR_LAG = YEAR), 
            by = c("ZIP5", "YEAR_LAG", "AGE_CAT","BENE_SEX_CD","BENE_RACE_CD","DUAL", "MA")) %>%  
  select(-YEAR_LAG)

data_nc <- data_nc %>%
  arrange(ZIP5, YEAR, AGE_CAT,BENE_SEX_CD,BENE_RACE_CD,DUAL, MA) %>%  
  mutate(YEAR_LAG = YEAR + 1) %>%
  left_join(data %>% 
              select(ZIP5, YEAR, AGE_CAT,BENE_SEX_CD, MA,
                     BENE_RACE_CD,DUAL, YEAR, PM25_US
              ) %>% 
              rename(PM25_US_lag = PM25_US, 
                     YEAR_LAG = YEAR), 
            by = c("ZIP5", "YEAR_LAG", "AGE_CAT","BENE_SEX_CD","BENE_RACE_CD","DUAL", "MA")) %>%  
  select(-YEAR_LAG)

covariates <- c( "NDVI_SUMMER","NatWalkInd", "RUCA1",
  "PCT_AGE65", "PCT_AGE65_POVERTY", "PCT_AGE65_LOWEDU", 
  "PCT_MINORITY",  "HOUSE_NO_INTERNET", "HOUSE_CROWD", "HOUSE_COST", 
  "UNEMPLOY", "POP_DENS",  "DRINK",  "SMOKE",  "OBESITY", 
  "PHYSINACT", "ANNUAL_CHECK", "CHOLSCREEN", "TMEAN_WARM", "RH_WARM", "TMEAN_COLD", "RH_COLD")
data_select <- data_nc[, c("ZIP5", "YEAR", "AGE_CAT","BENE_SEX_CD","BENE_RACE_CD", "MA",
                           "DUAL", "N_HTN",
                           "PERSON_YEAR_HTN", "N_HTN_lag","PERSON_YEAR_HTN_lag",
                           "PM25_US",
                           "PM25_US_lag",
                           "REGION",
                           covariates
)]
data_cc_1 <- data_select[complete.cases(data_select),]

# Analysis
outcomes <- c("HTN") 
results <- data.frame()

for (outcome in outcomes) {
  print(outcome)
  print(detectCores())
  cl = makeCluster(16,outfile='')
  registerDoParallel(cl)
  
  data_cc <- data_cc_1
  
  
  formula_1 <- as.formula(paste0("N_", outcome, "_lag ~ PM25_US + PM25_US_lag + RUCA1  +
                     NDVI_SUMMER + NatWalkInd + 
                     PCT_AGE65 + PCT_AGE65_POVERTY + 
                     PCT_AGE65_LOWEDU + PCT_MINORITY + HOUSE_NO_INTERNET + HOUSE_CROWD + HOUSE_COST + UNEMPLOY + POP_DENS + DRINK + SMOKE + OBESITY + PHYSINACT +
                     ANNUAL_CHECK + CHOLSCREEN + TMEAN_WARM + TMEAN_COLD + RH_WARM + RH_COLD + as.factor(YEAR) + as.factor(REGION) + 
                     as.factor(BENE_SEX_CD) + as.factor(BENE_RACE_CD) + as.factor(DUAL) + as.factor(AGE_CAT)  + as.factor(MA) + 
                                 offset(log(PERSON_YEAR_", outcome, "_lag))"))
  model_1 <- bam(formula_1,data = data_cc,family = quasipoisson(link="log"),chunk.size = 5000,
                 cluster = cl,control = gam.control(trace=TRUE))
  data_cc$expected_n_lag1 <- fitted(model_1)
  
  
  formula_2 <- as.formula(paste0("N_", outcome, " ~ expected_n_lag1 + PM25_US + RUCA1 +
                     NDVI_SUMMER + NatWalkInd + 
                     PCT_AGE65 + PCT_AGE65_POVERTY + 
                     PCT_AGE65_LOWEDU + PCT_MINORITY + HOUSE_NO_INTERNET + HOUSE_CROWD + HOUSE_COST + UNEMPLOY + POP_DENS + DRINK + SMOKE + OBESITY + PHYSINACT +
                     ANNUAL_CHECK + CHOLSCREEN + TMEAN_WARM + TMEAN_COLD + RH_WARM + RH_COLD + as.factor(YEAR) + as.factor(REGION) + 
                      as.factor(BENE_SEX_CD) + as.factor(BENE_RACE_CD) + as.factor(DUAL) + as.factor(AGE_CAT)  + as.factor(MA) + 
                                 offset(log(PERSON_YEAR_", outcome, "))"))
  model_2 <- bam(formula_2,data = data_cc,family = quasipoisson(link="log"),chunk.size = 5000,
                 cluster = cl,control = gam.control(trace=TRUE))
  
  beta_est <- summary(model_2)$p.coeff["PM25_US"]
  se_est <- summary(model_2)$se["PM25_US"]
  

  lower_CI <- exp(beta_est - 1.96 * se_est)
  upper_CI <- exp(beta_est + 1.96 * se_est)
  
  
  percent_change_lower <- (lower_CI - 1) * 100
  percent_change_upper <- (upper_CI - 1) * 100
  
  results <- rbind(results, data.frame(
    Outcome = outcome,
    Estimate = exp(beta_est), 
    Lower_CI = lower_CI, 
    Upper_CI = upper_CI, 
    Percent_Change = (exp(beta_est) - 1) * 100,
    Percent_Change_Lower = percent_change_lower,
    Percent_Change_Upper = percent_change_upper
  ))
}

results$Estimate <- round(results$Estimate, 2)
results$Lower_CI <- round(results$Lower_CI, 2)
results$Upper_CI <- round(results$Upper_CI, 2)
results$Percent_Change <- round(results$Percent_Change, 2)
results$Percent_Change_Lower <- round(results$Percent_Change_Lower, 2)
results$Percent_Change_Upper <- round(results$Percent_Change_Upper, 2)

results

