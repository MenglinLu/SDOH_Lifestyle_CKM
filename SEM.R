library(lavaan)
library(bruceR)

data <- read.csv('E:/lifestyle_ckm/data_UKB.csv')

for(fea_i in c(
  'eGFR',
  'HbA1c',
  'BMI',
  'NLR',
  'TG',
  'MAP',
  'Age','Townsend','Lifestyle_score',
  'prs','psrs')){
  z_scores <- scale(data[fea_i])
  outliers <- abs(z_scores)> 5 
  data[fea_i][outliers] <- NA
}

data$TG <- log(data$TG - min(data$TG) + 1)

for(fea_i in c('eGFR',
               'HbA1c',
               'BMI',
               'NLR',
               'TG',
               'MAP',
               'Age','Townsend','Lifestyle_score',
               'prs','psrs')){
  data[fea_i] <- scale(data[fea_i])
}

sem_model <- '
  RenalFunction =~ 1*eGFR
  ImmunoMetabolic =~ 1*BMI #TG HbA1c MAP NLR
  Lifestyle_ff =~ 1*Lifestyle_score
  SDOH =~ 1*psrs
  
  Lifestyle_ff ~ SDOH + Age + Sex + education_raw+Townsend+ckd_med+glucose_medication+chol_lower_medication+blood_pressure_medication
  ImmunoMetabolic ~SDOH + Lifestyle_ff + Age + Sex + education_raw+Townsend+ckd_med+glucose_medication+chol_lower_medication+blood_pressure_medication
  RenalFunction ~ SDOH + Lifestyle_ff + ImmunoMetabolic +Age + Sex + education_raw+Townsend+ckd_med+glucose_medication+chol_lower_medication+blood_pressure_medication
  
  label_CVD ~ SDOH+Lifestyle_ff+prs + ImmunoMetabolic+RenalFunction+Age +Sex+education_raw+Townsend+ckd_med+glucose_medication+chol_lower_medication+blood_pressure_medication
 '

fit <- sem(sem_model, data = data)

fitMeasures(fit,c("chisq","df","pvalue",
                  "gfi","cfi","rmr","srmr","rmsea"))

summary_fit <- lavaan_summary(fit)

p_values <- summary_fit$regression$pval

adjusted_p_values <- p.adjust(p_values, method = "BH")

summary_fit$regression$adjusted_pvalue <- adjusted_p_values