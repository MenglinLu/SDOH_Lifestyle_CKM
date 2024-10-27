library(ppcor)

data <- read.csv('E:/lifestyle_ckm/data_UKB.csv')

result_df = data.frame()
for(var1 in c('Lifestyle_score'
)){
  for(var2 in c('eGFR','Creatinine', 'CysC', 'Phosphate',
                'Urate', 'Urea','Albumin',
                'BMI','Waist','Cholesterol','HDLc','LDL','TG',
                'Glucose','HbA1c',
                'SBP','DBP','BP_minus','MAP',
                'CRP','NLR')){
    correlation_spearman <- pcor.test(data[,var1], data[,var2], data[,c('Age','Sex','education_raw','Townsend','ckd_med','glucose_medication','chol_lower_medication','blood_pressure_medication')],method = "spearman")
    result_df <- rbind(result_df, c(var1, var2, correlation_spearman$estimate,correlation_spearman$p.value))
    print(var2)
  }
}

colnames(result_df) <- c('Lifestyle','biomarker','Corr','pval')
result_df$P_FDR <- p.adjust(result_df$pval, method = "BH")