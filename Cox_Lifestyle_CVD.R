library("survival")
library("survminer")
library(epiDisplay) 
library(ggsurvfit)
library(survival)
library(survminer)
library(KMunicate)
library(ComparisonSurv)
library(stringr)

###lifestyle factor
data <- read.csv('E:/lifestyle_ckm/data_UKB.csv')

data[data['Lifestyle_score']<3,'Lifestyle_class'] = 'Unfavorable'
data[data['Lifestyle_score']>=3&data['Lifestyle_score']<=4,'Lifestyle_class'] = 'Intermediate'
data[data['Lifestyle_score']>4,'Lifestyle_class'] = 'Favorable'

task_li <- c('CVD') #'HF','AF','CHD','PAD','CH','CI'

lifestyle_li <- c('Smoking_healthy','Sedentary_healthy',
                  'Sleep_healthy','Activity_healthy','Alcohol_healthy','Diet_healthy')

result_df <- data.frame()
for(task_i in task_li){
  data_task = data
  for(lifestyle_i in lifestyle_li){
    data_task_reference_exp <- data_task[!is.na(data_task[,lifestyle_i]),]
    data_task_reference_exp$subgroup <- as.factor(data_task_reference_exp[,lifestyle_i])
    
    fit <- coxph(Surv(get(str_c('survival_day_',task_i)), get(str_c('label_',task_i)))~subgroup+Age+Sex+education_raw+Townsend+ckd_med+glucose_medication+chol_lower_medication+blood_pressure_medication,
                 data=data_task_reference_exp)
    
    cox_summary <- summary(fit)
    
    HR_stage <- cox_summary$conf.int['subgroup1','exp(coef)']
    HR_stage_upper <- cox_summary$conf.int['subgroup1','upper .95']
    HR_stage_lower <- cox_summary$conf.int['subgroup1','lower .95']
    HR_p <- cox_summary$coefficients['subgroup1','Pr(>|z|)']
    
    result_df <- rbind(result_df, c(task_i, lifestyle_i, sum(data_task_reference_exp[lifestyle_i] == 1 & data_task_reference_exp[str_c('label_',task_i)] == 1), sum(data_task_reference_exp[lifestyle_i]),HR_stage,HR_stage_lower,HR_stage_upper,HR_p))
  }
}

###lifestyle class
for(task_i in task_li){
  data_task <- data
  
  data_task_reference_exp <- data_task[!is.na(data_task['Lifestyle_class']),]
  data_task_reference_exp <- data_task_reference_exp[data_task_reference_exp$Lifestyle_class %in% c('Unfavorable','Intermediate','Favorable'),]
  
  data_task_reference_exp$Lifestyle_class <- ifelse(data_task_reference_exp$Lifestyle_class == "Unfavorable", 0,
                                                    ifelse(data_task_reference_exp$Lifestyle_class == "Intermediate", 1,
                                                           ifelse(data_task_reference_exp$Lifestyle_class == "Favorable", 2, NA)))
  
  data_task_reference_exp$subgroup <- as.factor(data_task_reference_exp[,'Lifestyle_class'])
  
  data_task_reference_exp01 <- data_task_reference_exp[data_task_reference_exp$subgroup == 0, ]
  data_task_reference_exp00 <- data_task_reference_exp01
  data_task_reference_exp00$subgroup = 1
  data_task_reference_exp0 <- rbind(data_task_reference_exp01,data_task_reference_exp00)
  data_task_reference_exp1 <- data_task_reference_exp[data_task_reference_exp$subgroup == 0 | data_task_reference_exp$subgroup == 1, ]
  data_task_reference_exp2 <- data_task_reference_exp[data_task_reference_exp$subgroup == 0 | data_task_reference_exp$subgroup == 2, ]
  
  cox_res <- coxph(Surv(get(str_c('survival_day_',task_i)), get(str_c('label_',task_i)))~subgroup+Age+Sex+education_raw+Townsend,
                   data=data_task_reference_exp0)
  
  cox_summary<- summary(cox_res)
  
  HR_stage <- cox_summary$conf.int['subgroup1','exp(coef)']
  HR_stage_upper <- cox_summary$conf.int['subgroup1','upper .95']
  HR_stage_lower <- cox_summary$conf.int['subgroup1','lower .95']
  HR_p <- cox_summary$coefficients['subgroup1','Pr(>|z|)']
  
  result_df <- rbind(result_df, c(task_i, 'Unfavorable', sum(data_task_reference_exp00$Lifestyle_class == 0 & data_task_reference_exp00[str_c('label_',task_i)] == 1), sum(data_task_reference_exp00$Lifestyle_class==0),HR_stage,HR_stage_lower,HR_stage_upper,HR_p))
  
  cox_res <- coxph(Surv(get(str_c('survival_day_',task_i)), get(str_c('label_',task_i)))~subgroup+Age+Sex+education_raw+Townsend,
                   data=data_task_reference_exp1)
  
  cox_summary<- summary(cox_res)
  
  HR_stage <- cox_summary$conf.int['subgroup1','exp(coef)']
  HR_stage_upper <- cox_summary$conf.int['subgroup1','upper .95']
  HR_stage_lower <- cox_summary$conf.int['subgroup1','lower .95']
  HR_p <- cox_summary$coefficients['subgroup1','Pr(>|z|)']
  
  result_df <- rbind(result_df, c(task_i, 'Intermediate', sum(data_task_reference_exp1$Lifestyle_class == 1 & data_task_reference_exp1[str_c('label_',task_i)] == 1), sum(data_task_reference_exp1$Lifestyle_class==1),HR_stage,HR_stage_lower,HR_stage_upper,HR_p))
  
  cox_res <- coxph(Surv(get(str_c('survival_day_',task_i)), get(str_c('label_',task_i)))~subgroup+Age+Sex+education_raw+Townsend,
                   data=data_task_reference_exp2)
  
  cox_summary<- summary(cox_res)
  
  HR_stage <- cox_summary$conf.int['subgroup2','exp(coef)']
  HR_stage_upper <- cox_summary$conf.int['subgroup2','upper .95']
  HR_stage_lower <- cox_summary$conf.int['subgroup2','lower .95']
  HR_p <- cox_summary$coefficients['subgroup2','Pr(>|z|)']
  
  result_df <- rbind(result_df, c(task_i, 'Favorable', sum(data_task_reference_exp2$Lifestyle_class == 2 & data_task_reference_exp2[str_c('label_',task_i)] == 1), sum(data_task_reference_exp2$Lifestyle_class==2),HR_stage,HR_stage_lower,HR_stage_upper,HR_p))
}

colnames(result_df) <- c('Disease','Lifestyle','Ncases','Total','HR','HR_lower','HR_upper','pval')
result_df$P_FDR <- p.adjust(result_df$pval, method = "BH")
