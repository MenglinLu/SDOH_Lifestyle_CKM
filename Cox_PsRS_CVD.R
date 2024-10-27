library("survival")
library("survminer")
library(epiDisplay) 
library(ggsurvfit)
library(survival)
library(survminer)
library(KMunicate)
library(ComparisonSurv)
library(stringr)

###psrs factor
data <- read.csv('E:/psrs_ckm/data_UKB.csv')

data[data['psrs']<6,'PSRS_class'] = 'Low'
data[data['psrs']>=6&data['psrs']<10,'PSRS_class'] = 'Intermediate'
data[data['psrs']>=10,'PSRS_class'] = 'High'

task_li <- c('CVD') #'HF','AF','CHD','PAD','CH','CI'

psrs_li <- c('psrs_1','psrs_2',
                  'psrs_3','psrs_4','psrs_5','psrs_6',
                  'psrs_7','psrs_8','psrs_9','psrs_10',
                  'psrs_11','psrs_12','psrs_13','psrs_14',
                  'psrs_15','psrs_16','psrs_17')

result_df <- data.frame()
for(task_i in task_li){
  data_task = data
  for(psrs_i in psrs_li){
    data_task_reference_exp <- data_task[!is.na(data_task[,psrs_i]),]
    data_task_reference_exp$subgroup <- as.factor(data_task_reference_exp[,psrs_i])
    
    fit <- coxph(Surv(get(str_c('survival_day_',task_i)), get(str_c('label_',task_i)))~subgroup+Age+Sex+education_raw+Townsend+ckd_med+glucose_medication+chol_lower_medication+blood_pressure_medication,
                 data=data_task_reference_exp)
    
    cox_summary <- summary(fit)
    
    HR_stage <- cox_summary$conf.int['subgroup1','exp(coef)']
    HR_stage_upper <- cox_summary$conf.int['subgroup1','upper .95']
    HR_stage_lower <- cox_summary$conf.int['subgroup1','lower .95']
    HR_p <- cox_summary$coefficients['subgroup1','Pr(>|z|)']
    
    result_df <- rbind(result_df, c(task_i, psrs_i, sum(data_task_reference_exp[psrs_i] == 1 & data_task_reference_exp[str_c('label_',task_i)] == 1), sum(data_task_reference_exp[psrs_i]),HR_stage,HR_stage_lower,HR_stage_upper,HR_p))
  }
}

###psrs class
for(task_i in task_li){
  data_task <- data
  
  data_task_reference_exp <- data_task[!is.na(data_task['PSRS_class']),]
  data_task_reference_exp <- data_task_reference_exp[data_task_reference_exp$PSRS_class %in% c('Low','Intermediate','High'),]
  
  data_task_reference_exp$PSRS_class <- ifelse(data_task_reference_exp$PSRS_class == "Low", 0,
                                               ifelse(data_task_reference_exp$PSRS_class == "Intermediate", 1,
                                                      ifelse(data_task_reference_exp$PSRS_class == "High", 2, NA)))
  
  data_task_reference_exp$subgroup <- as.factor(data_task_reference_exp[,'PSRS_class'])
  
  data_task_reference_exp01 <- data_task_reference_exp[data_task_reference_exp$subgroup == 0, ]
  data_task_reference_exp00 <- data_task_reference_exp01
  data_task_reference_exp00$subgroup = 1
  data_task_reference_exp0 <- rbind(data_task_reference_exp01,data_task_reference_exp00)
  data_task_reference_exp1 <- data_task_reference_exp[data_task_reference_exp$subgroup == 0 | data_task_reference_exp$subgroup == 1, ]
  data_task_reference_exp2 <- data_task_reference_exp[data_task_reference_exp$subgroup == 0 | data_task_reference_exp$subgroup == 2, ]
  
  cox_res <- coxph(Surv(get(str_c('survival_day_',task_i)), get(str_c('label_',task_i)))~subgroup+Age+Sex+education_raw+Townsend+ckd_med+glucose_medication+chol_lower_medication+blood_pressure_medication,
                   data=data_task_reference_exp0)
  
  cox_summary<- summary(cox_res)
  
  HR_stage <- cox_summary$conf.int['subgroup1','exp(coef)']
  HR_stage_upper <- cox_summary$conf.int['subgroup1','upper .95']
  HR_stage_lower <- cox_summary$conf.int['subgroup1','lower .95']
  HR_p <- cox_summary$coefficients['subgroup1','Pr(>|z|)']
  
  result_df <- rbind(result_df, c(task_i, 'Low', sum(data_task_reference_exp00$PSRS_class == 0 & data_task_reference_exp00[str_c('label_',task_i)] == 1), sum(data_task_reference_exp00$PSRS_class==0),HR_stage,HR_stage_lower,HR_stage_upper,HR_p))
  
  cox_res <- coxph(Surv(get(str_c('survival_day_',task_i)), get(str_c('label_',task_i)))~subgroup+Age+Sex+education_raw+Townsend+ckd_med+glucose_medication+chol_lower_medication+blood_pressure_medication,
                   data=data_task_reference_exp1)
  
  cox_summary<- summary(cox_res)
  
  HR_stage <- cox_summary$conf.int['subgroup1','exp(coef)']
  HR_stage_upper <- cox_summary$conf.int['subgroup1','upper .95']
  HR_stage_lower <- cox_summary$conf.int['subgroup1','lower .95']
  HR_p <- cox_summary$coefficients['subgroup1','Pr(>|z|)']
  
  result_df <- rbind(result_df, c(task_i, 'Intermediate', sum(data_task_reference_exp1$PSRS_class == 1 & data_task_reference_exp1[str_c('label_',task_i)] == 1), sum(data_task_reference_exp1$PSRS_class==1),HR_stage,HR_stage_lower,HR_stage_upper,HR_p))
  
  cox_res <- coxph(Surv(get(str_c('survival_day_',task_i)), get(str_c('label_',task_i)))~subgroup+Age+Sex+education_raw+Townsend+ckd_med+glucose_medication+chol_lower_medication+blood_pressure_medication,
                   data=data_task_reference_exp2)
  
  cox_summary<- summary(cox_res)
  
  HR_stage <- cox_summary$conf.int['subgroup2','exp(coef)']
  HR_stage_upper <- cox_summary$conf.int['subgroup2','upper .95']
  HR_stage_lower <- cox_summary$conf.int['subgroup2','lower .95']
  HR_p <- cox_summary$coefficients['subgroup2','Pr(>|z|)']
  
  result_df <- rbind(result_df, c(task_i, 'High', sum(data_task_reference_exp2$PSRS_class == 2 & data_task_reference_exp2[str_c('label_',task_i)] == 1), sum(data_task_reference_exp2$PSRS_class==2),HR_stage,HR_stage_lower,HR_stage_upper,HR_p))
}

colnames(result_df) <- c('Disease','SDOH','Ncases','Total','HR','HR_lower','HR_upper','pval')
result_df$P_FDR <- p.adjust(result_df$pval, method = "BH")
