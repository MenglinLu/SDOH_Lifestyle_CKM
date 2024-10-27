library("survival")
library("survminer")
library(epiDisplay) 
library(ggsurvfit)
library(survival)
library(survminer)
library(KMunicate)
library(ComparisonSurv)
library(stringr)

##生活方式分项
data <- read.csv('E:/lifestyle_ckm/data_UKB.csv')

data[data['Lifestyle_score']<3,'Lifestyle_class'] = 'Unfavorable'
data[data['Lifestyle_score']>=3&data['Lifestyle_score']<=4,'Lifestyle_class'] = 'Intermediate'
data[data['Lifestyle_score']>4,'Lifestyle_class'] = 'Favorable'

data[data['psrs']<6,'PSRS_class'] = 'Low'
data[data['psrs']>=6&data['psrs']<10,'PSRS_class'] = 'Intermediate'
data[data['psrs']>=10,'PSRS_class'] = 'High'

data$interactiongroup<-NA
data$interactiongroup[which(data$PSRS_class=='High'&data$Lifestyle_class=="Unfavorable")]<-7
data$interactiongroup[which(data$PSRS_class=='High'&data$Lifestyle_class=="Intermediate")]<-8
data$interactiongroup[which(data$PSRS_class=='High'&data$Lifestyle_class=="Favorable")]<-9
data$interactiongroup[which(data$PSRS_class=="Intermediate"&data$Lifestyle_class=="Unfavorable")]<-4
data$interactiongroup[which(data$PSRS_class=="Intermediate"&data$Lifestyle_class=="Intermediate")]<-5
data$interactiongroup[which(data$PSRS_class=="Intermediate"&data$Lifestyle_class=="Favorable")]<-6
data$interactiongroup[which(data$PSRS_class=='Low'&data$Lifestyle_class=="Unfavorable")]<-1
data$interactiongroup[which(data$PSRS_class=='Low'&data$Lifestyle_class=="Intermediate")]<-2
data$interactiongroup[which(data$PSRS_class=='Low'&data$Lifestyle_class=="Favorable")]<-3

task_li <- c('CVD') #'HF','AF','CHD','PAD','CH','CI'
result_df <- data.frame()
###lifestyle class
for(task_i in task_li){
  data_task <- data
  
  data_task_reference <- data_task[data_task$interactiongroup==1,]
  data_task_reference$subgroup <- 0
  for(group_i in c(1,2,3,4,5,6,7,8,9)){
    data_task_exp <- data_task[data_task$interactiongroup==group_i,]
    data_task_exp$subgroup <- 1
    data_task_reference_exp <- rbind(data_task_reference,data_task_exp)
    cox_res <- coxph(Surv(get(str_c('survival_day_',task_i)), get(str_c('label_',task_i)))~subgroup+Age+Sex+education_raw+Townsend+ckd_med+glucose_medication+chol_lower_medication+blood_pressure_medication,
                     data=data_task_reference_exp)
    
    cox_summary<- summary(cox_res)
    
    HR_stage <- cox_summary$conf.int['subgroup','exp(coef)']
    HR_stage_upper <- cox_summary$conf.int['subgroup','upper .95']
    HR_stage_lower <- cox_summary$conf.int['subgroup','lower .95']
    HR_p <- cox_summary$coefficients['subgroup','Pr(>|z|)']
    
    result_df <- rbind(result_df, c(task_i, group_i, sum(data_task_reference_exp$subgroup == 1 & data_task_reference_exp[,str_c('label_',task_i)] == 1), sum(data_task_reference_exp$subgroup),HR_stage,HR_stage_lower,HR_stage_upper,HR_p))
  }
}

colnames(result_df) <- c('Disease','interactiongroup','Ncases','Total','HR','HR_lower','HR_upper','pval')
result_df$P_FDR <- p.adjust(result_df$pval, method = "BH")


data1 <- data[data$interactiongroup %in% c(1,3,7,9),]
data1[data1$PSRS_class=='High','PSRS_class'] = 1
data1[data1$PSRS_class=='Low','PSRS_class'] = 0
data1[data1$Lifestyle_class=='Unfavorable','Lifestyle_class'] = 0
data1[data1$Lifestyle_class=='Favorable','Lifestyle_class'] = 1

library(interactionR)

fit1 <- coxph(Surv(survival_day_cvd,label_cvd)~PSRS_class*Lifestyle_class+Age+Sex+education_raw+Townsend+ckd_med+glucose_medication+chol_lower_medication+blood_pressure_medication,data=data1)
summary(fit1)

out<- interactionR(fit1,
                   exposure_names =c("PSRS_class", "Lifestyle_class"),
                   ci.type ="delta", ci.level = 0.95,em = F, recode = F)

out$dframe
