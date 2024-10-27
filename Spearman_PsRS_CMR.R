library(ppcor)

data <- read.csv('E:/lifestyle_ckm/data_UKB.csv')

result_df = data.frame()
for(var1 in c('psrs'
)){
  for(var2 in c('AAo_dis',
                'AAo_max', 'AAo_min', 'DAo_dis', 'DAo_max', 'DAo_min', 
                'LAEF',
                'LAV_max', 'LAV_min', 'LASV', 
                'LVCO', 'LVEF', 'LVEDV', 'LVESV', 'LVM','LVSV', 
                'RAEF', 'RAV_max', 'RAV_min', 'RASV', 
                'RVEF', 'RVEDV', 'RVESV','RVSV', 
                'WT_global','Ecc_global','Err_global', 'Ell_global', 
                'WT_AHA1', 'WT_AHA2',
                'WT_AHA3', 'WT_AHA4', 'WT_AHA5', 'WT_AHA6', 'WT_AHA7', 'WT_AHA8',
                'WT_AHA9', 'WT_AHA10', 'WT_AHA11', 'WT_AHA12', 'WT_AHA13', 'WT_AHA14',
                'WT_AHA15', 'WT_AHA16', 
                'Ecc_AHA1', 'Ecc_AHA2', 'Ecc_AHA3', 'Ecc_AHA4',
                'Ecc_AHA5', 'Ecc_AHA6', 'Ecc_AHA7', 'Ecc_AHA8', 'Ecc_AHA9', 'Ecc_AHA10',
                'Ecc_AHA11', 'Ecc_AHA12', 'Ecc_AHA13', 'Ecc_AHA14', 'Ecc_AHA15',
                'Ecc_AHA16',  
                'Err_AHA1', 'Err_AHA2', 'Err_AHA3',
                'Err_AHA4', 'Err_AHA5', 'Err_AHA6', 'Err_AHA7', 'Err_AHA8', 'Err_AHA9',
                'Err_AHA10', 'Err_AHA11', 'Err_AHA12', 'Err_AHA13', 'Err_AHA14',
                'Err_AHA15', 'Err_AHA16', 
                'Ell_seg1', 'Ell_seg2',
                'Ell_seg3', 'Ell_seg4', 'Ell_seg5', 'Ell_seg6')){
    data1 <- data[(!is.na(data[var1])&!is.na(data[var2])),]
    correlation_spearman <- pcor.test(data1[,var1], data1[,var2], data1[,c('Age','Sex','education_raw','Townsend', 'ckd_med','glucose_medication','chol_lower_medication','blood_pressure_medication')],method = "spearman")
    result_df <- rbind(result_df, c(var1, var2, correlation_spearman$estimate,correlation_spearman$p.value))
  }
}

colnames(result_df) <- c('PSRS','CMR','Corr','pval')
result_df$P_FDR <- p.adjust(result_df$pval, method = "BH")

