library(data.table)
library(RACER)
library(TwoSampleMR)
library(ggplot2)
effective_n <- function(ncase, ncontrol){ return(round(2 / (1/ncase + 1/ncontrol),0))} #taken from TWoSampleMR

######################################################################################

## MR neg cont exposures
##########################
#Finn r12, UKB, MVP HF
hf<-fread(paste("\\GWAS\\meta_analysis_mvp_ukbb_summary_stats_I9_HEARTFAIL_meta_out.tsv.gz",sep=""))
hf<-hf[hf$all_inv_var_meta_p < 5*10^-8,]; hf$pval<-hf$all_inv_var_meta_p
hf<-ieugwasr::ld_clump(  hf,  clump_kb = 10000,  clump_r2 = 0.001,  clump_p = 5*10^-8,  plink_bin = paste("C:/Users/",Benjamin,"/OneDrive - University of Bristol/Documents/Projects/active/clumping/plink_win64_20230116/plink.exe",sep=""),  bfile = paste("C:/Users/",Benjamin,"/OneDrive - University of Bristol/Documents/Projects/active/clumping/EUR",sep="") )
hf$pheno<-"HF"
hf<-format_data(as.data.frame(hf), type = "exposure",  phenotype_col = "pheno", snp_col = "rsid",  beta_col = "all_inv_var_meta_beta",  se_col = "all_inv_var_meta_sebeta",  eaf_col = "fg_af_alt",  effect_allele_col = "ALT",  other_allele_col = "REF",  pval_col = "all_inv_var_meta_p",  id_col = "pheno",  chr_col = "#CHR",  pos_col = "POS")
dat<-hf;hf<-NULL

#Finn r12, UKB, MVP IPF
ipf<-fread(paste("\\GWAS\\meta_analysis_mvp_ukbb_summary_stats_IPF_meta_out.tsv.gz",sep=""))
ipf<-ipf[ipf$all_inv_var_meta_p < 5*10^-8,]; ipf$pval<-ipf$all_inv_var_meta_p
ipf<-ieugwasr::ld_clump(  ipf,  clump_kb = 10000,  clump_r2 = 0.001,  clump_p = 5*10^-8,  plink_bin = paste("C:/Users/",Benjamin,"/OneDrive - University of Bristol/Documents/Projects/active/clumping/plink_win64_20230116/plink.exe",sep=""),  bfile = paste("C:/Users/",Benjamin,"/OneDrive - University of Bristol/Documents/Projects/active/clumping/EUR",sep="") )
ipf$pheno<-"IPF"
ipf<-format_data(as.data.frame(ipf), type = "exposure",  phenotype_col = "pheno", snp_col = "rsid",  beta_col = "all_inv_var_meta_beta",  se_col = "all_inv_var_meta_sebeta",  eaf_col = "fg_af_alt",  effect_allele_col = "ALT",  other_allele_col = "REF",  pval_col = "all_inv_var_meta_p",  id_col = "pheno",  chr_col = "#CHR",  pos_col = "POS")
dat<-rbind(dat,ipf);ipf<-NULL

#Finn r12, UKB, MVP PE
pe<-fread(paste("\\GWAS\\meta_analysis_mvp_ukbb_summary_stats_I9_PULMEMB_meta_out.tsv.gz",sep=""))
pe<-pe[pe$all_inv_var_meta_p < 5*10^-8,]; pe$pval<-pe$all_inv_var_meta_p
pe<-ieugwasr::ld_clump(  pe,  clump_kb = 10000,  clump_r2 = 0.001,  clump_p = 5*10^-8,  plink_bin = paste("C:/Users/",Benjamin,"/OneDrive - University of Bristol/Documents/Projects/active/clumping/plink_win64_20230116/plink.exe",sep=""),  bfile = paste("C:/Users/",Benjamin,"/OneDrive - University of Bristol/Documents/Projects/active/clumping/EUR",sep="") )
pe$pheno<-"PE"
pe<-format_data(as.data.frame(pe), type = "exposure",  phenotype_col = "pheno", snp_col = "rsid",  beta_col = "all_inv_var_meta_beta",  se_col = "all_inv_var_meta_sebeta",  eaf_col = "fg_af_alt",  effect_allele_col = "ALT",  other_allele_col = "REF",  pval_col = "all_inv_var_meta_p",  id_col = "pheno",  chr_col = "#CHR",  pos_col = "POS")
dat<-rbind(dat,pe);pe<-NULL

#Global biobank Meta COPD
copd<-fread(paste("\\GWAS\\COPD_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz",sep=""))
copd<-copd[copd$inv_var_meta_p < 5*10^-8,]; copd$pval<-copd$inv_var_meta_p
copd<-ieugwasr::ld_clump(  copd,  clump_kb = 10000,  clump_r2 = 0.001,  clump_p = 5*10^-8,  plink_bin = paste("C:/Users/",Benjamin,"/OneDrive - University of Bristol/Documents/Projects/active/clumping/plink_win64_20230116/plink.exe",sep=""),  bfile = paste("C:/Users/",Benjamin,"/OneDrive - University of Bristol/Documents/Projects/active/clumping/EUR",sep="") )
copd$pheno<-"COPD"
copd<-format_data(as.data.frame(copd), type = "exposure",  phenotype_col = "pheno", snp_col = "rsid",  beta_col = "inv_var_meta_beta",  se_col = "inv_var_meta_sebeta",  eaf_col = "all_meta_AF",  effect_allele_col = "ALT",  other_allele_col = "REF",  pval_col = "inv_var_meta_p",  id_col = "pheno",  chr_col = "#CHR",  pos_col = "POS")
dat<-rbind(dat,copd);copd<-NULL

#Finn R12 with UKB  Asthma
asth<-fread(paste("\\GWAS\\meta_analysis_ukbb_summary_stats_finngen_R12_J10_ASTHMA_MAIN_EXMORE_meta_out.tsv.gz",sep=""))
asth<-asth[asth$all_inv_var_meta_p < 5*10^-8,]; asth$pval<-asth$all_inv_var_meta_p
asth<-ieugwasr::ld_clump(  asth,  clump_kb = 10000,  clump_r2 = 0.001,  clump_p = 5*10^-8,  plink_bin = paste("C:/Users/",Benjamin,"/OneDrive - University of Bristol/Documents/Projects/active/clumping/plink_win64_20230116/plink.exe",sep=""),  bfile = paste("C:/Users/",Benjamin,"/OneDrive - University of Bristol/Documents/Projects/active/clumping/EUR",sep="") )
asth$pheno<-"ASTH"
asth<-format_data(as.data.frame(asth), type = "exposure",  phenotype_col = "pheno", snp_col = "rsid",  beta_col = "all_inv_var_meta_beta",  se_col = "all_inv_var_meta_sebeta",  eaf_col = "FINNGEN_af_alt",  effect_allele_col = "ALT",  other_allele_col = "REF",  pval_col = "all_inv_var_meta_p",  id_col = "pheno",  chr_col = "#CHR",  pos_col = "POS")
dat<-rbind(dat,asth);asth<-NULL

write.csv(dat, "PH sumstats\\biobank_negcont_expSNPs_Finnr12ukbMVP.csv")

################################################################################
# loading MVP-UKB-FinnGenR12 PAH data

PH<-fread(paste("\\meta_analysis_mvp_ukbb_summary_stats_I9_HYPTENSPUL_meta_out.tsv.gz",sep=""))

#3,302 cases 1,205,457 controls
#https://mvp-ukbb.finngen.fi/pheno/I9_HYPTENSPUL


#loading rhodes et al 
rhodes<-data.table::fread(paste("PH sumstats\\GCST007228_buildGRCh37_metaPAH.txt.gz",sep=""))


###############################
#checking Rhodes et al hits

PH[PH$rsid=="rs2856830"|
     PH$rsid=="rs13266183"|
     PH$rsid=="rs10103692",]

# #CHR      POS    REF    ALT            SNP fg_beta fg_sebeta fg_pval fg_af_alt fg_af_alt_cases
# <int>    <int> <char> <char>         <char>   <num>     <num>   <num>     <num>           <num>
#   1:     6 33073957      T      C 6:33073957:T:C  0.2890     0.104 0.00575  0.143689        0.186066
# 2:     8 54345567      A      G 8:54345567:A:G -0.0883     0.128 0.49000  0.114452        0.106293
# 3:     8 54355052      C      T 8:54355052:C:T -0.0286     0.103 0.78100  0.202426        0.196541
# fg_af_alt_controls MVP_EUR_beta MVP_EUR_sebeta MVP_EUR_pval MVP_EUR_af_alt MVP_EUR_r2 MVP_AFR_beta
# <num>        <num>          <num>        <num>          <num>      <num>        <num>
#   1:           0.143653      -0.0450         0.0512        0.380         0.1264     0.9998       0.0771
# 2:           0.114459       0.0408         0.0598        0.495         0.0907     0.9739      -0.0421
# 3:           0.202432       0.0178         0.0408        0.663         0.2515     0.9177      -0.2230
# MVP_AFR_sebeta MVP_AFR_pval MVP_AFR_af_alt MVP_AFR_r2 ukbb_beta ukbb_sebeta ukbb_pval ukbb_af_alt
# <num>        <num>          <num>      <num>     <num>       <num>     <num>       <num>
#   1:         0.1020      0.44800         0.0934     0.9996   0.00611      0.0888    0.9450    0.114485
# 2:         0.0765      0.58200         0.1822     0.9606   0.25500      0.0992    0.0102    0.090867
# 3:         0.0806      0.00564         0.1990     0.9033  -0.05160      0.0648    0.4260    0.264572
# all_meta_N all_inv_var_meta_beta all_inv_var_meta_sebeta all_inv_var_meta_p all_inv_var_meta_mlogp
# <int>                 <num>                   <num>              <num>                  <num>
#   1:          4                0.0252                  0.0379              0.506                   0.30
# 2:          4                0.0403                  0.0404              0.318                   0.50
# 3:          4               -0.0355                  0.0303              0.241                   0.62
# all_inv_var_het_p leave_fg_N leave_fg_inv_var_meta_beta leave_fg_inv_var_meta_sebeta
# <num>      <int>                      <num>                        <num>
#   1:            0.0361          3                    -0.0147                       0.0407
# 2:            0.0771          3                     0.0546                       0.0426
# 3:            0.0661          3                    -0.0362                       0.0317
# leave_fg_inv_var_meta_p leave_fg_inv_var_meta_mlogp leave_fg_inv_var_meta_het_p leave_MVP_EUR_N
# <num>                       <num>                       <num>           <int>
#   1:                   0.718                        0.14                      0.5440               3
# 2:                   0.200                        0.70                      0.0573               3
# 3:                   0.254                        0.60                      0.0275               3
# leave_MVP_EUR_inv_var_meta_beta leave_MVP_EUR_inv_var_meta_sebeta leave_MVP_EUR_inv_var_meta_p
# <num>                             <num>                        <num>
#   1:                          0.1100                            0.0563                       0.0508
# 2:                          0.0399                            0.0548                       0.4660
# 3:                         -0.1010                            0.0453                       0.0254
# leave_MVP_EUR_inv_var_meta_mlogp leave_MVP_EUR_inv_var_meta_het_p leave_MVP_AFR_N
# <num>                            <num>           <int>
#   1:                             1.29                           0.1110               3
# 2:                             0.33                           0.0327               3
# 3:                             1.60                           0.1850               3
# leave_MVP_AFR_inv_var_meta_beta leave_MVP_AFR_inv_var_meta_sebeta leave_MVP_AFR_inv_var_meta_p
# <num>                             <num>                        <num>
#   1:                         0.01680                            0.0408                        0.680
# 2:                         0.07210                            0.0476                        0.129
# 3:                        -0.00464                            0.0327                        0.887
# leave_MVP_AFR_inv_var_meta_mlogp leave_MVP_AFR_inv_var_meta_het_p leave_ukbb_N
# <num>                            <num>        <int>
#   1:                             0.17                           0.0163            3
# 2:                             0.89                           0.0730            3
# 3:                             0.05                           0.6440            3
# leave_ukbb_inv_var_meta_beta leave_ukbb_inv_var_meta_sebeta leave_ukbb_inv_var_meta_p
# <num>                          <num>                     <num>
#   1:                      0.02950                         0.0419                     0.482
# 2:                     -0.00228                         0.0442                     0.959
# 3:                     -0.03100                         0.0343                     0.366
# leave_ukbb_inv_var_meta_mlogp leave_ukbb_inv_var_meta_het_p       rsid
# <num>                         <num>     <char>
#   1:                          0.32                        0.0144  rs2856830
# 2:                          0.02                        0.5370 rs10103692
# 3:                          0.44                        0.0286 rs13266183


PH[PH$rsid=="rs2856830"|
     PH$rsid=="rs13266183"|
     PH$rsid=="rs10103692",c("rsid","all_inv_var_meta_p")]
# rsid all_inv_var_meta_p
# <char>              <num>
#   1:  rs2856830              0.506
# 2: rs10103692              0.318
# 3: rs13266183              0.241


###############################
#looking at FINN-UKB-MVP meta hits
rhodes[rhodes$rsid =="rs142623583"| #p = 2.8e-9 in biobanks
         rhodes$rsid =="rs10766879",c("rsid","p.value")] #p = 5.3e-8 in biobanks 

# rsid  p.value
# <char>    <num>
#   1: rs10766879 0.413992

# # rhodes[rhodes$`chr:pos(hg19)`=="17:79549007",] # gettign chr pos from https://www.ncbi.nlm.nih.gov/snp/?term=rs142623583
# SNP<-rhodes[rhodes$chromosome ==17 & rhodes$base_pair_location>79549007-10000 & rhodes$base_pair_location<79549007+10000,]
# ieugwasr::ld_matrix(c(SNP,"rs142623583"),  opengwas_jwt = jwt, plink_bin = paste("C:/Users/",Benjamin,"/OneDrive - University of Bristol/Documents/Projects/active/clumping/plink_win64_20230116/plink.exe",sep=""),  bfile = paste("C:/Users/",Benjamin,"/OneDrive - University of Bristol/Documents/Projects/active/clumping/EUR",sep=""))["rs142623583_T_A",]
# # no SNPs with r2>.25 within 10Mb so there probs are not any good LD proxies for rs142623583. 


##############################
# looking at NCE

PH2<-PH[PH$rsid%in%dat$SNP]
PH2<-format_data(as.data.frame(PH2), type = "outcome",  phenotype_col = "pheno", snp_col = "rsid",  beta_col = "all_inv_var_meta_beta",  se_col = "all_inv_var_meta_sebeta",  eaf_col = "fg_af_alt",  effect_allele_col = "ALT",  other_allele_col = "REF",  pval_col = "all_inv_var_meta_p",  id_col = "pheno",  chr_col = "#CHR",  pos_col = "POS")
mr(harmonise_data(dat,PH2), method_list="mr_ivw")

# id.exposure id.outcome outcome exposure                    method nsnp          b         se
# 1      167dfH     f207vs outcome      IPF Inverse variance weighted   21 0.03696198 0.03700030
# 2      9xxVGF     f207vs outcome       HF Inverse variance weighted  144 0.64968462 0.07300127
# 3      gFdQjZ     f207vs outcome     COPD Inverse variance weighted   19 0.31930752 0.11769121
# 4      MCZPjG     f207vs outcome     ASTH Inverse variance weighted   68 0.07084503 0.06771678
# 5      sru1l6     f207vs outcome       PE Inverse variance weighted   66 0.10357407 0.04469470
# pval
# 1 3.178119e-01
# 2 5.603090e-19
# 3 6.665783e-03
# 4 2.954705e-01
# 5 2.048370e-02

########################
#neg cont with rhodes

rhodes2<-rhodes[rhodes$rsid %in% dat$SNP,]
rhodes2$pheno<-"PH"
rhodes2<-format_data(as.data.frame(rhodes2), type = "outcome",  phenotype_col = "pheno", snp_col = "rsid",  beta_col = "logOR",  se_col = "OR_se",  eaf_col = "eaf",  effect_allele_col = "effect_allele",  other_allele_col = "other_allele",  pval_col = "p.value",  id_col = "pheno",  chr_col = "chromosome",  pos_col = "base_pair_location")
mr(harmonise_data(dat,rhodes2), method_list="mr_ivw")

# id.exposure id.outcome outcome exposure                    method nsnp            b         se
# 1      167dfH     LeE9mH      PH      IPF Inverse variance weighted   19  0.010931889 0.04854561
# 2      9xxVGF     LeE9mH      PH       HF Inverse variance weighted  140  0.190499736 0.11444769
# 3      gFdQjZ     LeE9mH      PH     COPD Inverse variance weighted   20  0.134712883 0.20136356
# 4      MCZPjG     LeE9mH      PH     ASTH Inverse variance weighted   64  0.168861980 0.10525801
# 5      sru1l6     LeE9mH      PH       PE Inverse variance weighted   63 -0.003072478 0.06692140
# pval
# 1 0.82183302
# 2 0.09600983
# 3 0.50349337
# 4 0.10865515
# 5 0.96338060

##########################################
#doing the analsyis 


AF<-fread(paste("\\GWAS\\meta_analysis_mvp_ukbb_summary_stats_I9_AF_meta_out.tsv.gz",sep=""))
AF<-AF[AF$all_inv_var_meta_p < 5*10^-8,]; AF$pval<-AF$all_inv_var_meta_p
AF<-ieugwasr::ld_clump(  AF,  clump_kb = 10000,  clump_r2 = 0.001,  clump_p = 5*10^-8,  plink_bin = paste("C:/Users/",Benjamin,"/OneDrive - University of Bristol/Documents/Projects/active/clumping/plink_win64_20230116/plink.exe",sep=""),  bfile = paste("C:/Users/",Benjamin,"/OneDrive - University of Bristol/Documents/Projects/active/clumping/EUR",sep="") )
AF$pheno<-"AF"
exp<-format_data(as.data.frame(AF), type = "exposure",  phenotype_col = "pheno", snp_col = "rsid",  beta_col = "all_inv_var_meta_beta",  se_col = "all_inv_var_meta_sebeta",  eaf_col = "fg_af_alt",  effect_allele_col = "ALT",  other_allele_col = "REF",  pval_col = "all_inv_var_meta_p",  id_col = "pheno",  chr_col = "#CHR",  pos_col = "POS")
AF<-NULL

#PH<-fread(paste("\\PH MR\\gwas\\biobank phenotyping\\research letter\\meta_analysis_mvp_ukbb_summary_stats_I9_HYPTENSPUL_meta_out.tsv.gz",sep=""))
PH2<-PH[PH$rsid%in%exp$SNP,]
PH2<-format_data(as.data.frame(PH2), type = "outcome",  phenotype_col = "pheno", snp_col = "rsid",  beta_col = "all_inv_var_meta_beta",  se_col = "all_inv_var_meta_sebeta",  eaf_col = "fg_af_alt",  effect_allele_col = "ALT",  other_allele_col = "REF",  pval_col = "all_inv_var_meta_p",  id_col = "pheno",  chr_col = "#CHR",  pos_col = "POS")
mr(harmonise_data(exp,PH2), method_list="mr_ivw")

##AF
# id.exposure id.outcome outcome exposure                    method nsnp         b         se
# 1      NKBtFB     yBMd6m outcome       AF Inverse variance weighted  349 0.2715708 0.03755704
# pval
# 1 4.798426e-13

#rhodes<-data.table::fread(paste("PH sumstats\\GCST007228_buildGRCh37_metaPAH.txt.gz",sep=""))
rhodes2<-rhodes[rhodes$rsid %in% exp$SNP,]
rhodes2<-format_data(as.data.frame(rhodes2), type = "outcome",  phenotype_col = "pheno", snp_col = "rsid",  beta_col = "logOR",  se_col = "OR_se",  eaf_col = "eaf",  effect_allele_col = "effect_allele",  other_allele_col = "other_allele",  pval_col = "p.value",  id_col = "pheno",  chr_col = "chromosome",  pos_col = "base_pair_location")
mr(harmonise_data(exp,rhodes2), method_list="mr_ivw")

## AF
# id.exposure id.outcome outcome exposure                    method nsnp          b         se
# 1      NKBtFB     U27Imm outcome       AF Inverse variance weighted  338 0.04270227 0.05741546
# pval
# 1 0.4570328


