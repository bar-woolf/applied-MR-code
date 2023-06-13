library(TwoStepCisMR)
#setwd("C:/Users/nf18217/OneDrive - University of Bristol/Documents/Projects/active/cis canabis")
library(data.table)
library(coloc)
library(TwoSampleMR)
library(ggplot2)
#library(gwasglue)
#library(gassocplot)
#library(MRSamePopTest)
#library(TwoStepCisMR)


FIQT_deflation <- function(B,SE, min.p=10^-300){
  z=B/SE
  pvals<-2*pnorm(abs(z),low=F)
  pvals[pvals<min.p]<- min.p
  adj.pvals<-p.adjust(pvals,method="fdr")
  mu.z<-sign(z)*qnorm(adj.pvals/2,low=F)
  mu.z[abs(z)>qnorm(min.p/2,low=F)]<-z[abs(z)>qnorm(min.p/2,low=F)]
  mu.z/z
}



# load GWASs and extract only SNPs that you will use
#####################################################

# exp<-read.csv("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\caff_exp_snps.csv")
# names(exp)[1]<-"SNP"

exp<-data.table::fread("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caff_share.gz")
a<-as.numeric(unlist(stringr::str_split(exp$`chr:pos`,":")))
exp$chr<-a [seq(from=1,to=length(a), by =2)] 
exp$pos<-a [seq(from=2,to=length(a), by =2)]
# CYP1A2  15:75041185-75048543
# AHR 7:17338246-17385776
window=100000
CYP1A2<-exp[exp$chr==15 & exp$pos>75041185-window & exp$pos<75048543+window & exp$P<5e-5,]
CYP1A2$gene<-"CYP1A2"
AHR <-exp[exp$chr==7 & exp$pos>17338246-window & exp$pos<17385776+window & exp$P<5e-5,]
AHR$gene<-"AHR"

exp<-rbind(AHR,CYP1A2)
exp$phenoscanner<-paste("chr",exp$`chr:pos`,"",sep="")
write.table(c( exp$phenoscanner  )
            ,"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\chrpos_em5.txt",row.names = F,col.names = F,quote = F)


id<-data.table::fread("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\chrpos_em5_PhenoScanner_SNP.tsv")
id<-as.data.frame(id)[c("snp","rsid","eur")]
colnames(id)<-c("phenoscanner","SNP","MAF")
exp<-merge(exp,id,by="phenoscanner")
exp<-exp[!(exp$MAF=="-"),]
exp$MAF<-as.numeric(exp$MAF)
exp$beta.outcome<-exp$Effect/(2*exp$MAF*(1-exp$MAF)*(9876 +exp$Effect^2))^0.5
exp$se.outcome<-1/(2*exp$MAF*(1-exp$MAF)*(9876 +exp$Effect^2))^0.5
#gene<-exp[,c("SNP","gene")]
exp <- format_data(exp, type = "exposure", gene_col = "gene", snp_col = "SNP", beta_col = "beta.outcome", eaf_col ="MAF",se_col = "se.outcome", effect_allele_col = "EA", other_allele_col = "OA",  pval_col = "P")
exp<-clump_data(exp, clump_r2 = .3)
write.csv(exp   ,"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caff_snps_em5_r.3.csv")



#n.b. build is 37



#Finngen data 
finn2=NULL
for (i in out$FinnGenID[!is.na(out$FinnGenID)]) {
  finn<-data.table::fread(paste("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\summary_stats_finngen_R8_",i,".gz",sep = ""))
  finn<-finn[finn$rsids %in% exp$SNP,]
  finn$outcome=i
  finn2=rbind(finn,finn2)
}
finn2$id.outcome=finn2$outcome
write.csv(finn2,"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\caff_finngenOUT_snps.csv")


# getting open GWAS data
ieu<-extract_outcome_data(exp$SNP,out$OpenGWAS_ID[!is.na(out$OpenGWAS_ID)],proxies = F)
#ieu<-extract_outcome_data(exp$SNP,c("ieu-b-4877",'ieu-a-1239','ukb-b-5779',"ukb-b-4424","ieu-b-5099","ieu-b-102","ukb-b-5238","ieu-a-1183","ieu-a-1185"),proxies = F)
write.csv(ieu,"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\caff_ieuOUT_snps.csv")


# other daa sources 
other2<-fread("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\GCST90104339_buildGRCh37.tsv")
other2<-other2[other2$rsid %in% exp$SNP,]
other2 <- format_data(other2, type = "outcome", snp_col = "rsid", beta_col = "BETA", se_col = "SE", effect_allele_col = "Allele2", other_allele_col = "Allele1", eaf_col = "Freq1", pval_col = "p_value")
other2$outcome<-"GCST90104339"
other2$id.outcome<-"GCST90104339"
other<-other[colnames(other) %in% colnames(other2)]
other2<-other2[colnames(other2) %in% colnames(other)]
other<-other2

other2<-fread("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\GCST90104341_buildGRCh37.tsv")
other2<-other2[other2$rsid %in% exp$SNP,]
other2 <- format_data(other2, type = "outcome", snp_col = "rsid", beta_col = "BETA", se_col = "SE", effect_allele_col = "Allele2", other_allele_col = "Allele1", eaf_col = "Freq1", pval_col = "p_value")
other2$outcome<-"GCST90104341"
other2$id.outcome<-"GCST90104341"
#other<-other[colnames(other) %in% colnames(other2)]
other2<-other2[colnames(other2) %in% colnames(other)]
other<-rbind(other,other2)

other2<-fread("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\pgcAN2.2019-07.vcf.tsv.gz")
other2<-other2[other2$ID %in% exp$SNP,]
other2 <- format_data(other2, type = "outcome", snp_col = "ID", beta_col = "BETA", se_col = "SE", effect_allele_col = "ALT", other_allele_col = "REF",  pval_col = "PVAL")
other2$eaf.outcome<-NA
other2$outcome<-"31308545"
other2$id.outcome<-"31308545"
#other<-other[colnames(other) %in% colnames(other2)]
other2<-other2[colnames(other2) %in% colnames(other)]
other<-rbind(other,other2)

other2<-fread("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\pgc-bip2021-all.vcf.tsv.gz")
other2<-other2[other2$ID %in% exp$SNP,]
other2 <- format_data(other2, type = "outcome", snp_col = "ID", beta_col = "BETA", se_col = "SE", effect_allele_col = "A1", other_allele_col = "A2",  pval_col = "PVAL")
#not sure if A1 and A2 are right
other2$eaf.outcome<-NA
other2$outcome<-"34002096"
other2$id.outcome<-"34002096"
#other<-other[colnames(other) %in% colnames(other2)]
other2<-other2[colnames(other2) %in% colnames(other)]
other<-rbind(other,other2)

write.csv(other,"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\caff_otherOUT_snps.csv")


#who cooked addam smiths dinner

#####################################
#muiptile test correction 
# 0.080806868	insomnia
# 0.646350038	schizophrenia
# 0.030635077	bipolar disorder
# 0.019428908	major depression
# 0.014682786	anorexia nervosa


#p.adjust(c(0.080806868,0.646350038,0.030635077,0.019428908	,0.014682786	),method="fdr")
p.adjust(c(0.646350038,0.030635077,0.019428908	,0.014682786	),method="fdr")


###########################
### psyc sup tabs
##############################
exp <-read.csv("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caff_snps_em5_r.3.csv",header=T)
write.csv(exp,"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\write up\\psyc_Stab_1exp.csv")

finn<-read.csv("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\caff_finngenOUT_snps.csv")
finn<-finn[finn$outcome=="F5_BIPO"|finn$outcome=="R18_ANOREXIA"|finn$outcome== "F5_DEPRESSIO"|finn$outcome==  "F5_SCHZPHR" , ]
write.csv(finn,"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\write up\\psyc_Stab_2finn.csv")

ieu<-read.csv("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\caff_ieuOUT_snps.csv")
other<-read.csv("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\caff_otherOUT_snps.csv")
ieu<-ieu[names(ieu)%in%names(other)]
other<-other[names(other)%in%names(ieu)]
finn<-rbind(ieu,other)
finn<-finn[finn$outcome=="34002096"|finn$outcome=="31308545"|finn$outcome== "Major depression || id:ieu-b-102"|finn$outcome==  "Schizophrenia || id:ieu-b-5099" , ]
write.csv(finn,"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\write up\\psyc_Stab_3out.csv")

###########
#######################################################################
#psyc
######################################################################


out<-as.data.frame(1:4)
out$pheno<-c("schizophrenia", "bipolar disorder", "major depression", "anorexia nervosa")
out$OpenGWAS_ID<-c("ieu-b-5099",NA,"ieu-b-102",NA)
out$PMIDorEBIid<-c(NA,"34002096",NA,"31308545") 
out$FinnGenID<-c("F5_SCHZPHR","F5_BIPO", "F5_DEPRESSIO", "R18_ANOREXIA")
out$outcome_units<-c("logOR","logOR","logOR","logOR")

out2<-as.data.frame(1:12)
out2$pheno<-c("major depression","major depression","major depression","bipolar disorder","bipolar disorder","bipolar disorder",  "anorexia nervosa", "anorexia nervosa", "anorexia nervosa","schizophrenia","schizophrenia","schizophrenia")
out2$OpenGWAS_ID<-c("ieu-b-102","ieu-b-102","ieu-b-102",NA,NA,NA,NA,NA,NA,"ieu-b-5099","ieu-b-5099","ieu-b-5099")
out2$PMIDorEBIid<-c(NA,NA,NA,"34002096","34002096","34002096","31308545","31308545","31308545",NA,NA,NA) 
out2$FinnGenID<-c("F5_DEPRESSIO", "F5_DEPRESSIO", "F5_DEPRESSIO","F5_BIPO","F5_BIPO","F5_BIPO",  "R18_ANOREXIA", "R18_ANOREXIA", "R18_ANOREXIA","F5_SCHZPHR","F5_SCHZPHR","F5_SCHZPHR")
out2$outcome_units<-c("logOR","logOR","logOR","logOR","logOR","logOR","logOR","logOR","logOR","logOR","logOR","logOR")
out2$gene<-c("AHR","CYP1A2","Combined","AHR","CYP1A2","Combined","AHR","CYP1A2","Combined","AHR","CYP1A2","Combined")



exp <-read.csv("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caff_snps_em5_r.3.csv",header=T)
exp$beta.exposure<-exp$beta.exposure*FIQT_deflation(B=exp$beta.exposure, SE = exp$se.exposure)
finn<-read.csv("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\caff_finngenOUT_snps.csv")

ieu<-read.csv("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\caff_ieuOUT_snps.csv")
other<-read.csv("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\caff_otherOUT_snps.csv")


matrix=ld_matrix(exp$SNP)
a<-unlist(stringr::str_split(colnames(matrix),"_"))
colnames(matrix)<-a [seq(from=1,to=length(a), by =3)] 
a<-unlist(stringr::str_split(rownames(matrix),"_"))
rownames(matrix)<-a [seq(from=1,to=length(a), by =3)] 

for (i in 1:nrow(out)) {
  if (!is.na(out$OpenGWAS_ID[i])){
    dat<-ieu[ieu$id.outcome==out$OpenGWAS_ID[i],]
    if (is.na(out$OpenGWAS_ID[i])){
      dat2<-other[other$id.outcome==out$PMIDorEBIid[i],]
      dat2<-convert_outcome_to_exposure(dat2)
      dat2<-harmonise_data(dat2,dat)
      dat2[c("beta.outcome","se.outcome","pval.outcome")]<-t(apply(dat2[c("beta.exposure","se.exposure","beta.outcome","se.outcome")],1,function(x){met<-meta::metagen(TE = unlist(x[c(1,3)]), seTE = unlist(x[c(2,4)]),random = F);return(c(met$TE.fixed,met$seTE.fixed,met$pval.fixed))}))
      dat<-dat2[colnames(dat)] 
    }
    if(!is.na(out$FinnGenID[i])){
      dat2<-finn[finn$outcome==out$FinnGenID[i],]
      dat2 <- format_data(dat2, type = "outcome",snp_col = "rsids", beta_col = "beta", se_col = "sebeta", effect_allele_col = "alt", other_allele_col = "ref", eaf_col = "af_alt", pval_col = "pval")
      
      dat2<-convert_outcome_to_exposure(dat2)
      dat2<-harmonise_data(dat2,dat)
      dat2[c("beta.outcome","se.outcome","pval.outcome")]<-t(apply(dat2[c("beta.exposure","se.exposure","beta.outcome","se.outcome")],1,function(x){met<-meta::metagen(TE = unlist(x[c(1,3)]), seTE = unlist(x[c(2,4)]),random = F);return(c(met$TE.fixed,met$seTE.fixed,met$pval.fixed))}))
      dat<-dat2[colnames(dat)]                                           
    }
  }
  
  if (is.na(out$OpenGWAS_ID[i])){
    dat<-other[other$id.outcome==out$PMIDorEBIid[i],]
    if(!is.na(out$FinnGenID[i])){
      dat2<-finn[finn$outcome==out$FinnGenID[i],]
      dat2 <- format_data(dat2, type = "outcome",snp_col = "rsids", beta_col = "beta", se_col = "sebeta", effect_allele_col = "alt", other_allele_col = "ref", eaf_col = "af_alt", pval_col = "pval")
      
      dat2<-convert_outcome_to_exposure(dat2)
      dat2<-harmonise_data(dat2,dat)
      dat2[c("beta.outcome","se.outcome","pval.outcome")]<-t(apply(dat2[c("beta.exposure","se.exposure","beta.outcome","se.outcome")],1,function(x){met<-meta::metagen(TE = unlist(x[c(1,3)]), seTE = unlist(x[c(2,4)]),random = F);return(c(met$TE.fixed,met$seTE.fixed,met$pval.fixed))}))
      dat<-dat2[colnames(dat)]                                           
    }
  } 
  
  
  dat<-harmonise_data(exp[exp$SNP %in% dat$SNP,],dat)
  
  
  #leave one out analysis 
  jpeg(paste("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\LOO",out$pheno[i],".jpg",sep = ""))
  cisLeave1out(SNP=dat$SNP,
               beta.exposure=dat$beta.exposure,
               se.exposure=dat$se.exposure,  
               beta.outcome=dat$beta.outcome,
               se.outcome=dat$se.outcome,
               sm="OR",
               gene=dat$gene.exposure,
               matrix=matrix)
  dev.off()  
  
  #leave one out
  jpeg(paste("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\LOO",out$pheno[i],"_combined.jpg",sep = ""))
  cisLeave1out(SNP=dat$SNP,
               beta.exposure=dat$beta.exposure,
               se.exposure=dat$se.exposure,  
               beta.outcome=dat$beta.outcome,
               se.outcome=dat$se.outcome,
               sm="OR",
               #gene=dat$gene.exposure,
               matrix=matrix)
  dev.off()
  
  #primary analysis 
  
  IVW<-IVWcorrel( #dbp on no children  
    betaYG=dat$beta.outcome,
    sebetaYG=dat$se.outcome,
    betaXG=dat$beta.exposure,
    sebetaXG=dat$se.exposure,
    rho=matrix[dat$SNP,dat$SNP]
  )
  out2$n_snp[out2$FinnGenID==out$FinnGenID[i]&out2$gene=="Combined"]<-as.numeric(IVW[2,4])#res$nsnp
  out2$beta[out2$FinnGenID==out$FinnGenID[i]&out2$gene=="Combined"]<-as.numeric(IVW[2,1])#res$b
  out2$se[out2$FinnGenID==out$FinnGenID[i]&out2$gene=="Combined"]<-as.numeric(IVW[2,2])#res$se
  out2$pval[out2$FinnGenID==out$FinnGenID[i]&out2$gene=="Combined"]<-as.numeric(IVW[2,3])#res$pval
  out2$MR_F[out2$FinnGenID==out$FinnGenID[i]&out2$gene=="Combined"]<-as.numeric(IVW[2,5])
  
  IVW<-IVWcorrel( #dbp on no children  
    betaYG=dat$beta.outcome[dat$gene.exposure=="AHR"],
    sebetaYG=dat$se.outcome[dat$gene.exposure=="AHR"],
    betaXG=dat$beta.exposure[dat$gene.exposure=="AHR"],
    sebetaXG=dat$se.exposure[dat$gene.exposure=="AHR"],
    rho=matrix[dat$SNP[dat$gene.exposure=="AHR"],dat$SNP[dat$gene.exposure=="AHR"]]
  )
  out2[out2$FinnGenID==out$FinnGenID[i]&out2$gene=="AHR","n_snp"]<-as.numeric(IVW[2,4])#res$nsnp
  out2[out2$FinnGenID==out$FinnGenID[i]&out2$gene=="AHR","beta"]<-as.numeric(IVW[2,1])#res$b
  out2[out2$FinnGenID==out$FinnGenID[i]&out2$gene=="AHR","se"]<-as.numeric(IVW[2,2])#res$se
  out2[out2$FinnGenID==out$FinnGenID[i]&out2$gene=="AHR","pval"]<-as.numeric(IVW[2,3])#res$pval
  out2[out2$FinnGenID==out$FinnGenID[i]&out2$gene=="AHR","MR_F"]<-as.numeric(IVW[2,5])

  IVW<-IVWcorrel( #dbp on no children  
    betaYG=dat$beta.outcome[dat$gene.exposure=="CYP1A2"],
    sebetaYG=dat$se.outcome[dat$gene.exposure=="CYP1A2"],
    betaXG=dat$beta.exposure[dat$gene.exposure=="CYP1A2"],
    sebetaXG=dat$se.exposure[dat$gene.exposure=="CYP1A2"],
    rho=matrix[dat$SNP[dat$gene.exposure=="CYP1A2"],dat$SNP[dat$gene.exposure=="CYP1A2"]]
  )
  out2[out2$FinnGenID==out$FinnGenID[i]&out2$gene=="CYP1A2","n_snp"]<-as.numeric(IVW[2,4])#res$nsnp
  out2[out2$FinnGenID==out$FinnGenID[i]&out2$gene=="CYP1A2","beta"]<-as.numeric(IVW[2,1])#res$b
  out2[out2$FinnGenID==out$FinnGenID[i]&out2$gene=="CYP1A2","se"]<-as.numeric(IVW[2,2])#res$se
  out2[out2$FinnGenID==out$FinnGenID[i]&out2$gene=="CYP1A2","pval"]<-as.numeric(IVW[2,3])#res$pval
  out2[out2$FinnGenID==out$FinnGenID[i]&out2$gene=="CYP1A2","MR_F"]<-as.numeric(IVW[2,5])
  
  
  # using lead SNP only 
  out2[out2$FinnGenID==out$FinnGenID[i]&out2$gene=="AHR","beta_leadSNP"]<-dat$beta.outcome[dat$SNP=="rs4410790"]/dat$beta.exposure[dat$SNP=="rs4410790"]
  out2[out2$FinnGenID==out$FinnGenID[i]&out2$gene=="AHR","se_leadSNP"]<-abs(dat$se.outcome[dat$SNP=="rs4410790"]/dat$beta.exposure[dat$SNP=="rs4410790"])
  out2[out2$FinnGenID==out$FinnGenID[i]&out2$gene=="CYP1A2","beta_leadSNP"]<-dat$beta.outcome[dat$SNP=="rs2472297"]/dat$beta.exposure[dat$SNP=="rs2472297"]
  out2[out2$FinnGenID==out$FinnGenID[i]&out2$gene=="CYP1A2","se_leadSNP"]<-abs(dat$se.outcome[dat$SNP=="rs2472297"]/dat$beta.exposure[dat$SNP=="rs2472297"])
  IVW<-IVWcorrel( #dbp on no children  
    betaYG=dat$beta.outcome[dat$SNP=="rs4410790"|dat$SNP=="rs2472297"],
    sebetaYG=dat$se.outcome[dat$SNP=="rs4410790"|dat$SNP=="rs2472297"],
    betaXG=dat$beta.exposure[dat$SNP=="rs4410790"|dat$SNP=="rs2472297"],
    sebetaXG=dat$se.exposure[dat$SNP=="rs4410790"|dat$SNP=="rs2472297"],
    rho=matrix[dat$SNP[dat$SNP=="rs4410790"|dat$SNP=="rs2472297"],dat$SNP[dat$SNP=="rs4410790"|dat$SNP=="rs2472297"]]
  )
  
  out2[out2$FinnGenID==out$FinnGenID[i]&out2$gene=="Combined","beta_leadSNP"]<-as.numeric(IVW[2,1])#summary(lm(beta.outcome~beta.exposure+0,weights = se.outcome^-2, data = dat[dat$SNP=="rs4410790"|dat$SNP=="rs2472297",]))$coef[1]
  out2[out2$FinnGenID==out$FinnGenID[i]&out2$gene=="Combined","se_leadSNP"]<-as.numeric(IVW[2,2])#summary(lm(beta.outcome~beta.exposure+0,weights = se.outcome^-2, data = dat[dat$SNP=="rs4410790"|dat$SNP=="rs2472297",]))$coef[2]
  
}
write.csv(out2,"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\caff_psyc_mrout.csv")


#plot figgure for the psyc ms
output<-read.csv("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\caff_psyc_mrout.csv")

#output<-output[output$pheno=="bipolar disorder"|output$pheno=="anorexia nervosa"|output$pheno== "major depression"|output$pheno==  "schizophrenia" , ]
output$pheno[output$pheno=="bipolar disorder"]<-"Bipolar disorder"
output$pheno[output$pheno=="anorexia nervosa"]<-"Anorexia nervosa"
output$pheno[output$pheno== "major depression"]<- "Major depression"
output$pheno[output$pheno==  "schizophrenia"]<- "Schizophrenia"
output$Outcome<-output$pheno
#meta::forest(meta::metagen(data=output,TE=beta,overall=FALSE, seTE =se , random = F,fixed=F ,studlab = pheno,  sm="OR",  tau.common = FALSE), lwd=1.5, col.square="black",plotwidth="12cm", squaresize = .2,  weight.study = "same",  leftcols = c("studlab"), leftlabs = c(""), rightlabs = c("Effect", "[ 95% CI ]"),smlab="")
meta::forest(meta::metagen(data=output,TE=beta,overall=FALSE, seTE =se ,  test.subgroup=FALSE, random = F,fixed=F ,studlab = gene,  sm="OR", subgroup = Outcome, tau.common = FALSE), lwd=1.5, col.square="black",plotwidth="12cm", squaresize = .2,  weight.study = "same",  leftcols = c("studlab"), leftlabs = c(""), rightlabs = c("Odds ratio", "[ 95% CI ]"),smlab="",sep.subgroup = "",subgroup.name ="")

#wiht lead SNP only
meta::forest(meta::metagen(data=output,TE=beta_leadSNP,overall=FALSE, seTE =se_leadSNP  ,  test.subgroup=FALSE, random = F,fixed=F ,studlab = gene,  sm="OR", subgroup = Outcome, tau.common = FALSE), lwd=1.5, col.square="black",plotwidth="12cm", squaresize = .2,  weight.study = "same",  leftcols = c("studlab"), leftlabs = c(""), rightlabs = c("Odds ratio", "[ 95% CI ]"),smlab="",sep.subgroup = "",subgroup.name ="")

#estimating hetrgenaty
meta::forest(meta::metagen(data=output[output$gene!="Combined",],TE=beta,overall=T, seTE =se  ,  test.subgroup=T, random = T,fixed=T ,studlab = gene,  sm="OR", subgroup = Outcome, tau.common = FALSE))


################
# study level het
###################
exp <-read.csv("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caff_snps_em5_r.3.csv",header=T)
exp$beta.exposure<-exp$beta.exposure*FIQT_deflation(B=exp$beta.exposure, SE = exp$se.exposure)
finn<-read.csv("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\caff_finngenOUT_snps.csv")
finn <- format_data(finn, type = "outcome",snp_col = "rsids", beta_col = "beta", se_col = "sebeta", effect_allele_col = "alt", other_allele_col = "ref", eaf_col = "af_alt", pval_col = "pval",  id_col= "id.outcome",phenotype_col ="id.outcome")

ieu<-read.csv("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\caff_ieuOUT_snps.csv")
other<-read.csv("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\caff_otherOUT_snps.csv")

matrix=ld_matrix(exp$SNP)
a<-unlist(stringr::str_split(colnames(matrix),"_"))
colnames(matrix)<-a [seq(from=1,to=length(a), by =3)] 
a<-unlist(stringr::str_split(rownames(matrix),"_"))
rownames(matrix)<-a [seq(from=1,to=length(a), by =3)] 

out<-as.data.frame(1:4)
out$pheno<-c("Anorexia nervosa", "Bipolar disorder", "Major depression","Schizophrenia" )
out$OpenGWAS_ID<-c(NA,NA,"ieu-b-102","ieu-b-5099")
out$PMIDorEBIid<-c("31308545","34002096",NA,NA) 
out$FinnGenID<-c( "R18_ANOREXIA","F5_BIPO", "F5_DEPRESSIO","F5_SCHZPHR")
out$outcome_units<-c("logOR","logOR","logOR","logOR")

dat<-harmonise_data(exp,finn)

for (i in unique(dat$outcome)){
  dat2<-dat[dat$outcome==i,]
IVW<-IVWcorrel( #dbp on no children  
  betaYG=dat2$beta.outcome,
  sebetaYG=dat2$se.outcome,
  betaXG=dat2$beta.exposure,
  sebetaXG=dat2$se.exposure,
  rho=matrix[dat2$SNP,dat2$SNP]
)
out$beta[out$FinnGenID==i]<-as.numeric(IVW[2,1])#res$b
out$se[out$FinnGenID==i]<-as.numeric(IVW[2,2])#res$se
}

out2<-out
out2$beta<-NA

dat<-harmonise_data(exp,ieu)
for (i in c("ieu-b-5099","ieu-b-102")){
  dat2<-dat[dat$id.outcome==i,]
  IVW<-IVWcorrel( #dbp on no children  
    betaYG=dat2$beta.outcome,
    sebetaYG=dat2$se.outcome,
    betaXG=dat2$beta.exposure,
    sebetaXG=dat2$se.exposure,
    rho=matrix[dat2$SNP,dat2$SNP]
  )
  out2$beta[out2$OpenGWAS_ID ==i]<-as.numeric(IVW[2,1])#res$b
  out2$se[out2$OpenGWAS_ID ==i]<-as.numeric(IVW[2,2])#res$se
}

dat<-harmonise_data(exp,other)
for (i in unique(dat$id.outcome)){
  dat2<-dat[dat$id.outcome==i,]
  IVW<-IVWcorrel( #dbp on no children  
    betaYG=dat2$beta.outcome,
    sebetaYG=dat2$se.outcome,
    betaXG=dat2$beta.exposure,
    sebetaXG=dat2$se.exposure,
    rho=matrix[dat2$SNP,dat2$SNP]
  )
  out2$beta[out2$PMIDorEBIid ==i]<-as.numeric(IVW[2,1])#res$b
  out2$se[out2$PMIDorEBIid ==i]<-as.numeric(IVW[2,2])#res$se
}
out<-rbind(out,out2)
out$GWAS<-c("FinnGen","FinnGen","FinnGen","FinnGen","PGC","PGC + UKB","PGC + 23andMe + UKB","PGC")

out$Outcome<-out$pheno
meta::forest(meta::metagen(data=out,TE=beta,overall=FALSE, seTE =se ,  test.subgroup=T, random = T,fixed=T ,studlab = GWAS,  sm="OR", subgroup = Outcome, tau.common = F), smlab="",sep.subgroup = "",subgroup.name ="")#, lwd=1.5, col.square="black",plotwidth="12cm", squaresize = .2,  weight.study = "same",  leftcols = c("studlab"), leftlabs = c(""), rightlabs = c("Odds ratio", "[ 95% CI ]"),smlab="",sep.subgroup = "",subgroup.name ="")

