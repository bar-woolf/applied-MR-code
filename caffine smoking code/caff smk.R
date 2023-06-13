library(data.table)
library(coloc)
library(TwoSampleMR)
library(ggplot2)
#library(gwasglue)
#library(gassocplot)
#library(MRSamePopTest)
library(TwoStepCisMR)


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



#  dtaa sources 
other<-read_outcome_data(  "C:/Users/nf18217/OneDrive - University of Bristol/Documents/Projects/active/passive smoking/Lifetime_Smoking_Sheet1.txt", snps=exp$SNP, chr_col = "CHR",  sep = " ",  snp_col = "SNP",  beta_col = "BETA",  se_col = "SE",  eaf_col = "EAF",  effect_allele_col = "EFFECT_ALLELE",  other_allele_col = "OTHER_ALLELE",pval_col = "P")
other$outcome<-"31689377"
other$id.outcome<-"31689377"
write.csv(other,"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\caff_otherOUT_snps.csv")


#who cooked addam smiths dinner

#################################################################
# stieger filtering for the smoking anlasyis 

exp <-read.csv("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caff_snps_em5_r.3.csv",header=T)
exp$beta.exposure<-exp$beta.exposure*FIQT_deflation(B=exp$beta.exposure, SE = exp$se.exposure)
other<-read.csv("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\caff_otherOUT_snps.csv")
other<-other[other$outcome=="31689377",]

dat<-harmonise_data(exp,other)
dat$samplesize.exposure<-9876 
dat$samplesize.outcome<-462690 
directionality_test(dat)




########################################################################
# smoking 
#####################################################################

#plot figure for smoking ms
other<-read.csv("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\caff_otherOUT_snps.csv")
other<-other[other$outcome=="31689377",]


exp <-read.csv("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caff_snps_em5_r.3.csv",header=T)
exp$beta.exposure<-exp$beta.exposure*FIQT_deflation(B=exp$beta.exposure, SE = exp$se.exposure)
dat<-harmonise_data(exp,other)
write.csv(dat,"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caffine and other outcomes\\write up\\Smk_Stab_1.csv")

matrix=ld_matrix(exp$SNP)
a<-unlist(stringr::str_split(colnames(matrix),"_"))
colnames(matrix)<-a [seq(from=1,to=length(a), by =3)] 
a<-unlist(stringr::str_split(rownames(matrix),"_"))
rownames(matrix)<-a [seq(from=1,to=length(a), by =3)] 

out<-as.data.frame(1:3)
out$pheno<-c("AHR","CYP1A2","Combined")

IVW<-IVWcorrel( #dbp on no children  
  betaYG=dat$beta.outcome[dat$gene.exposure=="AHR"],
  sebetaYG=dat$se.outcome[dat$gene.exposure=="AHR"],
  betaXG=dat$beta.exposure[dat$gene.exposure=="AHR"],
  sebetaXG=dat$se.exposure[dat$gene.exposure=="AHR"],
  rho=matrix[dat$SNP[dat$gene.exposure=="AHR"],dat$SNP[dat$gene.exposure=="AHR"]]
)
out[out$pheno=="AHR",IVW[1,]]<-as.numeric(IVW[2,])

IVW<-IVWcorrel( #dbp on no children  
  betaYG=dat$beta.outcome[dat$gene.exposure=="CYP1A2"],
  sebetaYG=dat$se.outcome[dat$gene.exposure=="CYP1A2"],
  betaXG=dat$beta.exposure[dat$gene.exposure=="CYP1A2"],
  sebetaXG=dat$se.exposure[dat$gene.exposure=="CYP1A2"],
  rho=matrix[dat$SNP[dat$gene.exposure=="CYP1A2"],dat$SNP[dat$gene.exposure=="CYP1A2"]]
)
out[out$pheno=="CYP1A2",IVW[1,]]<-as.numeric(IVW[2,])

IVW<-IVWcorrel( #dbp on no children  
  betaYG=dat$beta.outcome,
  sebetaYG=dat$se.outcome,
  betaXG=dat$beta.exposure,
  sebetaXG=dat$se.exposure,
  rho=matrix[dat$SNP,dat$SNP]
)
out[out$pheno=="Combined",IVW[1,]]<-as.numeric(IVW[2,])
out$se_IVWcorrel.random<-out$se_IVWcorrel.random/0.6940093 #standardising outcome
out$beta_IVWcorrel<-out$beta_IVWcorrel/0.6940093
meta::forest(meta::metagen(data=out,TE=beta_IVWcorrel,overall=FALSE, seTE =se_IVWcorrel.random , random = F,fixed=F ,studlab = pheno,  sm="md",  tau.common = FALSE), lwd=1.5, col.square="black",plotwidth="12cm", squaresize = .2,  weight.study = "same",  leftcols = c("studlab"), leftlabs = c(""), rightlabs = c("Effect", "[ 95% CI ]"),smlab="",digits=3)


#mr egger 

EGGERcorrel(
  betaYG=dat$beta.outcome,
  sebetaYG=dat$se.outcome,
  betaXG=dat$beta.exposure,
  sebetaXG=dat$se.exposure,
  rho=matrix[dat$SNP,dat$SNP]
)

#leave one out
dat2<-dat
dat2$gene.exposure<-"Combined"
dat2<-rbind(dat,dat2)
cisLeave1out(SNP=dat2$SNP,
             beta.exposure=dat2$beta.exposure,
             se.exposure=dat2$se.exposure,  
             beta.outcome=dat2$beta.outcome,
             se.outcome=dat2$se.outcome,
             sm="md",
             gene=dat2$gene.exposure,
             matrix=matrix)

#ever somking 
other<-extract_outcome_data(exp$SNP, "ieu-b-4877") #data extraction for ever smoking. Unhash to run this analsyis. 
dat<-harmonise_data(exp,other)
IVW<-IVWcorrel( #dbp on no children  
  betaYG=dat$beta.outcome,
  sebetaYG=dat$se.outcome,
  betaXG=dat$beta.exposure,
  sebetaXG=dat$se.exposure,
  rho=matrix[dat$SNP,dat$SNP]
)

