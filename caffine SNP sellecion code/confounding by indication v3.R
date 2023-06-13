





library(TwoSampleMR)
library(MVMR)

dat<-data.table::fread("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\caff_share.gz")
# beta = Z-score / sqrt(2 * alleleFreq * (1-alleleFreq) * (N + Z-score^2)
# se = 1 / sqrt(2 * alleleFreq * (1-alleleFreq) * (N + Z-score^2)

##################################
dat2<-dat[dat$P<5*10^-7,]
dat2$phenoscanner<-paste("chr",dat2$`chr:pos`,"",sep="")
#n.b. build is 37

write.table(dat2$phenoscanner,"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\chrpos.txt",row.names = F,col.names = F,quote = F)
id<-data.table::fread("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\chrpos_PhenoScanner_SNP.tsv")
id<-as.data.frame(id)[c("snp","rsid","eur")]
colnames(id)<-c("phenoscanner","SNP","MAF")
dat2<-merge(dat2,id,by="phenoscanner")
dat2$pval.exposure<-dat2$P

write.csv(dat2,"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\sig_caff_snps.csv")

######################################################################################
#########################################################################################
# tea 
######################################################################################
######################################################################################


#####################################################
# stieger filtering for cis SNPs
#################################


# drink vs caf dir
################################


#getting the data for UVMR
#tea<-extract_instruments("ukb-b-5237", p1 = 5e-07) #tea ukb-b-6066, coffe ukb-b-5237
tea<-extract_instruments("ukb-b-6066", p1 = 5e-07) #tea ukb-b-6066, coffe ukb-b-5237
tea<-harmonise_data(tea,extract_outcome_data(tea$SNP,"ieu-b-40",proxies = F)) #t2d ebi-a-GCST006867, bmi  ieu-b-40, az ebi-a-GCST90012877

tea$`chr:pos`<-paste("",tea$chr.exposure,":",tea$pos.exposure,"",sep="")
dat_tea<-dat[dat$`chr:pos` %in% tea$`chr:pos`,]
dat_tea<-merge(tea,dat_tea,by="chr:pos")
dat_tea$beta.caff<-dat_tea$Effect/(2*dat_tea$eaf.exposure*(1-dat_tea$eaf.exposure)*(9876 +dat_tea$Effect^2))^0.5
dat_tea$se.caff<-1/(2*dat_tea$eaf.exposure*(1-dat_tea$eaf.exposure)*(9876 +dat_tea$Effect^2))^0.5
dat_tea$EA[dat_tea$EA=="a"]<-"A"
dat_tea$EA[dat_tea$EA=="g"]<-"G"
dat_tea$EA[dat_tea$EA=="t"]<-"T"
dat_tea$EA[dat_tea$EA=="c"]<-"C"
dat_tea$beta.caff[dat_tea$effect_allele.exposure!=dat_tea$EA]<- -1 *dat_tea$beta.caff[dat_tea$effect_allele.exposure!=dat_tea$EA]

#stieger filtering
dat_tea$r2_tea<-2*dat_tea$eaf.exposure*(1-dat_tea$eaf.exposure)*dat_tea$beta.exposure^2
dat_tea$r2_caff<- 2*dat_tea$eaf.exposure*(1-dat_tea$eaf.exposure)*dat_tea$beta.caff^2# dat_tea$Effect^2/(9876-2+dat_tea$Effect^2)#
nrow(dat_tea)
nrow(dat_tea[dat_tea$r2_tea<dat_tea$r2_caff,]) # only 21 out of 41 SNPs primerily target consumption for tea.  for coffe it is 21 out of 56
# write.table(dat_tea$SNP[dat_tea$r2_tea>dat_tea$r2_caff],"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\snp sellecion\\coffe_snps_stigerdrinking.txt",sep="\t",quote=F,col.names = F,row.names = F)
# write.table(dat_tea$SNP[dat_tea$r2_tea<dat_tea$r2_caff],"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\snp sellecion\\coffe_snps_stigerplasma.txt",sep="\t",quote=F,col.names = F,row.names = F)
write.table(dat_tea$SNP[dat_tea$r2_tea>dat_tea$r2_caff],"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\snp sellecion\\tea_snps_stigerdrinking.txt",sep="\t",quote=F,col.names = F,row.names = F)
write.table(dat_tea$SNP[dat_tea$r2_tea<dat_tea$r2_caff],"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\snp sellecion\\tea_snps_stigerplasma.txt",sep="\t",quote=F,col.names = F,row.names = F)


#tea on BMI 
mean(dat_tea$beta.exposure^2/dat_tea$beta.outcome^2)
summary(lm(beta.outcome~0+beta.exposure,weights = se.outcome^-2, data = dat_tea)) #no stieger filtering

#cis mr
summary(lm(beta.outcome~0+beta.caff,weights = se.outcome^-2, data = dat_tea[dat_tea$SNP=="rs2472297"|dat_tea$SNP=="rs4410790",])) #no stieger filtering
mean(dat_tea$beta.caff[dat_tea$SNP=="rs2472297"|dat_tea$SNP=="rs4410790"]^2/dat_tea$se.caff[dat_tea$SNP=="rs2472297"|dat_tea$SNP=="rs4410790"]^2)
sum(dat_tea$r2_tea[dat_tea$SNP=="rs2472297"|dat_tea$SNP=="rs4410790"])
sum(dat_tea$r2_caff[dat_tea$SNP=="rs2472297"|dat_tea$SNP=="rs4410790"])

#IVW in stigher tea drinking
summary(lm(beta.caff~0+beta.exposure,weights = se.caff^-2, data = dat_tea[dat_tea$r2_tea>dat_tea$r2_caff,])) #stieger filtering
SNPdrink<-dat_tea$SNP[dat_tea$r2_tea>dat_tea$r2_caff]
summary(lm(beta.outcome~0+beta.exposure,weights = se.outcome^-2, data = dat_tea[dat_tea$SNP%in%SNPdrink,])) #stieger filtering

#IVW in stiger plasam caff
summary(lm(beta.exposure~0+beta.caff,weights = se.exposure^-2, data = dat_tea[dat_tea$r2_tea<dat_tea$r2_caff,])) #stieger filtering
SNPcaff<-dat_tea$SNP[dat_tea$r2_tea<dat_tea$r2_caff]
summary(lm(beta.outcome~0+beta.exposure,weights = se.outcome^-2, data = dat_tea[dat_tea$SNP%in%SNPcaff,])) #stieger filtering, n.b. higher scores of these SNPs predict lower circ caffine so this is in agreement with the cis study
summary(lm(beta.outcome~0+beta.caff,weights = se.outcome^-2, data = dat_tea[dat_tea$SNP%in%SNPcaff,])) #stieger filtering, n.b. higher scores of these SNPs predict lower circ caffine so this is in agreement with the cis study
mean(dat_tea$beta.caff[dat_tea$r2_tea<dat_tea$r2_caff]^2/dat_tea$se.caff[dat_tea$r2_tea<dat_tea$r2_caff]^2)
mean(dat_tea$beta.outcome[dat_tea$r2_tea<dat_tea$r2_caff]^2/dat_tea$se.outcome[dat_tea$r2_tea<dat_tea$r2_caff]^2)
mean(dat_tea$beta.exposure[dat_tea$r2_tea<dat_tea$r2_caff]^2/dat_tea$se.exposure[dat_tea$r2_tea<dat_tea$r2_caff]^2)

#hetrgneaty in the behavoural SNPs
res<-mr(dat_tea[dat_tea$SNP%in%SNPdrink,])
res<-rbind(res,mr(dat_tea))
res$SNP<-c("Steiger tea SNPs","Steiger tea SNPs","Steiger tea SNPs","Steiger tea SNPs","Steiger tea SNPs","no Steiger filter","no Steiger filter","no Steiger filter","no Steiger filter","no Steiger filter")
write.csv(res,"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\snp sellecion\\SuppTableS3a.csv")
#Isq(bmi$beta.outcome[bmi$SNP%in%SNPdrink]/bmi$beta.exposure[bmi$SNP%in%SNPdrink],bmi$se.outcome[bmi$SNP%in%SNPdrink]/abs(bmi$beta.exposure[bmi$SNP%in%SNPdrink]))



############################################################

#MVMR
##################

#get data above

#extract tea data from caff sig SNPs, convert to outcome, repeate above, merge thsi and the above data sets. 

id<-"ukb-b-6066"#tea ukb-b-6066, coffe ukb-b-5237
tea<-extract_instruments(id, p1 = 5e-07, clump = F) 
tea<-tea[c("SNP","pval.exposure")]

dat2<-read.csv("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\sig_caff_snps.csv")
dat2$SNP<-dat2$SNP.x
dat2<-dat2[dat2$P<5e-07,]
dat2$pval.exposure<-dat2$P
dat2<-dat2[c("SNP","pval.exposure")]

dat2<-rbind(dat2,tea)
dat2<-clump_data(dat2,  clump_p1 = 5e-07, clump_r2 = 0.001, clump_kb = 10000,)

tea<-extract_outcome_data(dat2$SNP,id,proxies = F) 
tea<-convert_outcome_to_exposure(tea)
tea$`chr:pos`<-paste("",tea$chr.exposure,":",tea$pos.exposure,"",sep="")

dat_tea<-dat[dat$`chr:pos` %in% tea$`chr:pos`,]
dat_tea<-merge(tea,dat_tea,by="chr:pos")
dat_tea$beta.caff<-dat_tea$Effect/(2*dat_tea$eaf.exposure*(1-dat_tea$eaf.exposure)*(9876 +dat_tea$Effect^2))^0.5
dat_tea$se.caff<-1/(2*dat_tea$eaf.exposure*(1-dat_tea$eaf.exposure)*(9876 +dat_tea$Effect^2))^0.5
dat_tea$EA[dat_tea$EA=="a"]<-"A"
dat_tea$EA[dat_tea$EA=="g"]<-"G"
dat_tea$EA[dat_tea$EA=="t"]<-"T"
dat_tea$EA[dat_tea$EA=="c"]<-"C"
dat_tea$beta.caff[dat_tea$effect_allele.exposure!=dat_tea$EA]<- -1 *dat_tea$beta.caff[dat_tea$effect_allele.exposure!=dat_tea$EA]

bmi<-harmonise_data(dat_tea[colnames(tea)],extract_outcome_data(dat_tea$SNP,"ieu-b-40",proxies = F)) #t2d ebi-a-GCST006867, bmi  ieu-b-40, az ebi-a-GCST90012877
bmi<-bmi[c("SNP","beta.outcome","se.outcome")]
dat_tea<-merge(dat_tea,bmi,by="SNP")

cor(dat_tea$beta.exposure,dat_tea$beta.caff)^2
summary(lm(beta.outcome~0+beta.exposure+beta.caff,weights = se.outcome^-2, data = dat_tea))
# summary(lm(beta.outcome~beta.exposure+beta.caff,weights = se.outcome^-2, data = dat_tea)) 
# 
# MendelianRandomization::mr_mvmedian(MendelianRandomization::mr_mvinput(bx = cbind(dat_tea$beta.exposure, dat_tea$beta.caff), bxse = cbind(dat_tea$se.exposure, dat_tea$se.caff),
#                                                                        by = dat_tea$beta.outcome, byse = dat_tea$se.outcome), iterations = 1000)
# MVMRmode::mv_mrmode(Bout=dat_tea$beta.outcome, 
#                     Bexp =cbind(dat_tea$beta.exposure, dat_tea$beta.caff),
#                     SEout =dat_tea$se.outcome, 
#                     SEexp =cbind(dat_tea$se.exposure, dat_tea$se.caff), CIMin = -10, CIMax =10, CIStep=0.001  )


F.data<-format_mvmr( BXGs = dat_tea[,c("beta.exposure","beta.caff")], BYG = dat_tea$beta.outcome,
                     seBXGs = dat_tea[,c("se.exposure","se.caff")],   seBYG = dat_tea$se.outcome,  RSID = dat_tea$SNP)

sres<-strength_mvmr(r_input=F.data) # F-statistic  8.262705  3.560528
qhet_mvmr(r_input = F.data, pcor=matrix(c(1,1 ,1 ,1), ncol = 2, nrow = 2), CI=T, iterations=50) #Exposure 1       0.27961733  0.073-0.568
qhet_mvmr(r_input = F.data, pcor=matrix(c(1,.5 ,.5 ,1), ncol = 2, nrow = 2), CI=F, iterations=1000)
qhet_mvmr(r_input = F.data, pcor=matrix(c(1,0 ,0 ,1), ncol = 2, nrow = 2), CI=F, iterations=1000)
#giving the same estiamtes under differnt assumptions


######################################################################################
#########################################################################################
# coffee 
######################################################################################
######################################################################################


#####################################################
# stieger filtering for cis SNPs
#################################


# drink vs caf dir
################################
tea<-extract_instruments("ukb-b-5237", p1 = 5e-07) #tea ukb-b-6066, coffe ukb-b-5237
#tea<-extract_instruments("ukb-b-6066", p1 = 5e-08) #tea ukb-b-6066, coffe ukb-b-5237
tea<-harmonise_data(tea,extract_outcome_data(tea$SNP,"ieu-b-40",proxies = F)) #t2d ebi-a-GCST006867, bmi  ieu-b-40, az ebi-a-GCST90012877

tea$`chr:pos`<-paste("",tea$chr.exposure,":",tea$pos.exposure,"",sep="")
dat_tea<-dat[dat$`chr:pos` %in% tea$`chr:pos`,]
dat_tea<-merge(tea,dat_tea,by="chr:pos")
dat_tea$beta.caff<-dat_tea$Effect/(2*dat_tea$eaf.exposure*(1-dat_tea$eaf.exposure)*(9876 +dat_tea$Effect^2))^0.5
dat_tea$se.caff<-1/(2*dat_tea$eaf.exposure*(1-dat_tea$eaf.exposure)*(9876 +dat_tea$Effect^2))^0.5
dat_tea$EA[dat_tea$EA=="a"]<-"A"
dat_tea$EA[dat_tea$EA=="g"]<-"G"
dat_tea$EA[dat_tea$EA=="t"]<-"T"
dat_tea$EA[dat_tea$EA=="c"]<-"C"
dat_tea$beta.caff[dat_tea$effect_allele.exposure!=dat_tea$EA]<- -1 *dat_tea$beta.caff[dat_tea$effect_allele.exposure!=dat_tea$EA]

#stieger filtering
dat_tea$r2_tea<-2*dat_tea$eaf.exposure*(1-dat_tea$eaf.exposure)*dat_tea$beta.exposure^2
dat_tea$r2_caff<- 2*dat_tea$eaf.exposure*(1-dat_tea$eaf.exposure)*dat_tea$beta.caff^2# dat_tea$Effect^2/(9876-2+dat_tea$Effect^2)#
nrow(dat_tea)
nrow(dat_tea[dat_tea$r2_tea<dat_tea$r2_caff,]) # only 21 out of 41 SNPs primerily target consumption for tea.  for coffe it is 21 out of 56
write.table(dat_tea$SNP[dat_tea$r2_tea>dat_tea$r2_caff],"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\snp sellecion\\coffe_snps_stigerdrinking.txt",sep="\t",quote=F,col.names = F,row.names = F)
write.table(dat_tea$SNP[dat_tea$r2_tea<dat_tea$r2_caff],"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\snp sellecion\\coffe_snps_stigerplasma.txt",sep="\t",quote=F,col.names = F,row.names = F)


#tea on BMI 
mean(dat_tea$beta.exposure^2/dat_tea$beta.outcome^2)
summary(lm(beta.outcome~0+beta.exposure,weights = se.outcome^-2, data = dat_tea)) #no stieger filtering

# #cis mr
# summary(lm(beta.outcome~0+beta.caff,weights = se.outcome^-2, data = dat_tea[dat_tea$SNP=="rs2472297"|dat_tea$SNP=="rs4410790",])) #no stieger filtering
# mean(dat_tea$beta.caff[dat_tea$SNP=="rs2472297"|dat_tea$SNP=="rs4410790"]^2/dat_tea$se.caff[dat_tea$SNP=="rs2472297"|dat_tea$SNP=="rs4410790"]^2)
# sum(dat_tea$r2_tea[dat_tea$SNP=="rs2472297"|dat_tea$SNP=="rs4410790"])
# sum(dat_tea$r2_caff[dat_tea$SNP=="rs2472297"|dat_tea$SNP=="rs4410790"])

#IVW in stigher tea drinking
summary(lm(beta.caff~0+beta.exposure,weights = se.caff^-2, data = dat_tea[dat_tea$r2_tea>dat_tea$r2_caff,])) #stieger filtering
SNPdrink<-dat_tea$SNP[dat_tea$r2_tea>dat_tea$r2_caff]
summary(lm(beta.outcome~0+beta.exposure,weights = se.outcome^-2, data = dat_tea[dat_tea$SNP%in%SNPdrink,])) #stieger filtering

#IVW in stiger plasam caff
summary(lm(beta.exposure~0+beta.caff,weights = se.exposure^-2, data = dat_tea[dat_tea$r2_tea<dat_tea$r2_caff,])) #stieger filtering
SNPcaff<-dat_tea$SNP[dat_tea$r2_tea<dat_tea$r2_caff]
summary(lm(beta.outcome~0+beta.exposure,weights = se.outcome^-2, data = dat_tea[dat_tea$SNP%in%SNPcaff,])) #stieger filtering, n.b. higher scores of these SNPs predict lower circ caffine so this is in agreement with the cis study
summary(lm(beta.outcome~0+beta.caff,weights = se.outcome^-2, data = dat_tea[dat_tea$SNP%in%SNPcaff,])) #stieger filtering, n.b. higher scores of these SNPs predict lower circ caffine so this is in agreement with the cis study
mean(dat_tea$beta.caff[dat_tea$r2_tea<dat_tea$r2_caff]^2/dat_tea$se.caff[dat_tea$r2_tea<dat_tea$r2_caff]^2)
mean(dat_tea$beta.outcome[dat_tea$r2_tea<dat_tea$r2_caff]^2/dat_tea$se.outcome[dat_tea$r2_tea<dat_tea$r2_caff]^2)
mean(dat_tea$beta.exposure[dat_tea$r2_tea<dat_tea$r2_caff]^2/dat_tea$se.exposure[dat_tea$r2_tea<dat_tea$r2_caff]^2)

mean(dat_tea$beta.exposure[dat_tea$r2_tea>dat_tea$r2_caff]^2/dat_tea$se.exposure[dat_tea$r2_tea>dat_tea$r2_caff]^2)
mean(dat_tea$beta.caff[dat_tea$r2_tea>dat_tea$r2_caff]^2/dat_tea$se.caff[dat_tea$r2_tea>dat_tea$r2_caff]^2)

#hetrgneaty in the behavoural SNPs
res<-mr(dat_tea[dat_tea$SNP%in%SNPdrink,])
res<-rbind(res,mr(dat_tea))
res$SNP<-c("Steiger coffee SNPs","Steiger coffee SNPs","Steiger coffee SNPs","Steiger coffee SNPs","Steiger coffee SNPs","no Steiger filter","no Steiger filter","no Steiger filter","no Steiger filter","no Steiger filter")

write.csv(res,"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\snp sellecion\\SuppTableS3b.csv")
#Isq(bmi$beta.outcome[bmi$SNP%in%SNPdrink]/bmi$beta.exposure[bmi$SNP%in%SNPdrink],bmi$se.outcome[bmi$SNP%in%SNPdrink]/abs(bmi$beta.exposure[bmi$SNP%in%SNPdrink]))



############################################################

#MVMR
##################

#get data above

#extract tea data from caff sig SNPs, convert to outcome, repeate above, merge thsi and the above data sets. 

id<-"ukb-b-5237"#tea ukb-b-6066, coffe ukb-b-5237
tea<-extract_instruments(id, p1 = 5e-07, clump = F) 
tea<-tea[c("SNP","pval.exposure")]

dat2<-read.csv("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\sig_caff_snps.csv")
#dat2$SNP<-dat2$SNP.x
dat2<-dat2[dat2$P<5e-07,]
dat2$pval.exposure<-dat2$P
dat2<-dat2[c("SNP","pval.exposure")]

dat2<-rbind(dat2,tea)
dat2<-clump_data(dat2,  clump_p1 = 5e-07, clump_r2 = 0.001, clump_kb = 10000,)

tea<-extract_outcome_data(dat2$SNP,id,proxies = F) 
tea<-convert_outcome_to_exposure(tea)
tea$`chr:pos`<-paste("",tea$chr.exposure,":",tea$pos.exposure,"",sep="")

dat_tea<-dat[dat$`chr:pos` %in% tea$`chr:pos`,]
dat_tea<-merge(tea,dat_tea,by="chr:pos")
dat_tea$beta.caff<-dat_tea$Effect/(2*dat_tea$eaf.exposure*(1-dat_tea$eaf.exposure)*(9876 +dat_tea$Effect^2))^0.5
dat_tea$se.caff<-1/(2*dat_tea$eaf.exposure*(1-dat_tea$eaf.exposure)*(9876 +dat_tea$Effect^2))^0.5
dat_tea$EA[dat_tea$EA=="a"]<-"A"
dat_tea$EA[dat_tea$EA=="g"]<-"G"
dat_tea$EA[dat_tea$EA=="t"]<-"T"
dat_tea$EA[dat_tea$EA=="c"]<-"C"
dat_tea$beta.caff[dat_tea$effect_allele.exposure!=dat_tea$EA]<- -1 *dat_tea$beta.caff[dat_tea$effect_allele.exposure!=dat_tea$EA]

bmi<-harmonise_data(dat_tea[colnames(tea)],extract_outcome_data(dat_tea$SNP,"ieu-b-40",proxies = F)) #t2d ebi-a-GCST006867, bmi  ieu-b-40, az ebi-a-GCST90012877
bmi<-bmi[c("SNP","beta.outcome","se.outcome")]
dat_tea<-merge(dat_tea,bmi,by="SNP")

cor(dat_tea$beta.exposure,dat_tea$beta.caff)^2
summary(lm(beta.outcome~0+beta.exposure+beta.caff,weights = se.outcome^-2, data = dat_tea))
# summary(lm(beta.outcome~beta.exposure+beta.caff,weights = se.outcome^-2, data = dat_tea)) 
# 
# MendelianRandomization::mr_mvmedian(MendelianRandomization::mr_mvinput(bx = cbind(dat_tea$beta.exposure, dat_tea$beta.caff), bxse = cbind(dat_tea$se.exposure, dat_tea$se.caff),
#                                                                        by = dat_tea$beta.outcome, byse = dat_tea$se.outcome), iterations = 1000)
# MVMRmode::mv_mrmode(Bout=dat_tea$beta.outcome, 
#                     Bexp =cbind(dat_tea$beta.exposure, dat_tea$beta.caff),
#                     SEout =dat_tea$se.outcome, 
#                     SEexp =cbind(dat_tea$se.exposure, dat_tea$se.caff), CIMin = -10, CIMax =10, CIStep=0.001  )

F.data<-format_mvmr( BXGs = dat_tea[,c("beta.exposure","beta.caff")], BYG = dat_tea$beta.outcome,
                     seBXGs = dat_tea[,c("se.exposure","se.caff")],   seBYG = dat_tea$se.outcome,  RSID = dat_tea$SNP)

sres<-strength_mvmr(r_input=F.data) # F-statistic  41.96133  5.001634
qhet_mvmr(r_input = F.data, pcor=matrix(c(1,1 ,1 ,1), ncol = 2, nrow = 2), CI=T, iterations=50) #Exposure 1       0.55611841  0.296-1.829
qhet_mvmr(r_input = F.data, pcor=matrix(c(1,.5 ,.5 ,1), ncol = 2, nrow = 2), CI=F, iterations=1000)
qhet_mvmr(r_input = F.data, pcor=matrix(c(1,0 ,0 ,1), ncol = 2, nrow = 2), CI=F, iterations=1000)
#giving the same estiamtes under differnt assumptions

# library("simex")  #load simex package
# F.data$egger_gy<-F.data$betaYG*sign(F.data$betaX1)#MR egger requires that all gene--exposure estimates are positive. flipping the sign of the g-y if it is negative
# F.data$egger_gx2<-F.data$betaX2*sign(F.data$betaX1)
# F.data$egger_gx<-abs(F.data$betaX1) #this is then making sure that the bx are positive. in MVMR we have to make this decision for one exposure, this is assuming it is the first (as with code above) 
# summary(egger_simex<-simex(egger<-lm(F.data$egger_gy~F.data$egger_gx+F.data$egger_gx2, x=TRUE,y=TRUE, weights=1/F.data$sebetaYG^2),B=1000, measurement.error = cbind(F.data$sebetaX1,F.data$sebetaX2), SIMEXvariable=c("F.data$egger_gx","F.data$egger_gx2"),fitting.method ="quad",asymptotic="FALSE") )

#############################################################################################
###########################################################################################
#########################################################################################

# control outcomes
##################
cont<-extract_instruments(c("ukb-b-4078","ukb-b-12558","ukb-b-8553"), p1 = 5e-07) ##green tea ukb-b-4078 decaf tea ukb-b-12558, decaf coffe ukb-b-8553
bmi<-harmonise_data(cont,extract_outcome_data(cont$SNP,"ieu-b-40",proxies = F)) #t2d ebi-a-GCST006867, bmi  ieu-b-40, az ebi-a-GCST90012877
mr(bmi, method_list=c("mr_ivw"))


# negative controls for pop structre
###################################

tea<-extract_instruments(c("ukb-b-5237","ukb-b-6066"), p1 = 5e-07) #tea ukb-b-6066, coffe ukb-b-5237
out<- extract_outcome_data(
  #exp_em6$SNP,
  tea$SNP, #list of snps to include
  c("ukb-d-1747_1", "ukb-d-1747_2", "ukb-d-1747_3","ukb-d-1747_4","ukb-d-1747_5","ukb-d-1747_6"),
  proxies = F, # Look for LD tags? Default is TRUE.
  rsq = 0.8, # Minimum LD rsq value (if proxies = 1). Default = 0.8.
  align_alleles = 1, # Try to align tag alleles to target alleles (if proxies = 1). 1 = yes, 0 = no. The default is 1.
  palindromes = 1, # Allow palindromic SNPs (if proxies = 1). 1 = yes, 0 = no. The default is 1.
  maf_threshold = 0.3,
  access_token = ieugwasr::check_access_token(),
  splitsize = 10000,
  proxy_splitsize = 500
)
dat2 <- harmonise_data(tea, out, action = 2)
mr_results2 <- mr(dat2, method_list=c("mr_ivw"))
write.csv(mr_results2,"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\snp sellecion\\SuppTableS1.csv")


#understanding the cuasla pathwy
################################
snp<-data.table::fread("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\snp sellecion\\tea_snps_stigerplasma_PhenoScanner_GWAS.tsv.gz")
snp<-data.table::fread("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\snp sellecion\\coffe_snps_stigerplasma_PhenoScanner_GWAS.tsv.gz")



snp<-data.table::fread("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\snp sellecion\\tea_snps_stigerdrinking_PhenoScanner_GWAS.tsv.gz")

mat<-matrix(ncol = length(unique(snp$rsid)), nrow = length(unique(snp$trait)))
colnames(mat) = unique(snp$rsid)
rownames(mat) = unique(snp$trait)
for (i in 1:nrow(snp)){
  mat[snp$trait[i],snp$snp[i]]<-T
}
View(mat)
write.csv(mat,"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\snp sellecion\\tea_snps_stigerdrinking_PhenoScanner_GWAS.csv")


snp<-data.table::fread("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\snp sellecion\\coffe_snps_stigerdrinking_PhenoScanner_GWAS.tsv.gz")
mat<-matrix(ncol = length(unique(snp$rsid)), nrow = length(unique(snp$trait)))
colnames(mat) = unique(snp$rsid)
rownames(mat) = unique(snp$trait)
for (i in 1:nrow(snp)){
  mat[snp$trait[i],snp$snp[i]]<-T
}
View(mat)
write.csv(mat,"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\snp sellecion\\coffe_snps_stigerdrinking_PhenoScanner_GWAS.csv")

snp<-data.table::fread("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\snp sellecion\\tea_snps_stigerdrinking_PhenoScanner_GWAS.tsv.gz")
snp2<-data.table::fread("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\snp sellecion\\coffe_snps_stigerdrinking_PhenoScanner_GWAS.tsv.gz")
mat<-matrix(ncol = length(unique(c(snp$rsid,snp2$rsid))), nrow = length(unique(c(snp$trait,snp2$trait))))
colnames(mat) = unique(c(snp$rsid,snp2$rsid))
rownames(mat) = unique(c(snp$trait,snp2$trait))
for (i in colnames(mat)){
  for (j in rownames(mat))
  if ((!(is.na(snp$snp[snp$snp==i])|is.na(snp$trait[snp$trait==j])))[1]){
  mat[j,i]<-T
  }
  if ((!(is.na(snp2$snp[snp$snp==i])|is.na(snp2$trait[snp$trait==j])))[1]){
    mat[j,i]<-T
  }}
View(mat)
write.csv(mat,"C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\snp sellecion\\coffe_snps_stigerdrinking_PhenoScanner_GWAS.csv")



snp<-data.table::fread("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\snp sellecion\\tea_snps_stigerplasma_PhenoScanner_GWAS.tsv.gz")
snp<-data.table::fread("C:\\Users\\nf18217\\OneDrive - University of Bristol\\Documents\\Projects\\active\\AA done or published\\caffine bmi t2d\\snp sellecion\\coffe_snps_stigerplasma_PhenoScanner_GWAS.tsv.gz")
