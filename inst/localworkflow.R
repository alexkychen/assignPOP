#Local workflow 
set.seed(124)
library(assignPOP);library(ggplot2)

# Data analysis in tutorial v1.1 (simulated genetic and morphormetric data)
## Evaluate baseline using genetic data alone
genin <- read.genpop( "inst/extdata/simGenepop.txt", pop.names=c("pop_A","pop_B","pop_C") )
genin_rd <- reduce.allele(genin, p = 0.95)

assign.MC(genin_rd, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1), 
          loci.sample="fst", iterations=30, dir="Result-folder_svm_simgenetics/", model="svm" )
assign.MC(genin_rd, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1), 
          loci.sample="fst", iterations=30, dir="Result-folder_lda_simgenetics/", model="lda" )

accuMC_svm <- accuracy.MC(dir = "Result-folder_svm_simgenetics/") ##Used this in the tutorial
accuMC_svm <- read.table("inst/extdata/AA-MC-svm-simgenetics.txt", header=T)
accuMC_lda <- accuracy.MC(dir = "Result-folder_lda_simgenetics/")
accuMC_lda <- read.table("inst/extdata/AA-MC-lda-simgenetics.txt", header=T)

accuracy.plot(accuMC_svm, pop=c("all","pop_A","pop_B","pop_C"))+ylim(0,1)+ggtitle("genetics SVM") ##Used this in tutorial
accuracy.plot(accuMC_lda, pop=c("all","pop_A","pop_B","pop_C"))+ylim(0,1)+ggtitle("genetics LDA")

check.loci(dir="Result-folder_svm_simgenetics/", top.loci = 20)

# Evaluate baseline using integrated data
comin <- compile.data(genin_rd, "inst/extdata/morphData.csv")

assign.MC(comin, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1), 
          loci.sample="fst", iterations=30, dir="Result-folder_svm_simgen-morph/", model="svm" )
assign.MC(comin, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1), 
          loci.sample="fst", iterations=30, dir="Result-folder_lda_simgen-morph/", model="lda" )

accuMCcom_svm <- accuracy.MC(dir = "Result-folder_svm_simgen-morph/") ##Used this in tutorial
accuMCcom_lda <- accuracy.MC(dir = "Result-folder_lda_simgen-morph/")

accuracy.plot(accuMCcom_svm, pop=c("all","pop_A","pop_B","pop_C"))+ylim(0,1)+ggtitle("integrated SVM") ##Used it in tutorial
accuracy.plot(accuMCcom_lda, pop=c("all","pop_A","pop_B","pop_C"))+ylim(0,1)+ggtitle("integrated LDA")

# Evaluate baseline usning morphormetric data
morphdf <- read.csv( "inst/extdata/morphData.csv", header=TRUE )
pop_label <- c( rep("pop_A", 30), rep("pop_B", 30), rep("pop_C", 30) ) 
morphdf_pop <- cbind(morphdf, pop_label)

assign.MC(morphdf_pop, train.inds=c(0.5, 0.7, 0.9), pca.method=T, iterations=30, dir="Result-folder_svm_PCAT_simorph/", model="svm" )
assign.MC(morphdf_pop, train.inds=c(0.5, 0.7, 0.9), pca.method=F, iterations=30, dir="Result-folder_svm_PCAF_simorph/", model="svm" )

accuMCmorph_svm_PT <- accuracy.MC(dir="Result-folder_svm_PCAT_simorph/") ##Used this in tutorial
accuMCmorph_svm_PF <- accuracy.MC(dir="Result-folder_svm_PCAF_simorph/") ##Used this in tutorial

accuracy.plot(accuMCmorph_svm_PT, pop=c("all","pop_A","pop_B","pop_C"))+ylim(0,1)+ggtitle("morph SVM PCAT") ##Used it in tutorial
accuracy.plot(accuMCmorph_svm_PF, pop=c("all","pop_A","pop_B","pop_C"))+ylim(0,1)+ggtitle("morph SVM PCAF") ##Used it in tutorial


#MCCV on weakly differentiated 2 pops (genetic-only)
geninfile <- read.genpop("inst/extdata/Sim2pop_n100m3000Fst0.012.gen");geninfile_rd <- reduce.allele(geninfile, p=0.95)
assign.MC(geninfile_rd, iterations = 30, model="lda", dir = "Result_MC_lda_fst_Noscaled_genetic/")

cominfile3 <- compile.data(geninfile_rd, "inst/extdata/simOtoData3.csv") 
cominfile5 <- compile.data(geninfile_rd, "inst/extdata/simOtoData5.csv")
assign.MC(cominfile3, train.loci=1, iterations = 30, model="lda", dir="Result_MC_lda_fst_Noscaled_gen-oto3/")
assign.MC(cominfile5, train.loci=1, iterations = 30, model="lda", dir="Result_MC_lda_fst_Noscaled_gen-oto5/")

otoinfile3 <- read.csv("inst/extdata/simOtoData3wPop.csv", header = T)
otoinfile5 <- read.csv("inst/extdata/simOtoData5wPop.csv", header = T)
assign.MC(otoinfile3, iterations = 30, model="lda", dir = "Result_MC_lda_Noscaled_oto3/")
assign.MC(otoinfile5, iterations = 30, model="lda", dir = "Result_MC_lda_Noscaled_oto5/")

accuRes_gen1 <- accuracy.MC(dir="Result_MC_lda_fst_Noscaled_genetic/") 
accuRes_com1 <- accuracy.MC(dir="Result_MC_lda_fst_Noscaled_gen-oto3/")
accuRes_oto1 <- accuracy.MC(dir="Result_MC_lda_Noscaled_oto3/")

accuracy.plot(accuRes_gen1) + ylim(0,1) + ggtitle("Genetic-only Unscaled (overall)")
accuracy.plot(accuRes_com1) + ylim(0,1) + ggtitle("Genetic-otolith Unscaled (overall)")
accuracy.plot(accuRes_oto1) + ylim(0,1) + ggtitle("Otolith-only Unscaled (overall)")

###Test kfold
otodf <- read.csv("inst/extdata/simOtoData5wPop.csv", header=T)
assign.kfold(geninfile_rd, model="lda", dir="Result_KF_lda_fst_Noscaled_genetic/")
assign.kfold(cominfile, model = "lda", dir="Result_KF_lda_fst_Noscaled_gen-oto5/")
assign.kfold(otodf, model="lda", dir="Result_KF_Noscaled_oto5/")

membership.plot(dir="Result_KF_lda_fst_Noscaled_genetic/", style = 1)
membership.plot(dir="Result_KF_lda_fst_Noscaled_gen-oto5/", style = 1)
membership.plot(dir="Result_KF_Noscaled_oto5/", style = 1)

#####
geninfile <- read.genpop("inst/extdata/simGenepop.txt", pop.names = c("popA","popB","popC"))
geninfile_rd <- reduce.allele(geninfile)

assign.MC(geninfile_rd, train.inds = c(0.5,0.7), train.loci=c(0.25,0.5,1), loci.sample="fst",iterations = 20, model="lda",
          scaled = F, dir="ResultTest_MC_fst_lda_NoScaled_3popGen/")
assign.MC(geninfile_rd, train.inds = c(0.5,0.7), train.loci=c(0.25,0.5,1), loci.sample="fst",iterations = 20, model="lda",
          scaled = T, dir="ResultTest_MC_fst_lda_Scaled_3popGen/")
assign.kfold(geninfile_rd, k.fold = c(3,4), train.loci = c(0.25,0.5,1), loci.sample = "fst", model = "lda",
             scaled = F, dir="ResultTest_KF_lda_Scaled_3popGen/")

accuResult1 <- accuracy.MC(dir="ResultTest_MC_fst_lda_NoScaled_3popGen/")
accuResult2 <- accuracy.MC(dir="ResultTest_MC_fst_lda_Scaled_3popGen/")
accuResult3 <- accuracy.kfold(dir="ResultTest_KF_lda_Scaled_3popGen/")

accuracy.plot(accuResult1, pop="all")+ylim(0,1)+ggtitle("3pop genetics-NoScaled (overall)")
accuracy.plot(accuResult2, pop="all")+ylim(0,1)+ggtitle("3pop genetics-Scaled (overall)")
accuracy.plot(accuResult3, pop="all")+ylim(0,1)+ggtitle("3pop genetics-NoScaled (overall)")

accuracy.plot(accuResult1, pop="popA")+ylim(0,1)+ggtitle("3pop genetics-NoScaled (pop A)")
accuracy.plot(accuResult2, pop="popA")+ylim(0,1)+ggtitle("3pop genetics-Scaled (pop A)")

accuracy.plot(accuResult1, pop = "popB")+ylim(0,1)+ggtitle("3pop genetics-NoScaled (pop B)")
accuracy.plot(accuResult2, pop = "popB")+ylim(0,1)+ggtitle("3pop genetics-Scaled (pop B)")

accuracy.plot(accuResult1, pop = "popC")+ylim(0,1)+ggtitle("3pop genetics-NoScaled (pop C)")
accuracy.plot(accuResult2, pop = "popC")+ylim(0,1)+ggtitle("3pop genetics-Scaled (pop C)")

membership.plot(dir="ResultTest_KF_lda_Scaled_3popGen/", style = 1)

#Analysis integrated dataset
cominfile <- compile.data(geninfile_rd, "inst/extdata/morphData.csv")

assign.MC(cominfile, train.inds = c(0.5,0.7), train.loci=c(0.25,0.5,1), loci.sample="fst",iterations = 20, model="lda",
          scaled = F, pca.method = "mixed", dir="ResultTest_MC_fst_lda_Noscaled_mixed_3popCom/" )
assign.MC(cominfile, train.inds = c(0.5,0.7), train.loci=c(0.25,0.5,1), loci.sample="fst",iterations = 20, model="lda",
          scaled = T, pca.method = "mixed", dir="ResultTest_MC_fst_lda_Scaled_mixed_3popCom/" )

assign.MC(cominfile, train.inds = c(0.5,0.7), train.loci=c(0.25,0.5,1), loci.sample="fst",iterations = 20, model="lda",
          scaled = F, pca.method = "independent", dir="ResultTest_MC_fst_lda_Noscaled_indep_3popCom/" )
assign.MC(cominfile, train.inds = c(0.5,0.7), train.loci=c(0.25,0.5,1), loci.sample="fst",iterations = 20, model="lda",
          scaled = T, pca.method = "independent", dir="ResultTest_MC_fst_lda_Scaled_indep_3popCom/" )

assign.MC(cominfile, train.inds = c(0.5,0.7), train.loci=c(0.25,0.5,1), loci.sample="fst",iterations = 20, model="lda",
          scaled = F, pca.method = "original", dir="ResultTest_MC_fst_lda_Noscaled_origin_3popCom/" )
assign.MC(cominfile, train.inds = c(0.5,0.7), train.loci=c(0.25,0.5,1), loci.sample="fst",iterations = 20, model="lda",
          scaled = T, pca.method = "original", dir="ResultTest_MC_fst_lda_Scaled_origin_3popCom/" )

accuResMix1 <- accuracy.MC(dir="ResultTest_MC_fst_lda_Noscaled_mixed_3popCom/")
accuResMix2 <- accuracy.MC(dir="ResultTest_MC_fst_lda_Scaled_mixed_3popCom/")

accuResIndep1 <- accuracy.MC(dir="ResultTest_MC_fst_lda_Noscaled_indep_3popCom/")
accuResIndep2 <- accuracy.MC(dir="ResultTest_MC_fst_lda_Scaled_indep_3popCom/")

accuResOrig1 <- accuracy.MC(dir="ResultTest_MC_fst_lda_Noscaled_origin_3popCom/")
accuResOrig2 <- accuracy.MC(dir="ResultTest_MC_fst_lda_Scaled_origin_3popCom/")

accuracy.plot(accuResMix1, pop="all")+ylim(0,1)+ggtitle("3pop genetic-morph Noscaled Mixed (overall)")
accuracy.plot(accuResMix2, pop="all")+ylim(0,1)+ggtitle("3pop genetic-morph Scaled Mixed (overall)")
accuracy.plot(accuResIndep1, pop="all")+ylim(0,1)+ggtitle("3pop genetic-morph Noscaled Independent (overall)")
accuracy.plot(accuResIndep2, pop="all")+ylim(0,1)+ggtitle("3pop genetic-morph Scaled Independent (overall)")
accuracy.plot(accuResOrig1, pop="all")+ylim(0,1)+ggtitle("3pop genetic-morph Noscaled Original (overall)")
accuracy.plot(accuResOrig2, pop="all")+ylim(0,1)+ggtitle("3pop genetic-morph Scaled Original (overall)")

accuracy.plot(accuResMix1, pop="popA")+ylim(0,1)+ggtitle("3pop genetic-morph Noscaled Mixed (pop A)")
accuracy.plot(accuResMix2, pop="popA")+ylim(0,1)+ggtitle("3pop genetic-morph Scaled Mixed (pop A)")
accuracy.plot(accuResIndep1, pop="popA")+ylim(0,1)+ggtitle("3pop genetic-morph Noscaled Independent (pop A)")
accuracy.plot(accuResIndep2, pop="popA")+ylim(0,1)+ggtitle("3pop genetic-morph Scaled Independent (pop A)")
accuracy.plot(accuResOrig1, pop="popA")+ylim(0,1)+ggtitle("3pop genetic-morph Noscaled Original (pop A)")
accuracy.plot(accuResOrig2, pop="popA")+ylim(0,1)+ggtitle("3pop genetic-morph Scaled Original (pop A)")

accuracy.plot(accuResMix1, pop="popB")+ylim(0,1)+ggtitle("3pop genetic-morph Noscaled Mixed (pop B)")
accuracy.plot(accuResMix2, pop="popB")+ylim(0,1)+ggtitle("3pop genetic-morph Scaled Mixed (pop B)")
accuracy.plot(accuResIndep1, pop="popB")+ylim(0,1)+ggtitle("3pop genetic-morph Noscaled Independent (pop B)")
accuracy.plot(accuResIndep2, pop="popB")+ylim(0,1)+ggtitle("3pop genetic-morph Scaled Independent (pop B)")
accuracy.plot(accuResOrig1, pop="popB")+ylim(0,1)+ggtitle("3pop genetic-morph Noscaled Original (pop B)")
accuracy.plot(accuResOrig2, pop="popB")+ylim(0,1)+ggtitle("3pop genetic-morph Scaled Original (pop B)")

accuracy.plot(accuResMix1, pop="popC")+ylim(0,1)+ggtitle("3pop genetic-morph Noscaled Mixed (pop C)")
accuracy.plot(accuResMix2, pop="popC")+ylim(0,1)+ggtitle("3pop genetic-morph Scaled Mixed (pop C)")
accuracy.plot(accuResIndep1, pop="popC")+ylim(0,1)+ggtitle("3pop genetic-morph Noscaled Independent (pop C)")
accuracy.plot(accuResIndep2, pop="popC")+ylim(0,1)+ggtitle("3pop genetic-morph Scaled Independent (pop C)")
accuracy.plot(accuResOrig1, pop="popC")+ylim(0,1)+ggtitle("3pop genetic-morph Noscaled Original (pop C)")
accuracy.plot(accuResOrig2, pop="popC")+ylim(0,1)+ggtitle("3pop genetic-morph Scaled Original (pop C)")

#MCCV on non-genetic data (two otolith elements)
otodf <- read.csv("inst/extdata/simOtoData5wPop.csv", header=T)

assign.MC(otodf, train.inds=c(0.5,0.7), iter=20, scaled = F, pca.method = T, model = "lda", 
           dir="ResultTest_MC_lda_Noscaled_PCATrue_otoData5/")
assign.MC(otodf, train.inds=c(0.5,0.7), iter=20, scaled = T, pca.method = T, model = "lda", 
           dir="ResultTest_MC_lda_Scaled_PCATrue_otoData5/")
assign.MC(otodf, train.inds=c(0.5,0.7), iter=20, scaled = F, pca.method = F, model = "lda", 
           dir="ResultTest_MC_lda_Noscaled_PCAFalse_otoData5/")
assign.MC(otodf, train.inds=c(0.5,0.7), iter=20, scaled = T, pca.method = F, model = "lda", 
           dir="ResultTest_MC_lda_Scaled_PCAFalse_otoData5/")

accuResOto5 <- accuracy.MC(dir="ResultTest_MC_lda_Noscaled_PCATrue_otoData5/")
accuResOto5_2 <- accuracy.MC(dir="ResultTest_MC_lda_Scaled_PCATrue_otoData5/")
accuResOto5_3 <- accuracy.MC(dir="ResultTest_MC_lda_Noscaled_PCAFalse_otoData5/")
accuResOto5_4 <- accuracy.MC(dir="ResultTest_MC_lda_Scaled_PCAFalse_otoData5/")

accuracy.plot(accuResOto5, pop = "all")+ylim(0,1)+ggtitle("oto5 Noscaled PCA-T (overall)")
accuracy.plot(accuResOto5_2, pop = "all")+ylim(0,1)+ggtitle("oto5 Scaled PCA-T (overall)")
accuracy.plot(accuResOto5_3, pop = "all")+ylim(0,1)+ggtitle("oto5 Noscaled PCA-F (overall)")
accuracy.plot(accuResOto5_4, pop = "all")+ylim(0,1)+ggtitle("oto5 Scaled PCA-F (overall)")

#K-folde cross validation



# PCA analysis
morphdf <- read.csv("inst/extdata/morphData.csv", header=T)
mean(morphdf$D1.2);sd(morphdf$D1.2);var(morphdf$D1.2)
mean(morphdf$D2.3);sd(morphdf$D2.3);var(morphdf$D2.3)
mean(morphdf$D3.4);sd(morphdf$D3.4);var(morphdf$D3.4)
mean(morphdf$D1.4);sd(morphdf$D1.4);var(morphdf$D1.4)
morphdff <- morphdf[,-1] 
pca_result1 <- prcomp(morphdff, scale=T, center=T);summary(pca_result1)
pca_result2 <- prcomp(morphdff, scale=F, center=T);summary(pca_result2)
pca_result3 <- prcomp(morphdff, cor=T);summary(pca_result3)
#Biplot results of pca_result2 and result3 are the same

morphdf_scale <- as.data.frame(scale(morphdf[,2:5])) 
mean(morphdf_scale$D1.2);sd(morphdf_scale$D1.2);var(morphdf_scale$D1.2)
mean(morphdf_scale$D2.3);sd(morphdf_scale$D2.3);var(morphdf_scale$D2.3)
mean(morphdf_scale$D3.4);sd(morphdf_scale$D3.4);var(morphdf_scale$D3.4)
mean(morphdf_scale$D1.4);sd(morphdf_scale$D1.4);var(morphdf_scale$D1.4)
pca_result4 <- prcomp(morphdf_scale, scale=T, center=T);summary(pca_result4)
pca_result5 <- prcomp(morphdf_scale, scale=F, center=T);summary(pca_result5)
pca_result6 <- prcomp(morphdf_scale, cor = T); summary(pca_result6)
#Above biplot results are the same as pca_result1

#Test how scale and center influence scores of PCs
#No centered and scaled => consistent results (scores are identical after multiplying the loadings)
pca_res1 <- prcomp(morphdff, scale=F, center=F);head(pca_res1$x)
head(as.matrix(morphdff) %*% pca_res1$rotation)

#Centered but not scaled => inconsistent results
pca_res2 <- prcomp(morphdff, scale=F, center=T);head(pca_res2$x)
head(as.matrix(morphdff) %*% pca_res2$rotation)

#Centered and scaled => inconsistent results
pca_res3 <- prcomp(morphdff, scale=T, center=T);head(pca_res3$x)
head(as.matrix(morphdff) %*% pca_res3$rotation)

pca_res4 <- prcomp(morphdff, cor=T); head(pca_res4$x)
head(as.matrix(morphdff) %*% pca_res4$rotation)

#Scale entire dataset first
morphdf_scale <- as.data.frame(scale(morphdff))
pca_res5 <- prcomp(morphdf_scale, scale=F, center=F);head(pca_res5$x)
head(as.matrix(morphdf_scale) %*% pca_res5$rotation)

pca_res6 <- prcomp(morphdf_scale, scale=F, center=T);head(pca_res6$x)
head(as.matrix(morphdf_scale) %*% pca_res6$rotation)

pca_res7 <- prcomp(morphdf_scale, scale=T, center=T);head(pca_res7$x)
head(as.matrix(morphdf_scale) %*% pca_res7$rotation)


