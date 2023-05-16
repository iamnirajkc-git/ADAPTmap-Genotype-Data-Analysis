install.packages("tidyverse", dependencies = TRUE) #for data modification, visualization
library(tidyverse)


#set working directory
setwd("C:\\Users\\Eigenaar\\Desktop\\ADAPTmap_genotypeTOP_20160222_full")


#run PLINK: bped to ped
system("plink --bfile ADAPTmap_genotypeTOP_20160222_full  --recode --out ADAPTmap_genotypeTOP_20160222_full")

# goat and cow has equal number of autosomes.
system("plink --bfile ADAPTmap_genotypeTOP_20160222_full --cow  --recode --out ADAPTmap_genotypeTOP_20160222_full")

#run PLINK: bped to vcf
#Converting from bped to VCF allows genotype data to be easily shared between 
#different software tools and platforms that support VC
system("plink --bfile ADAPTmap_genotypeTOP_20160222_full --cow  --recode vcf --out ADAPTmap_genotypeTOP_20160222_full")

#run PLINK: ped to vcf
system("plink --file ADAPTmap_genotypeTOP_20160222_full --cow  --recode vcf --out ADAPTmap_genotypeTOP_20160222_full")

#Quality control analysis
# Missingness per SNP: 0.1 --geno
# Missingness per indivisual: 0.1 --mind
# Minor allele frequency: 0.05 --maf
# Run PLINK QC
system("plink --bfile ADAPTmap_genotypeTOP_20160222_full --cow --autosome --geno 0.1 --mind 0.1 --maf 0.05 --nonfounders --allow-no-sex --recode --out ADAPTmap_TOP")

######################################
# Principal Component Analysis - PCA #
######################################

## Genetic distances between individuals
system("plink --cow --allow-no-sex --nonfounders --file ADAPTmap_TOP --distance-matrix --out dataForPCA")

## Load data
dist_populations<-read.table("dataForPCA.mdist",header=F)
### Extract breed names
fam <- data.frame(famids=read.table("dataForPCA.mdist.id")[,1])
### Extract individual names 
famInd <- data.frame(IID=read.table("dataForPCA.mdist.id")[,2])

## Perform PCA using the cmdscale function 
# Time intensive step - takes a few minutes with the 4.5K animals
# eig	indicates whether eigenvalues should be returned, k is dimension of space
mds_populations <- cmdscale(dist_populations,eig=T,5)

## Extract the eigen vectors
eigenvec_populations <- cbind(fam,famInd,mds_populations$points)

## Proportion of variation captured by each eigen vector
eigen_percent <- round(((mds_populations$eig)/sum(mds_populations$eig))*100,2)

# PCA plot
ggplot(data = eigenvec_populations) +
  geom_point(mapping = aes(x = `1`, y = `2`,color = famids), show.legend = FALSE ) + 
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = "PCA of wordwide goat populations",
       x = paste0("Principal component 1 (",eigen_percent[1]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[2]," %)")) + 
  theme_minimal()



