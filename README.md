# ADAPTmap-Genotype-Data-Analysis
This repository contains code for analyzing the ADAPTmap genotype data. The analysis includes data quality , visualization, and performing Principal Component Analysis 

# Steps involved and R codes.
# Installation
To run the analysis, make sure you have the tidyverse package installed. If not, you can install it by running the following command:
install.packages("tidyverse", dependencies = TRUE)

# Data Preparation
Set the working directory to the location where the ADAPTmap genotype data is stored. For example:
setwd("C:\\Users\\Eigenaar\\Desktop\\ADAPTmap_genotypeTOP_20160222_full")

# Convert the genotype data from .bped to .ped format using PLINK:
system("plink --bfile ADAPTmap_genotypeTOP_20160222_full --recode --out ADAPTmap_genotypeTOP_20160222_full")

# Convert the genotype data from .bped to VCF format using PLINK:
system("plink --bfile ADAPTmap_genotypeTOP_20160222_full --recode vcf --out ADAPTmap_genotypeTOP_20160222_full")

# Perform quality control analysis using PLINK:
system("plink --bfile ADAPTmap_genotypeTOP_20160222_full --autosome --geno 0.1 --mind 0.1 --maf 0.05 --nonfounders --allow-no-sex --recode --out ADAPTmap_TOP")

# Principal Component Analysis (PCA)
Calculate genetic distances between individuals using PLINK:
system("plink --cow --allow-no-sex --nonfounders --file ADAPTmap_TOP --distance-matrix --out dataForPCA")

# Load the data and extract breed and individual names:
dist_populations <- read.table("dataForPCA.mdist", header = FALSE)
fam <- data.frame(famids = read.table("dataForPCA.mdist.id")[, 1])
famInd <- data.frame(IID = read.table("dataForPCA.mdist.id")[, 2])

# Perform PCA using the cmdscale function:
mds_populations <- cmdscale(dist_populations, eig = TRUE, k = 5)

# Extract the eigen vectors and calculate the proportion of variation captured by each eigen vector:
eigenvec_populations <- cbind(fam, famInd, mds_populations$points)
eigen_percent <- round((mds_populations$eig / sum(mds_populations$eig)) * 100, 2)

# Visualize the PCA results using a scatter plot:
ggplot(data = eigenvec_populations) +
  geom_point(mapping = aes(x = `1`, y = `2`, color = famids), show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  labs(
    title = "PCA of worldwide goat populations",
    x = paste0("Principal component 1 (", eigen_percent[1], " %)"),
    y = paste0("Principal component 2 (", eigen_percent[2], " %)")
  ) +
  theme_minimal()

# Results
The PCA plot shows the distribution of worldwide goat populations based on genetic distances.
