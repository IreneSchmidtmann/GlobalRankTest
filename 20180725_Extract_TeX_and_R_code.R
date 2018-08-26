library(knitr)

knit2pdf("SchmidtmannKonstantinidesBinderBiomJ-R2.rnw")
purl("SchmidtmannKonstantinidesBinderBiomJ-R2.rnw", documentation = 1, output = "GlobalRankTest.R")
sessionInfo()

