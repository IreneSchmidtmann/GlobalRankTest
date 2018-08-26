library(knitr)

knit2pdf("SchmidtmannKonstantinidesBinderBiomJ-R3.rnw", output = "SchmidtmannKonstantinidesBinderBiomJ-R3.tex")
purl("SchmidtmannKonstantinidesBinderBiomJ-R3.rnw", documentation = 1, output = "GlobalRankTest.R")
sessionInfo()