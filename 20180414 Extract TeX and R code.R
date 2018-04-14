library(knitr)

knit("SchmidtmannKonstantinidesBinderBiomJ-R1.rnw")
purl("SchmidtmannKonstantinidesBinderBiomJ-R1.rnw",
     documentation = 1, output = "SchmidtmannKonstantinidesBinderBiomJ-R1.R")
sessionInfo()