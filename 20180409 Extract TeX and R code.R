library(knitr)

knit("SchmidtmannKonstantinidesBinderBiomJ-R1.rnw")
purl("SchmidtmannKonstantinidesBinderBiomJ-R1.rnw", documentation = 0, output = "SchmidtmannKonstantinidesBinderBiomJ-R10.R")
purl("SchmidtmannKonstantinidesBinderBiomJ-R1.rnw", documentation = 1, output = "SchmidtmannKonstantinidesBinderBiomJ-R11.R")
purl("SchmidtmannKonstantinidesBinderBiomJ-R1.rnw", documentation = 2, output = "SchmidtmannKonstantinidesBinderBiomJ-R12.R")

