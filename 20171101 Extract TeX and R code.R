library(knitr)

knit("201700901_GMDS_GlobalRankNoninferiority.rnw")
purl("201700901_GMDS_GlobalRankNoninferiority.rnw", documentation = 0, output="test0.R")
purl("201700901_GMDS_GlobalRankNoninferiority.rnw", documentation = 1, output="test1.R")
purl("201700901_GMDS_GlobalRankNoninferiority.rnw", documentation = 2, output="test2.R")

