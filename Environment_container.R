
setwd("/mnt/backedup/home/pantehaK/Templates/EBV_MS_manuscripts/")
renv::init()
renv::snapshot(type = "all")
renv::dependencies(".")

setwd("/mnt/backedup/home/pantehaK/Templates/EBV_MS_manuscripts/")
renv::status()

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")


