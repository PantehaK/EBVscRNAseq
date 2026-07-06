setwd("/mnt/backedup/home/pantehaK/Templates/Multi_codes/Test/Github/")
renv::init()
renv::snapshot(type = "all")
renv::dependencies(".")

setwd("/mnt/backedup/home/pantehaK/Templates/Multi_codes/Test/Github/")
renv::status()

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

