args = commandArgs(trailingOnly=TRUE)

WD <- getwd()
if (!is.null(WD)) setwd(WD)
################# Load master file #################
master_file <- read.delim(args[1], sep = "\t")
####################################################

############ Segregate Combinations ################
args_pass <- args[2]
combinations <- unlist(strsplit(args_pass, split = "+", fixed = T))
####################################################

############# Create Combination files #################
for (i in combinations) {
  groups <- strsplit(i, split = "_", fixed = T) 
  print(groups)
  ctrl <- as.data.frame(groups)[1,]
  trtd <- as.data.frame(groups)[2,]
  comb_split <- as.data.frame(rbind(master_file[master_file$condition == ctrl,], master_file[master_file$condition == trtd,]))
  write.table(comb_split, file = file.path("2_Combinations", paste(ctrl, trtd, sep = "_")), sep = "\t", row.names = F, quote = F)
}
######################################################
