library(UpSetR)
args = commandArgs(trailingOnly=TRUE)
data <- read.csv(args[1], header = T, sep = "\t")
jpeg(args[2], width = 3000, height = 1000)
upset(data, sets = c("mapq.60","multimapped","mapped_fraction","mapq_fraction","insert_not_fount","in_reference_sample","sc_length","read_not_found","sc_not_found","flank_size"), mb.ratio = c(0.55, 0.45), order.by = "freq", text.scale=2, nintersects=50)
dev.off()

