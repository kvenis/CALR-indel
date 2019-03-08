getPeaks <- function(mat){
  temp_reads <- mat[mat$V1>=100,]
  temp_reads$V3 <- round(temp_reads$V1*100/sum(temp_reads$V1), digits=2)
  lengths <- head(temp_reads[order(temp_reads$V1, decreasing=TRUE),],)
  #  print ("possible CALR lengths:")
  lengths$V4 <- lengths$V2-240
  return(lengths[lengths$V2 > (240+2) 
                 | lengths$V2 < (240-2),])
  
}



R129_bc065_s5xlA_reads <- read.delim(file="Bioinformatics_Projects/CALR_indels/out_dir/R129_bc065_s5XLA_calr_readLenths.txt",
                                    sep="", header=FALSE, skip=1)
plot(R129_bc065_s5xlA_reads$V2,R129_bc065_s5xlA_reads$V1,type="l",xlab="read length",ylab="occurences",col="blue", main="R129_bc065-s5xlA")

getPeaks(R129_bc065_s5xlA_reads)

R40_bc063_s5xlB_reads <- read.delim(file="/Test Development/Bioinformatics_Projects/CALR_indels/out_dir/R40_MPCP_BC063/R40_bc063_s5xlb_calr_readLengths.txt",
                                    sep="", header=FALSE, skip=1)
plot(R40_bc063_s5xlB_reads$V2,R40_bc063_s5xlB_reads$V1,type="l",xlab="read length",ylab="occurences",col="blue", main="R40_bc063_s5xlB")

getPeaks(R40_bc063_s5xlB_reads)
data <- getPeaks(R40_bc063_s5xlB_reads)
data
lengths <- head(R40_bc063_s5xlB_reads[order(R40_bc063_s5xlB_reads$V1, decreasing=TRUE),],)
lengths
write.txt(data, file = "/Bioinformatics_Projects/CALR_indels/out_dir/R40_MPCP_BC063/R40_bc063_s5xlB.txt")
names(data) <- c("No_of_reads","read_length","percentage_of_reads","indel_length")
data

##R209-bc041-S5xlA
R209_bc041_s5xlA_reads <- read.delim(file="/Bioinformatics_Projects/CALR_indels/out_dir/R209-S5XLA/R209-bc041-S5xla-calr-exon8_readLength.txt",
                                    sep="", header=FALSE, skip=1)
plot(R209_bc041_s5xlA_reads$V2,R209_bc041_s5xlA_reads$V1,type="l",xlab="read length",ylab="occurences",col="blue", main="R209_bc041_s5xlA")

getPeaks(R209_bc041_s5xla_reads)
data <- getPeaks(R209_bc041_s5xlA_reads)
data

lengths <- head(R209_bc041_s5xlA_reads[order(R209_bc041_s5xlA_reads$V1, decreasing=TRUE),],)
lengths
names(data) <- c("No_of_reads","read_length","percentage_of_reads","indel_length")
data
write.table(data, file = "/Bioinformatics_Projects/CALR_indels/out_dir/R209-S5XLA/R209_bc041_s5xlA.txt", sep = "\t", row.names = FALSE)
