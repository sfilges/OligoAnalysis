path = getwd()
library(stringr)
library(tis)

#file.names <- c("Phusion-80ng-4_S12_m.bamPositionComposition.cons10.txt")
file.names <- list.files(pattern = "cons10.txt")

df = data.frame()

#---------------------------------// Reference sequences //---------------------------------
primer_1 ="GTGGTGAGGCTCCCCTTT"
primer_2 = "CAAAGCTGTTCCGTCCCAGT"

forward_insert = "ATACAGAATATCTGTTCGCACTCCGAGTGCGGCTTGCGGAGATTCTCTTCCTCTGTGCGCCG"
reverse_insert = "CTTGCGGAGATTCTCTTCCTCTGTGCGCCGATACAGAATATCTGTTCGCACTCCGAGTGCGG"
#-------------------------------------------------------------------------------------------


for(i in 1:length(file.names)){
  out.file<-""
  file <- read.csv(file.names[i],header=TRUE, sep="\t", stringsAsFactors=FALSE)
  colnames(file) <- c("position","barcode", "counts","sequence")
  
  file <- file[file$counts >= 10,]
  file <- file[file$barcode != "",]
  
  total_reads <- length(file$position)
  
  perfect_reads_fwd <- length(grep("ATACAGAATATCTGTTCGCACTCCGAGTGCGGCTTGCGGAGATTCTCTTCCTCTGTGCGCCG",file$sequence))/total_reads
  perfect_reads_rev <- length(grep("CTTGCGGAGATTCTCTTCCTCTGTGCGCCGATACAGAATATCTGTTCGCACTCCGAGTGCGG",file$sequence))/total_reads
  
  file$n_D <- str_count(file$sequence, "D")
  nol_del <- nrow(file[file$n_D == 0,])/total_reads
  adin_del <- nrow(file[file$n_D == 1,])/total_reads
  dwa_del <- nrow(file[file$n_D == 2,])/total_reads
  tri_del <- nrow(file[file$n_D == 3,])/total_reads
  chetire_del <- nrow(file[file$n_D == 4,])/total_reads
  pyat_del <- nrow(file[file$n_D == 5,])/total_reads
  schest_del <- nrow(file[file$n_D == 6,])/total_reads
  
  row_new <- as.data.frame(t(c(file.names[i],total_reads,nol_del,adin_del,dwa_del,tri_del,chetire_del,pyat_del,schest_del)))
  
  df <- rbind(df,row_new)
}

errors <- read.csv("180417_tp1_synth_error.csv",header=TRUE, sep="\t")
errors$mean_idt_desalted <- RowMeans(errors[,2:10])
errors$mean_idt_page <- RowMeans(errors[,11:19])

plot(errors[,2:5])
abline(0,1)


cors <- cor(errors[,2:11])