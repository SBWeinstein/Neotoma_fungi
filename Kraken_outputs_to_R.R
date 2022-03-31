#Some really ugly R code to get stuff out of mpa formatted kraken reports.

setwd("<folder of kraken2 outputs>")  #contains kraken2 .report files and the slurm output file

#packages
library("tidyr")

#get amounts classified from slurm output (single file from whole run)
slurm <- read.table("slurm-4195701.out", 
                      sep ="\t", header = FALSE, dec =".")

dim(slurm) #each sample gets 4 lines in output (loading, total seqs, classified, unclasssified)
slurm1<-data.frame(slurm)
col1<-c(rep(1:4, 66) )
slurm1$V3<-col1
slurm2 <- slurm1[ which(slurm1$V3 != 1),]

#a way to get some numbers 
m <- gregexpr('[0-9]+\\s+s',slurm2$V1)
b<- regmatches(slurm2$V1 ,m)
slurm2$V4<-unlist(b)
slurm2$V5<-as.numeric(sub('\\s+s','',slurm2$V4))
slurm2$V6<- c(rep(c("all_reads", "classified", "unclassified"), 66))

fileNames <- Sys.glob("*.report")
slurm2$V7<-rep(fileNames, each = 3)
df_reads<-data.frame(sample = slurm2$V7, read_type = slurm2$V6, reads = slurm2$V5)

#use mpa converted outputs from krakentools; because...tab delimited numbers

metafiles <- Sys.glob("*.report") 

df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(df) <- c("sample", "read_type", "reads")

for (fileName in metafiles) {
  data2 <- read.delim(fileName, 
                        sep ="\t", header = FALSE )
  Fungi_Rows <- data2[grep("k__Fungi", data2$V1), ]
  Neo_Rows <- data2[grep("c__Neocallimastigomycetes", data2$V1), ]
  Bac_Rows <- data2[grep("k__Bacteria", data2$V1), ]
  plant_Rows <- data2[grep("k__Viridiplantae", data2$V1), ]
  
 #create new dataframe
 sample1 <- data.frame(matrix(ncol = 3, nrow = 4))
 
 #provide column names
 colnames(sample1) <- c("sample", "read_type", "reads")
  sample1$sample = c(rep(fileName, 4))
  sample1$read_type= c("Fungi", "Neo", "Bac", "Plant")
  sample1$reads = c(Fungi_Rows[1,2], Neo_Rows[1,2], Bac_Rows[1,2],  plant_Rows[1,2])
  
  df <- rbind(df,sample1)
  
}

#combine df and df_reads into one dataframe
df2<-rbind(df_reads, df)
df3<-spread(df2, key=read_type, value=reads)

#calculate proportion of reads classified, etc
df3$per_class<-df3$classified/df3$all_reads
mean(df3$per_class)
sd(df3$per_class)

#prop classified reads that are bacteria
df3$Bac[is.na(df3$Bac)] <- 0
df3$pCBac<-df3$Bac/df3$classified
mean(df3$pCBac)
sd(df3$pCBac)

#Fungi
df3$Fungi[is.na(df3$Fungi)] <- 0
df3$pCFun<-df3$Fungi/df3$classified
mean(df3$pCFun )
sd(df3$pCFun)

#Plants
df3$Plant[is.na(df3$Plant)] <- 0
df3$pCplant<-df3$Plant/df3$classified
mean(df3$pCplant)
sd(df3$pCplant)

#Neo reads
df3$Neo[is.na(df3$Neo)] <- 0
mean(df3$Neo)
sd(df3$Neo)

write.csv(df3, "kraken_bigDB_Conf0.05_outputs_30Mar22.csv")