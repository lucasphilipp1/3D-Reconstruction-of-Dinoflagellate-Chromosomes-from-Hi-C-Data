rm(list=ls())
install.packages("remotes")

#this code prepares .hic files for CSynth
#.hic data files are formatted strangely to allow for dynamic zoom
#with .hic files, only data in the field of view at an appropriate resolution is loaded in real time
#download .hic file for Breviolum minutum here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4658994

remotes::install_github("aidenlab/straw/R")
nchr<-100

str_array1<-array(,nchr)
str_array2<-array(,nchr)

.hic_to_.txt <- function(arg_1, arg_2) {
  hic.data.frame <- strawr::straw("NONE", "/Users/lucasphilipp/Desktop/Research/GitHub/straw/GSM4658994_L142+L533+L534.inter.hic", arg_1, arg_1, "BP", 5000)
  hic.data.frame[hic.data.frame==0] <- NA #remove monomer at 0bp. for some reason CSynth does not recognize the contact probability matrix when there is a monomer at 0bp
  hic.data.frame<-hic.data.frame[complete.cases(hic.data.frame),]
  temp<-paste("/Users/lucasphilipp/Downloads/hic contact probabilities/", arg_2, sep="")
  write.table(hic.data.frame, paste(temp, ".txt", sep=""), append = FALSE, sep = "\t", eol = "\n", dec = ".", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

for (x in 1:nchr) {
  str_array1[x]<-paste("HIC_SCAFFOLD_", toString(x), sep="")
  str_array2[x]<-paste("bminutum_pseudochromosome_", toString(x), sep="")
  .hic_to_.txt(str_array1[x],str_array2[x])
}

#this code prepares HiC .txt files (output from cooler dump which converts .mcool -> .txt ) for CSynth
#download .mcool files for Symbiodinium microadriaticum here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152150
#in a new terminal window:
#pip install cooler
#cooler ls GSE152150_Dino-HiC-cD1plus-R1.smic1.0.mapq_30.1000.mcool
#cooler dump -t pixels -r chr1 -r2 chr1 -b --no-balance -o symbiodinium_microadriaticum_coccoid.txt --join GSE152150_Dino-HiC-cD1plus-R1.smic1.0.mapq_30.1000.mcool::/resolutions/1000
#cooler dump -t pixels -r chr1 -r2 chr1 -b --no-balance -o symbiodinium_microadriaticum_mastigote.txt --join GSE152150_Dino-HiC-mD1plus-R1.smic1.0.mapq_30.1000.mcool::/resolutions/1000

library(dplyr)

chr_from_cooldump <- function(x,lifestage) {
  lifestage_suffix <- paste(lifestage, ".txt", sep="")
  microadriaticum<-read.table(paste('/Users/lucasphilipp/Downloads/hic contact probabilities microadriaticum/symbiodinium_microadriaticum_',lifestage_suffix, sep=""), sep = "\t", header=FALSE)
  temp <- paste("chr", x, sep="")
  temp_symbol <- paste(paste("^", temp, sep=""), "$", sep="") 
  chr <- microadriaticum %>% filter(grepl(temp_symbol, V1))
  #^ symbol denotes the start of the string
  #$ symbol denotes the end of the string
  chr <- na.omit(chr)
  chr$V1 <- NULL
  chr <- filter(chr, V2 > 0, V3 > 0) #remove monomer at 0bp. for some reason UCSF Chimera cannot be used with a 0bp start
  filename<-paste("/Users/lucasphilipp/Downloads/hic contact probabilities microadriaticum/", paste(paste(lifestage,"_", sep=""),temp, sep=""), sep="")
  write.table(chr, paste(filename, ".txt", sep=""), append = FALSE, sep = "\t", eol = "\n", dec = ".", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

for (x in 1:94) {
  chr_from_cooldump(x,"coccoid")
}

for (x in 1:94) {
  chr_from_cooldump(x,"mastigote")
}

lifestage_suffix <- paste("coccoid", ".txt", sep="")
microadriaticum<-read.table(paste('/Users/lucasphilipp/Downloads/hic contact probabilities microadriaticum/symbiodinium_microadriaticum_',lifestage_suffix, sep=""), sep = "\t", header=FALSE)
microadriaticum<-microadriaticum[,-1]
microadriaticum<-microadriaticum[,-1]
microadriaticum<-microadriaticum[,-1]

temp <- paste("chr", 1, sep="")
temp_symbol <- paste(paste("^", temp, sep=""), "$", sep="") 
chr <- microadriaticum %>% filter(grepl(temp_symbol, V4))
#^ symbol denotes the start of the string
#$ symbol denotes the end of the string
chr <- na.omit(chr)
chr$V1 <- NULL
chr <- filter(chr, V2 > 0, V3 > 0) #remove monomer at 0bp. for some reason UCSF Chimera cannot be used with a 0bp start
filename<-paste("/Users/lucasphilipp/Downloads/hic contact probabilities microadriaticum/", paste(paste(lifestage,"_", sep=""),temp, sep=""), sep="")
write.table(chr, paste(filename, ".txt", sep=""), append = FALSE, sep = "\t", eol = "\n", dec = ".", row.names = FALSE, col.names = FALSE, quote = FALSE)

