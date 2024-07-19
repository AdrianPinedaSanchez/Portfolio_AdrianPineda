
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DECIPHER")

library("Biostrings")
library("DECIPHER")
library("seqinr")

#Lectura
omiFrancia <- read.fasta("Omicron_Francia.fasta")
delIndia <- read.fasta("Delta_India.fasta")
omiIndia <- read.fasta("Omicron_India.fasta")
alphaUSA <- read.fasta("Alpha_USA.fasta")
omiUSA <- read.fasta("Omicron_USA.fasta")
omiOCEA <- read.fasta("Omicron_OCEA.fasta")
alphaSDA <- read.fasta("Alpha_SDA.fasta")
delSDA <- read.fasta("Delta_SDA.fasta")

#Longitud del genoma
omiFranciaLength <- length(omiFrancia[[1]])
delIndiaLength <- length(delIndia[[1]])
omiIndiaLength <- length(omiIndia[[1]])
alphaUSALength <- length(alphaUSA[[1]])
omiUSALength <- length(omiUSA[[1]])
omiOCEALength <- length(omiOCEA[[1]])
alphaSDALength <- length(alphaSDA[[1]])
delSDALength <- length(delSDA[[1]])

genomasLongitudes <- c(omiFranciaLength, delIndiaLength, omiIndiaLength, alphaUSALength, omiUSALength, omiOCEALength, alphaSDALength, delSDALength)

#Contenido de GC
omiFranciaGC <- GC(omiFrancia[[1]])
delIndiaGC <- GC(delIndia[[1]])
omiIndiaGC <- GC(omiIndia[[1]])
alphaUSAGC <- GC(alphaUSA[[1]])
omiUSAGC <- GC(omiUSA[[1]])
omiOCEAGC <- GC(omiOCEA[[1]])
alphaSDAGC <- GC(alphaSDA[[1]])
delSDAGC <- GC(delSDA[[1]])

genomasGC <-c(omiFranciaGC, delIndiaGC, omiIndiaGC, alphaUSAGC, omiUSAGC, omiOCEAGC, alphaSDAGC, delSDAGC)

aPorcentajeLista = c()
cPorcentajeLista = c()
tPorcentajeLista = c()
gPorcentajeLista = c()
nombres = c()

composicion <- function(genomaVirus, nombre){
  
  genoma = readDNAStringSet(genomaVirus)
  genoma
  
  ancho <- width(genoma)
  frecuencia <- alphabetFrequency(genoma, baseOnly = TRUE)
  
  aFrec <- frecuencia[1]
  cFrec <- frecuencia[2]
  gFrec <- frecuencia[3]
  tFrec <- frecuencia[4]
  
  cat("A: ", aFrec,
      "\nC: ", cFrec,
      "\nG: ", gFrec,
      "\nT: ", tFrec
  )
  
  aPorcentaje <- (aFrec/ancho) * 100
  cPorcentaje <- (cFrec/ancho) * 100
  gPorcentaje <- (gFrec/ancho) * 100
  tPorcentaje <- (tFrec/ancho) * 100
  
  cat("\n% de A: ", aPorcentaje,
      "\n% de C: ", cPorcentaje,
      "\n% de G: ", gPorcentaje,
      "\n% de T: ", tPorcentaje)
  
  
  aPorcentajeLista <<- c(aPorcentajeLista, aPorcentaje)
  cPorcentajeLista <<- c(cPorcentajeLista, cPorcentaje)
  tPorcentajeLista <<- c(tPorcentajeLista, tPorcentaje)
  gPorcentajeLista <<- c(gPorcentajeLista, gPorcentaje)
  nombres <<- c(nombres, nombre)
  
}

composicion("Omicron_Francia.fasta", "Omicron-Francia")
composicion("Delta_India.fasta", "Delta-India")
composicion("Omicron_India.fasta", "Omicron-India")
composicion("Alpha_USA.fasta", "Alpha-USA")
composicion("Omicron_USA.fasta", "Omicron-USA")
composicion("Omicron_OCEA.fasta", "Omicron-Oceanía")
composicion("Alpha_SDA.fasta", "Alpha-Sudáfrica")
composicion("Delta_SDA.fasta", "Delta-Sudáfrica")



par(mar=c(7,4,4,4))
barplot(aPorcentajeLista, 
        horiz = FALSE, 
        names.arg = nombres, 
        axes = TRUE, 
        xlab="",
        ylab="Porcentaje", 
        las = 2, 
        cex.names = 0.8, 
        cex.axis = 0.8, 
        ylim = c(29.8, 29.95), 
        xpd = FALSE, 
        col="#ff4040",
        main = "Porcentaje de Adenina", 
        border="black")

par(mar=c(7,4,4,4))
barplot(cPorcentajeLista, 
        horiz = FALSE, 
        names.arg = nombres, 
        axes = TRUE, 
        xlab="",
        ylab="Porcentaje", 
        las = 2, 
        cex.names = 0.8, 
        cex.axis = 0.8, 
        ylim = c(18.27, 18.38), 
        xpd = FALSE, 
        col="#9cff45",
        main = "Porcentaje de Citocina", 
        border="black")

par(mar=c(7,4,4,4))
barplot(tPorcentajeLista, 
        horiz = FALSE, 
        names.arg = nombres, 
        axes = TRUE, 
        xlab="",
        ylab="Porcentaje", 
        las = 2, 
        cex.names = 0.8, 
        cex.axis = 0.8, 
        ylim = c(32.07, 32.23), 
        xpd = FALSE, 
        col="#ffd045",
        main = "Porcentaje de Timina", 
        border="black")

par(mar=c(7,4,4,4))
barplot(gPorcentajeLista, 
        horiz = FALSE, 
        names.arg = nombres, 
        axes = TRUE, 
        xlab="",
        ylab="Porcentaje", 
        las = 2, 
        cex.names = 0.8, 
        cex.axis = 0.8, 
        ylim = c(19.58, 19.65), 
        xpd = FALSE, 
        col="#3eb8f0",
        main = "Porcentaje de Guanina", 
        border="black")

par(mar=c(7,4,4,4))
barplot(genomasLongitudes, 
        horiz = FALSE, 
        names.arg = nombres, 
        axes = TRUE, 
        xlab="",
        ylab="Longitud", 
        las = 2, 
        cex.names = 0.8, 
        cex.axis = 0.8, 
        ylim = c(29500, 30000), 
        xpd = FALSE, 
        col="#ff4040",
        main = "Longitudes", 
        border="black")

par(mar=c(7,4,4,4))
barplot(genomasGC, 
        horiz = FALSE, 
        names.arg = nombres, 
        axes = TRUE, 
        xlab="",
        ylab="GC", 
        las = 2, 
        cex.names = 0.8, 
        cex.axis = 0.8, 
        ylim = c(0.378, 0.381), 
        xpd = FALSE, 
        col="#0096ff",
        main = "Contenido GC", 
        border="black")

corona_virus <- c("ON287427","OK189630","MZ336026","ON134749","ON188698","OM737996","OM765457","OM773365")
virus_sequences <- read.GenBank(corona_virus)
virus_sequences

write.dna(virus_sequences, file="coronavirus_secuencias.fasta", format="fasta")

virus_seq_not_align <- readDNAStringSet("coronavirus_secuencias.fasta", format="fasta")

virus_seq_align <- AlignSeqs(virus_seq_not_align)
BrowseSeqs(virus_seq_align)

writeXStringSet(virus_seq_align, file="coronavirus_sec_align.fasta")

virus_aligned <- read.alignment("coronavirus_sec_align.fasta", format="fasta")

matriz_distancia <- dist.alignment(virus_aligned, matrix = "similarity")
matriz_distancia_frame <- as.data.frame(as.matrix(matriz_distancia))
tabla_grises <- table.paint(matriz_distancia_frame, cleg=0, clabel.row=.5, clabel.col=.5) + scale_color_viridis()

virus_tree <- nj(matriz_distancia)

virus_tree <- ladderize(virus_tree)
plot(virus_tree)

