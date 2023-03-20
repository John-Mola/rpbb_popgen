# install.packages(c("fields","RColorBrewer","mapplots"))
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.16")
# 
# BiocManager::install("LEA")
# # source("http://bioconductor.org/biocLite.R")
# # biocLite("LEA")
# 
# source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
# source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")
# 
# input.file = "http://membres-timc.imag.fr/Olivier.Francois/secondary_contact.str"
# struct2geno(file = input.file, TESS = TRUE, diploid = TRUE, FORMAT = 2,
#             extra.row = 0, extra.col = 0, output = "secondary_contact.geno")
# library(LEA)
# obj.snmf = snmf("secondary_contact.geno", K = 3, alpha = 100, project = "new")
# qmatrix = Q(obj.snmf, K = 3)
# 
# barplot(t(qmatrix), col = c("orange","violet","lightgreen"), border = NA, space = 0,
#         xlab = "Individuals", ylab = "Admixture coefficients")
# 
# coord = read.table("coordinates.coord")
# pop = rep(1:60, each = 10)
# 
# K = 3
# Npop = length(unique(pop))
# qpop = matrix(NA, ncol = K, nrow = Npop)
# coord.pop = matrix(NA, ncol = 2, nrow = Npop)
# for (i in unique(pop)){
#   qpop[i,] = apply(qmatrix[pop == i,], 2, mean)
#   coord.pop[i,] = apply(coord[pop == i,], 2, mean)}
# 
# 
# library(mapplots)
# plot(coord, xlab = "Longitude", ylab = "Latitude", type = "n")
# map(add = T, col = "grey90", fill = TRUE)
# for (i in 1:Npop){
#   add.pie(z = qpop[i,], x = coord.pop[i,1], y = coord.pop[i,2], labels = "",
#           col = c("orange","violet","lightgreen"))}
# 
# pop = scan("mypop.txt")
# 
# ####
# 
# # example use: 
# # data(nancycats)
# # genind2structure(nancycats, file="nancy_structure.txt", pops=TRUE)
# 
# genind2structure <- function(obj, file="", pops=FALSE){
#   if(!"genind" %in% class(obj)){
#     warning("Function was designed for genind objects.")
#   }
#   
#   # get the max ploidy of the dataset
#   pl <- max(obj@ploidy)
#   # get the number of individuals
#   S <- adegenet::nInd(obj)
#   # column of individual names to write; set up data.frame
#   tab <- data.frame(ind=rep(indNames(obj), each=pl))
#   # column of pop ids to write
#   if(pops){
#     popnums <- 1:adegenet::nPop(obj)
#     names(popnums) <- as.character(unique(adegenet::pop(obj)))
#     popcol <- rep(popnums[as.character(adegenet::pop(obj))], each=pl)
#     tab <- cbind(tab, data.frame(pop=popcol))
#   }
#   loci <- adegenet::locNames(obj) 
#   # add columns for genotypes
#   tab <- cbind(tab, matrix(-9, nrow=dim(tab)[1], ncol=adegenet::nLoc(obj),
#                            dimnames=list(NULL,loci)))
#   
#   # begin going through loci
#   for(L in loci){
#     thesegen <- obj@tab[,grep(paste("^", L, "\\.", sep=""), 
#                               dimnames(obj@tab)[[2]]), 
#                         drop = FALSE] # genotypes by locus
#     al <- 1:dim(thesegen)[2] # numbered alleles
#     for(s in 1:S){
#       if(all(!is.na(thesegen[s,]))){
#         tabrows <- (1:dim(tab)[1])[tab[[1]] == indNames(obj)[s]] # index of rows in output to write to
#         tabrows <- tabrows[1:sum(thesegen[s,])] # subset if this is lower ploidy than max ploidy
#         tab[tabrows,L] <- rep(al, times = thesegen[s,])
#       }
#     }
#   }
#   
#   # export table
#   write.table(tab, file=file, sep="\t", quote=FALSE, row.names=FALSE)
# }
# 
# 
# data(nancycats)
# genind2structure(nancycats, file="nancy_structure.txt", pops=TRUE)
# 
# 
# struct2geno(file = "./nancy_structure.txt", TESS = TRUE, diploid = TRUE, FORMAT = 2,
#             extra.row = 1, extra.col = 2, output = "nancy_cats.geno")
# library(LEA)
# obj.snmf = snmf("nancy_cats.geno", K = 3, alpha = 100, project = "new")
# qmatrix = Q(obj.snmf, K = 3)
# 
# barplot(t(qmatrix), col = c("orange","violet","lightgreen"), border = NA, space = 0,
#         xlab = "Individuals", ylab = "Admixture coefficients")

install.packages("dartR")
library(dartR)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SNPRelate")

gl.install.vanilla.dartR()


genind_rpbb <- readRDS("./analyses/analyses_output/03a_rpbb_femaleNOknown_NOSibs_genind.Rdata")

genlight_rpbb <- gi2gl(genind_rpbb)

stc_rpbb <- gl.run.structure(genlight_rpbb, exec = "~/Downloads/console/structure", k.range = 2:5, num.k.rep = 3)

#sr <- gl.run.structure(bc, k.range = 2:5, num.k.rep = 3, 
# exec = './structure.exe')
ev <- gl.evanno(stc_rpbb)
ev

qmat <- gl.plot.structure(stc_rpbb, K=3, colors_clusters = c("red", "blue", "green", "yellow", "black"))


###########

# bc <- bandicoot.gl[,1:100]
# sr <- gl.run.structure(bc, k.range = 2:5, num.k.rep = 3, exec = '~/Downloads/console/structure')
# ev <- gl.evanno(sr)
# ev
# qmat <- gl.plot.structure(sr, K=3)
# head(qmat)
# gl.map.structure(qmat, K=3, bc, scalex=1, scaley=0.5)

# make sure you have Rtools installed
if (!require('devtools')) install.packages('devtools')
# install from GitHub
devtools::install_github('ericarcher/strataG', build_vignettes = TRUE)
