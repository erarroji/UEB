#########################BIBLIOTECAS#########################
#install.packages("devtools")
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
#installed.packages()
#########################MATRIZ#########################
paraelheatFinal<-read.table("HeatMapFinal.txt",sep="\t", header = T, quote="") #tabla con alteraciones de Intervar
paraelheat<- paraelheatFinal %>% filter(type == "Mut")
View(paraelheat)
#paraelheat$ExonicFunc.refGene<-gsub("nonsynonymous SNV", "nonsynonymousSNV", paraelheat$ExonicFunc.refGene)
#paraelheat$ExonicFunc.refGene<-gsub("frameshift deletion", "frameshiftdeletion", paraelheat$ExonicFunc.refGene)
#paraelheat$ExonicFunc.refGene<-gsub("nonframeshift deletion", "nonframeshiftdeletion", paraelheat$ExonicFunc.refGene)
#paraelheat$ExonicFunc.refGene<-gsub("frameshift insertion", "frameshiftinsertion", paraelheat$ExonicFunc.refGene)
#paraelheat$ExonicFunc.refGene<-gsub("nonframeshift insertion", "nonframeshiftinsertion", paraelheat$ExonicFunc.refGene)
nameValsSample <- sort(unique(unlist(paraelheat[,"Sample"])))
nameValsGene <- sort(unique(unlist(paraelheat[,"Gene.refGene"])))
# construct 0 matrix of correct dimensions with row and column names
myMat <- matrix("", length(nameValsSample),
                length(nameValsGene), dimnames = list(nameValsSample, nameValsGene))
# fill in the matrix with matrix indexing on row and column names
myMat[as.matrix(paraelheat[c("Sample",
                             "Gene.refGene")])] <- as.character(paraelheat[["consequence"]])

myMat
newdata<- paraelheat[order(paraelheat$gene),]
levels(paraelheat$consequence)
#########################FUNCIONES###################################
get_type_fun_onoco = function(x) strsplit(x, ";")[[1]]

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w*0.7, h*0.7, gp = gpar(fill = "gray95", col = "gray99"))
  },
  Gain = function(x, y, w, h) {
    grid.rect(x, y, w*0.7, h*0.7, gp = gpar(fill = "#DE4D4D", col = NA)) #fill es el color de relleno y 
    #col es el color de delineado, NO MOVER COL, SOLO MODIFICAR FILL
  },
  Loss = function(x, y, w, h) {
    grid.rect(x, y, w*0.7, h*0.7, gp = gpar(fill = "blue", col = NA))
  },
  Frameshift = function(x, y, w, h) {
    grid.rect(x, y, w*0.7, h*0.7, gp = gpar(fill = "#42B540B2", col = NA))
  },
  Missense = function(x, y, w, h) {
    grid.rect(x, y, w*0.7, h*0.7, gp = gpar(fill = "#048ABF", col = NA))
  },
  UTR = function(x, y, w, h) {
    grid.rect(x, y, w*0.7, h*0.7, gp = gpar(fill = "#F2D027", col = NA))
  },
#  HotSpot = function(x, y, w, h) {
#    grid.rect(x, y, w*0.3, h*0.3, gp = gpar(col = "red", lwd = 20 #col= NA#lwd = 19
 #                                           ))
#  }, 
  HotSpot = function(x, y, w, h) 
    grid.points(x, y, pch = 16)
#  Tx = function(x, y, w, h) {
#    grid.rect(x, y, w*0.6, h*0.6, gp = gpar(fill = NA, lwd = 2
#    ))
#  }
)



########################COLORES###########################
colores=c("Gain"="#DE4D4D",
          "Frameshift"="#42B540B2",
          "Loss" = "blue",
          "Missense" = "#048ABF",
          "UTR" = "#F2D027"
          )
#
########################ANOTACIONESFILAS########################
#lab13<-cbind(paraelheat["Gene.refGene"], paraelheat["N19_Lab13"])
#lab13<-lab13[order(lab13$Gene.refGene),]
#lab13<-lab13[!duplicated(lab13), ]
#lab13<-(lab13[,"N19_Lab13"])
geneRole<-cbind(newdata["Gene.refGene"],newdata["gene_role"])
geneRole<-geneRole[!duplicated(geneRole), ]
geneRole<-(geneRole[,"gene_role"])

AF<-cbind(newdata["Gene.refGene"],newdata["AF"])
AF<-AF[!duplicated(AF), ]
AF<-(AF[,"AF"])

La<-rowAnnotation(GeneRole=geneRole,
                  col=list(GeneRole=c(ambiguous="#FFE66D",OG="#FF6B6B",TSG="#0191D4")),
                           annotation_height = unit(c(5,7,2,2,10,8), "mm"),
                           gp = gpar(col = "gray"),
                           annotation_name_side = "top",
                  na_col = "gray95",
                  border=T)
                  
                  
                  
Ra<-rowAnnotation("Low AF"= AF,
                  col=list(Angiogenesis= c(Low= "#17B57A")),
                           #Death_evasion= c(Death_evasion= "#17B57A"),
                           #Cell_cycle= c(Cell_cycle= "#17B57A"),
                           #Metabolism= c(Metabolism= "#17B57A"),
                           #DNA_Repair= c(DNA_Repair= "#17B57A"),
                           #Proliferation= c(Proliferation= "#17B57A"),
                           #Migration= c(Migration= "#17B57A"),
                           #Immunology= c(Immunology= "#17B57A")),
                  annotation_height = unit(c(5,7,2,2,10,8), "mm"),
                  show_legend = c(F),#F,F,F,F,F,F,F,F,F),
                  border =T,
                  gp = gpar(col = "gray"),
                  annotation_name_side = "top",
                  na_col = "gray95")
levels(immuneact)

#ha = rowAnnotation(TumSupp=sup,
#                   OncoGene=oncog,
#                   col = list(TumSupp = c(VERDADERO= "#685C79",FALSO="#DA727E"),
#                              OncoGene = c(VERDADERO= "#685C79",FALSO="#DA727E")),
#                   annotation_height = unit(c(5,7,2,2,10,8), "mm"), gp = gpar(col = "white"),
#                   annotation_name_side = "top",
#                   annotation_legend_param = list(TumSupp =list(labels = "TumSupp"),
#                                                  OncoGene =list(labels = "OncoGene"),
 #                                                 nrow = 3, title_position = "topcenter", legend.cex=0.7))

#

########################ANOTACIONESCOLUMNAS########################
clinicos<-read.table("ConsolidadoDataBase.tsv",sep = "\t", quote = "", header = T)
clinicos<-subset(clinicos, !is.na(Oncoprint))
clinicos<-clinicos[order(clinicos$ID),]
Edad<-(clinicos[,"EDADALMOMENTODEDIAGNOSTICO"])
ECDx<-(clinicos[,"ECALMOMENTODEDIAGNOSTICO"])
Burden<-(clinicos[,"TMBMb"])
SG<-(clinicos[,"SGAB"])
OS2<-(clinicos[,"OS"])
SG<-gsub("BAJA","Low",SG)
SG<-gsub("ALTA","High",SG)
summary(clinicos$TMB)
Sig3 <-(clinicos[,"Sig3"])
Sig3<-gsub("FALSO","absent",Sig3)
Sig3<-gsub("VERDADERO","present",Sig3)

PcR <-(clinicos[,"PcR"])
PcR<-gsub("NO","incomplete", PcR)
PcR<-gsub("SI","complete", PcR)

Tx <-(clinicos[,"TraslationalTx"])
#Tx<-gsub("absence","absent", Tx)
#Tx<-gsub("SI","available", Tx)

platneo <-(clinicos[,"platneo"])
platneo<- gsub("NO", "CCT", platneo)
platneo<- gsub("SI", "Ch_PLA", platneo)

Immunoact <-(clinicos[,"Immunology_Signatures"])
TILs<-(clinicos[,"TILSesag"])
levels(Immunoact)
#Top annotation
ta<-HeatmapAnnotation(TMB = anno_barplot(Burden, width = unit(2, "cm"), gp = gpar(fill="pink")),Age=anno_points(Edad),
                      OS = anno_barplot(OS2, width = unit(2, "cm"), gp = gpar(fill="#FFE66D")))
decorate_annotation("annotation_Age_1", {
  grid.lines(c(0, 1), unit(c(40, 40), "native"), gp = gpar(col = "red"))
})
#bottom annotation
ba = HeatmapAnnotation(Stage=ECDx,
                       OS=SG,
                       Sig3=Sig3,
                       PcR=PcR,
                       Target_Treatment=Tx,
                       AdjCh=platneo,
                       #Imm_Signature=Immunoact,
                       TILs=TILs,
                       col = list(Stage=c(IIB="#0455BF", IIIA="#033E8C",
                                          IIIB="#F24405",IIIC="#F21905"),
                                  OS=c(High="#E5446D",Low="#0B4F6C"),
                                  Sig3=c(absent="#01BAEF",present="#FF8966"),
                                  PcR=c(incomplete="#7768AE",complete="#E1BC29"),
                                  Target_Treatment=c(absent="#7768AE",available="#3BB273"),
                                  AdjCh=c(CCT="#FF5964",Ch_PLA="#6BF178"),
                                  #Imm_Signature=c(Interferon_g="#d11d34",T_cell_inflamed="#4272ca",
                                   #               Th1_response="#17b57a", X="white"),
                                  TILs= c(L="#4272ca", H="#d11d34")),
                       annotation_height = unit(c(5,7,2,2,10,8), "mm"), gp = gpar(col = "gray"),
                       annotation_name_side = "left",
                       border =T# c(Stage = T, OS = T)#,
                       #annotation_legend_param = list(EdoDx =list(labels = "Edo Dx"),
                       #nrow = 3, legend.cex=0.7)
)
#########################ONCOPRINT###################################
opt<-oncoPrint(t(myMat), get_type = get_type_fun_onoco, alter_fun = alter_fun, col=colores,
          show_column_names = TRUE,
          column_names_gp = gpar(fontsize = 10), # No funciona
          row_names_gp = gpar(fontsize=10, col= "black",fontface="italic"),
          pct_gp = gpar(fontsize = 10),
          left_annotation = La,
          #top_annotation = ta,
          #bottom_annotation = ba,
          #right_annotation = Ra
          )
opt
#          bottom_annotation = ba)
#png(file = "graph25.png", width= 16000 ,height= 10000, res=1200)
#draw(opt, merge_legend = TRUE, ht_gap = unit(0.1, "cm"))
#draw(opt)
draw(opt,heatmap_legend_list = list(
  Legend(labels = "HotSpot", type = "points", pch = 20)
))

#9.67 x 10.63 Lanscape

  #¢dev.off()
#10.27*14.67
########################HEATMAPSignatures###############
Firmas<-read.table("ProyectoColombia/Firmas.txt", sep="\t", quote = "", header = T, row.names = 1)
#CosmicSig<-Firmas$Cosmic
FirmasM = subset(Firmas, select = -c(Cosmic) )
library(circlize)
col_fun <- colorRamp2(c(0, 5), c("gray99", "red"))

#RAF<-rowAnnotation(FirmasTotal=anno_barplot(firmasfrec))
#CAF<-HeatmapAnnotation(SampleSig = anno_barplot(Muestrasfrec, width = unit(2, "cm")))


ClustersFirmas <- data.frame(cutree(fit, k = 5))
ColoresClusters <- ClustersFirmas[,1]
ColoresClusters <- gsub(1, "red", ColoresClusters)
ColoresClusters <- gsub(2, "blue", ColoresClusters)
ColoresClusters <- gsub(3, "springgreen4", ColoresClusters)
ColoresClusters <- gsub(4, "orange", ColoresClusters)
ColoresClusters <- gsub(5, "purple", ColoresClusters)
  
  
Heatmap(Firmas,col = col_fun, name = "Proportion",
        cluster_rows = FALSE,
        column_dend_height = unit(2.5, "cm"),
        cluster_columns = fit,
        show_heatmap_legend = F,
        column_names_side = "top",
        column_names_gp = gpar(fontsize=10, 
                               col = ColoresClusters#c(rep("red", 27), rep("blue", 23), rep("green", 30))
        ),
        #column_km = 10,
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
  grid.text(round(Firmas[i, j],2), x, y)
}
#, right_annotation = RAF, top_annotation = CAF
)

install.packages("pvclust")
library(pvclust)
#METHODS <- c("euclidean", "maximum", "manhattan", "canberra", 
#             "binary", "minkowski")

#d<-pvclust(FirmasM, method.dist="cor", method.hclust="average", nboot=1000)
#plot(d)
#pvrect(d, alpha=0.95)
#seplot(d)

View(d)
d<-dist(t(Firmas), method = "euclidian")
fit <- hclust(d, method="ward.D") 
plot(fit) # display dendogram
install.packages("mclust")


library(mclust)
#fit <- Mclust(FirmasM)
#plot(fit) # plot results 
#summary(fit)
#

CNV<-read.table("CNVHeatMap.txt", sep="\t", quote = "", header = T)
summary(CNV$cn)
nameValsSampleCNV <- sort(unique(unlist(CNV[,"Sample"])))
nameValsGeneCNV <- sort(unique(unlist(CNV[,"gene"])))
# construct 0 matrix of correct dimensions with row and column names
myMatCNV <- matrix(0, length(nameValsSampleCNV),
                length(nameValsGeneCNV), dimnames = list(nameValsSampleCNV, nameValsGeneCNV))
# fill in the matrix with matrix indexing on row and column names
myMatCNV[as.matrix(CNV[c("Sample",
                             "gene")])] <- (CNV[["cn"]])
str(CNV)
myMatCNV



Heatmap(myMatCNV,col = col_fun, name = "Proportion",
        cluster_rows = FALSE,
        column_dend_height = unit(2.5, "cm"),
        #cluster_columns = fit,
        show_heatmap_legend = F,
        column_names_side = "top",
        column_names_gp = gpar(fontsize=10#, 
                               #col = ColoresClusters#c(rep("red", 27), rep("blue", 23), rep("green", 30))
        ),
        #column_km = 10,
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(round(myMatCNV[i, j],2), x, y)
        }
        #, right_annotation = RAF, top_annotation = CAF
)


paraelheatFinal


nameValsSampleFinal <- sort(unique(unlist(paraelheatFinal[,"Sample"])))
nameValsGeneFinal <- sort(unique(unlist(paraelheatFinal[,"Gene.refGene"])))
# construct 0 matrix of correct dimensions with row and column names
myMatFinal <- matrix("", length(nameValsSampleFinal),
                   length(nameValsGeneFinal), dimnames = list(nameValsSampleFinal, nameValsGeneFinal))
# fill in the matrix with matrix indexing on row and column names
myMatFinal[as.matrix(paraelheatFinal[c("Sample",
                         "Gene.refGene")])] <- as.character(paraelheatFinal[["consequence2"]])






optFinal<-oncoPrint(t(myMatFinal), get_type = get_type_fun_onoco, alter_fun = alter_fun, col=colores,
               show_column_names = TRUE,
               column_names_gp = gpar(fontsize = 10), # No funciona
               row_names_gp = gpar(fontsize=10, col= "black",fontface="italic"),
               pct_gp = gpar(fontsize = 10)#,
               #left_annotation = La,
               #top_annotation = ta,
               #bottom_annotation = ba,
               #right_annotation = Ra
)
optFinal
#          bottom_annotation = ba)
#png(file = "graph25.png", width= 16000 ,height= 10000, res=1200)
#draw(opt, merge_legend = TRUE, ht_gap = unit(0.1, "cm"))
#draw(opt)
draw(opt,heatmap_legend_list = list(
  Legend(labels = "HotSpot", type = "points", pch = 20)
))



TMB<-read.table("TMB.txt", sep="\t", quote = "", header = T)
TMB <- data.frame(table(TMB$Sample))
str(TMB)
TMB$KB <- TMB$Freq*0.8
median(TMB$Freq)
ggplot(TMB, aes(reorder(Var1, -KB), KB)) + geom_col(fill = "blue")+ theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Mutations by Mb") +
  geom_hline(yintercept=348.5 , color='coral', size=1) + xlab("")+
  annotate(geom="text", x=25, y=360, label="Mean", color="black")






