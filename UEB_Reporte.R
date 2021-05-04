##############-----------Bibliotecas-----------##############
library(ggplot2)
library(reshape2)
library(dplyr)
library(plotly)
library(tidyr)
library(gridExtra)
library(ggpubr)
library(lemon)
library(GGally)
library(fmsb)
library(plotrix)
library(lookup)
library(ComplexHeatmap)
library(pvclust)
library(mclust)
library(circlize)
#
##############-----------FastQc-----------##############
#------------------------Data--------------------------#
FastQc <- read.table("mqc_fastqc_per_base_sequence_quality_plot_1.csv", header = T,
                     quote = "")
head((FastQc))
#------------------------Data Wrangling------------------------#
FastQcM <- FastQc %>%
    pivot_longer(-Sample,
                 names_to = 'SeqReads',
                 values_to = 'Qc') %>%
                  separate(SeqReads,
                           into = c("Letter", "SeqReads"),
                           remove = T,
                           sep = "(?<=[A-Za-z])(?=[0-9])" # regex
                           )

#------------------------Plot------------------------#
FastQcbySample <-  ggplot(FastQcM, aes(Sample ,Qc, fill= Sample)) +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = c(29, 20, -Inf), ymax = c(Inf, 29, 20),
           fill = c("#c3e6c3", "#e6dcc3", "#e6c3c3"), colour = "black", alpha = 0.65) +
  geom_boxplot(outlier.shape = NA) +
  theme_linedraw() +
    ylab("Phred score") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.9, hjust=1)) +
  expand_limits(y=c(19, 40))
  #stat_summary(fun = mean, geom = "line",
  #             aes(group = 1),
  #             position = position_dodge(width = 0.9))  #+ 
  #stat_summary(fun=mean, geom="point")


FastQcbyRead <- ggplot(FastQcM, aes(as.numeric(SeqReads), Qc, group = Sample)) +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = c(29, 20, -Inf), ymax = c(Inf, 29, 20),
           fill = c("#c3e6c3", "#e6dcc3", "#e6c3c3"), colour = "black", alpha = 0.65) +
  geom_line() +
  theme_linedraw() +
  ylab("Phred score") +
  geom_hline(yintercept = mean(as.numeric(FastQcM$SeqReads))) +
  ylim(19, 41) +
  xlab("Position (bp)") #+
  #expand_limits(y = c(19, 40), x = c(0, 150))

ggarrange(FastQcbyRead, FastQcbySample, ncol = 1, nrow = 2)

##############-----------Hs Metrics-----------##############
#------------------------Data--------------------------#
AlMet <- read.table("UEB.Metrics.Merge.txt", header = T, sep = "\t",
                    quote = "")
head(AlMet)
#------------------------Data Wrangling----------------#
AlMet.1 <- melt(AlMet[, c("Sample", "TissueType", 
                     "MEAN_TARGET_COVERAGE", 
                     "PCT_EXC_DUPE", "PCT_EXC_OFF_TARGET")])
AlMet.2 <- melt(AlMet[, c("Sample", "TissueType", 
                          "PCT_TARGET_BASES_1X", 
                          "PCT_TARGET_BASES_2X", 
                          "PCT_TARGET_BASES_10X",
                          "PCT_TARGET_BASES_20X",
                          "PCT_TARGET_BASES_30X",
                          "PCT_TARGET_BASES_40X",
                          "PCT_TARGET_BASES_50X",
                          "PCT_TARGET_BASES_100X")])
head(AlMet.1)
head(AlMet.2)
#------------------------Plot--------------------------#
AlMet.1$Sample <- factor(AlMet.1$Sample,
                         levels = AlMet.1$Sample[order(AlMet$MEAN_TARGET_COVERAGE, decreasing = T)])



variable_names <- list(
  "MEAN_TARGET_COVERAGE" = "Mean Target Coverage",
  "PCT_EXC_DUPE" = "Duplicates (%)",
  "PCT_EXC_OFF_TARGET" = "Off Target (%)"
)
variable_labeller <- function(variable,value){
  return(variable_names[value])
}

Alin.1 <- 
ggplot(AlMet.1, aes(Sample, value, fill = TissueType)) +
  geom_col() + facet_wrap(~ variable, ncol = 1, scales = "free_y",
                          labeller=variable_labeller
                          ) +
  theme_linedraw() + ylab("") +
  theme(#legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.9, hjust=1))+
  scale_fill_manual("Tissue type",
                    values = c("Normal" = "#D9A404", "Tumor" = "#D90B1C"))

Alin.1 <- Alin.1 + theme(legend.position="bottom")
max(AlMet.2$value)
Alin.2 <- 
ggplot(AlMet.2, aes(variable, (value * 100), group = Sample, color = TissueType)) +
  geom_point() + geom_line() +
  scale_x_discrete(labels=c("PCT_TARGET_BASES_1X" = "1 X",
                            "PCT_TARGET_BASES_2X" = "2 X",
                            "PCT_TARGET_BASES_10X" = "10 X",
                            "PCT_TARGET_BASES_20X" = "20 X",
                            "PCT_TARGET_BASES_30X" = "30 X",
                            "PCT_TARGET_BASES_40X" = "40 X",
                            "PCT_TARGET_BASES_50X" = "50 X",
                            "PCT_TARGET_BASES_100X" = "100 X"))+
  xlab("Depth") + ylab("Coverage depth (%)") +
  theme_linedraw() +
  scale_color_manual("Tissue type",
                    values = c("Normal" = "#D9A404", "Tumor" = "#D90B1C"))

Alin.2 <- Alin.2 + theme(legend.position="bottom")

head(AlMet.2)
AlMet.3 <- AlMet.2 %>%
    group_by(variable, TissueType) %>%
    summarise(mean = mean(value))
AlMet.3$mean <- AlMet.3$mean * 100

AlMet.3 <- as.data.frame(matrix(AlMet.3$mean, ncol = 8))
colnames(AlMet.3) <- c("d1_X", "d2_X", "d10_X", "d20_X",
                       "d30_X", "d40_X", "d50_X", "d100_X") 
rownames(AlMet.3) <- c("Normal", "Tumor")

AlMet.3 <- rbind(rep(100, 8), rep(0, 8), AlMet.3)

radarchart(AlMet.3, axistype = 1,
            #custom polygon
            pcol = c("#D9A404", "#D90B1C"),
            #pfcol = c("#D9A404", "#D90B1C"),
            plwd = 4,
            plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey",
            #caxislabels=seq(0, 100, 50), 
            cglwd=0.8,
            #custom labels
            #vlcex=0.8 
)

# Add a legend
legend(x = 1.5, y=1, legend = rownames(AlMet.3[-c(1,2),]),
       bty = "n", pch=20 , col = c("#D9A404", "#D90B1C"),
       text.col = "grey", cex=1.2, pt.cex=3)

##############-----------Filtrado-----------##############
#------------------------Data--------------------------#
Filtrado <- data.frame(Parametros = c("Total", "Paso1Pass", "Paso2Depth",
                                      "Paso3AF", "Paso4GnomAD",
                                      "Paso5NoInter", "Paso7NotSyn",
                                      "Paso9Final", 
                                      "Paso9Final"),
                       values = c(2082690,
                                  1474940,
                                  1017252,
                                  122409,
                                  22347,
                                  11768,
                                  11724,
                                  25,
                                  18),
                       Grupo = c(rep("PasoInd", 7), "Paso9a", "Paso9b")
                       )
Filtrado$Parametros <- 
  factor(Filtrado$Parametros, levels = unique(Filtrado$Parametros))

Filtrado.1 <- Filtrado %>%
    filter(values > 11723)
Filtrado.2 <- Filtrado %>%
    filter(values < 11720)
Filtrado.2[,3] <- as.character(Filtrado.2[,3])
Filtrado.2[1,3] <- "Predicted"
Filtrado.2[2,3] <- "ClinVar"
Filtrado.2[,3] <- factor(Filtrado.2[,3], levels = c("Predicted", "ClinVar"))
#------------------------Plot--------------------------#
FiltradpPlot <- 
    ggplot(Filtrado.1, aes(Parametros, values, fill = Grupo)) +
      geom_col(fill = "#90CBFB") +
  geom_text(aes(label = values), position=position_dodge(width=0.9), vjust=-0.25) +
      theme_linedraw() +
    theme(#legend.position = "none",
      axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
      scale_fill_manual("Variables",
                        values = c("Paso9a" = "#D9A404", "Paso9b" = "#D90B1C",
                                   "PasoInd" = "blue"))# +
    
        #ylim(0, 80)

FiltradpPlot.1 <- 
  ggplot(Filtrado.2, aes(Parametros, values, fill = Grupo)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = values), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_linedraw() +
  #theme(#legend.position = "none",
    #axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
  scale_fill_manual("Variables",
                    values = c("Predicted" = "#D9A404", "ClinVar" = "#D90B1C")) +
  xlab("") + ylab("# mutaciones")


FiltradpPlot + annotation_custom(ggplotGrob(FiltradpPlot.1), xmin = 4, xmax = 6.5, 
                                 ymin = 990000, ymax = 2200000)
##############-----------TMB-----------##############
#------------------------Data--------------------------#
TMB <- read.table("TMB.txt", header = T, quote = "",
                  sep = "\t")
#------------------------Data Wrangling----------------#
TMB <- as.data.frame(table(TMB))
TMB$Freq <- TMB$Freq / 800
#------------------------Plot--------------------------#
TMBplot <- 
ggplot(TMB, aes(reorder(TMB, -Freq), Freq)) +
  geom_col(aes(fill = Freq > 0.7549662), alpha = 0.8
    ) +
  geom_hline(yintercept = mean(TMB$Freq), linetype = 2, color = "red") +
  annotate("text", x = 32, y = 0.78, label = "Mean 0.7549662") +
  theme_linedraw() +
  theme(legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.9, hjust=1)) +
  xlab("Tumors") + ylab("TMB/Kb") +
  scale_fill_manual(values = c("purple", "brown"))

##############-----------Frecuencia Mutaciones-----------##############
#------------------------Data--------------------------#
GenesFreq <- read.table("mutation_analysis.tsv", sep = "\t", quote = "",
                        header = T)
#------------------------Data Wrangling----------------#
GenesFreq.1 <- 
  GenesFreq %>%
    filter(gene_role != "") %>%
    group_by(gene, gene_role) %>%
    summarise(sum = n())
#------------------------Plot--------------------------#
GenesFreq.plot <- ggplot(GenesFreq.1, aes(reorder(gene, -sum), sum, fill = gene_role)) +
    geom_col() +
  theme_linedraw() +
  theme(#legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.9, hjust=1)) +
  xlab("") + ylab("Frequency") +
  scale_fill_manual("Gene role", values = c("Act" = "#D90404", ambiguous = "#F29422",
                               LoF = "#0477BF"))
GenesFreq.plot <- GenesFreq.plot + theme(legend.position = "bottom")

##############-----------Oncoprint-----------##############
#-----------------------Data--------------------------#
GenesFreq.2 <- 
  GenesFreq %>%
  filter(gene_role != "") 

View(GenesFreq)
#------------------------Data Wrangling----------------#
nameCol <- (unique(unlist(GenesFreq.2[,"sample"])))
nameRow <- (unique(unlist(GenesFreq.2[,"gene"])))
#nameRow <- factor(FDA$Gene.refgene, levels = OrdenGenes)

# construct 0 matrix of correct dimensions with row and column names
myMat <- matrix("", length(nameRow),
                length(nameCol), dimnames = list(nameRow, nameCol))
# fill in the matrix with matrix indexing on row and column names
myMat[as.matrix(GenesFreq.2[c("gene",
                      "sample")])] <- as.character(GenesFreq.2[["consequence"]])
#write.table(myMat,"MatTx.txt", row.names = T, col.names = T, quote = F, sep = "\t")

#myMat <- read.table("MatTx.txt", sep = '\t',
#                    quote = '', header = T, row.names = 1)

#------------------------Plot--------------------------#
#------------------------Function--------------------------#
get_type_fun_onoco = function(x) strsplit(x, ";")[[1]]

levels(GenesFreq.2$consequence)

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w*0.7, h*0.7, gp = gpar(fill = "gray95", col = "gray99"))
  },
  "InFrameInsertion" = function(x, y, w, h) {
    grid.rect(x, y, w*0.7, h*0.7, gp = gpar(fill = "#1D495E", col = NA)) #fill es el color de relleno y 
    #col es el color de delineado, NO MOVER COL, SOLO MODIFICAR FILL
  },
  "Missense"  = function(x, y, w, h) {
    grid.rect(x, y, w*0.7, h*0.7, gp = gpar(fill = "#8FAD39", col = NA))
  },
  "IntronicSNV"  = function(x, y, w, h) {
    grid.rect(x, y, w*0.7, h*0.7, gp = gpar(fill = "#E1934E", col = NA))
  },
  "Frameshift"  = function(x, y, w, h) {
    grid.rect(x, y, w*0.7, h*0.7, gp = gpar(fill = "blue", col = NA))
  },
  "Nonsense" = function(x, y, w, h) {
    grid.rect(x, y, w*0.7, h*0.7, gp = gpar(fill = "red", col = NA))
  },
  "Nonframeshift substitution" = function(x, y, w, h) {
    grid.rect(x, y, w*0.7, h*0.7, gp = gpar(fill = "#BE311B", col = NA))
  },
  "Nonframeshift insertion" = function(x, y, w, h) {
    grid.rect(x, y, w*0.7, h*0.7, gp = gpar(fill = "purple", col = NA))
  },
  "Frameshift substitution" = function(x, y, w, h) {
    grid.rect(x, y, w*0.7, h*0.7, gp = gpar(fill = "#681B16", col = NA))
  },
  "Stoploss" = function(x, y, w, h) {
    grid.rect(x, y, w*0.7, h*0.7, gp = gpar(fill = "yellow", col = NA))
  },
  "Del 9-12" = function(x, y, w, h) {
    grid.rect(x, y, w*0.7, h*0.7, gp = gpar(fill = "#F6C3CB", col = NA))
  }
)

colores=c("InFrameInsertion" = "#1D495E",
          "Missense" = "#8FAD39",
          "IntronicSNV" = "#E1934E",
          "Frameshift" = "blue",
          "Nonsense" = "red",
          "Nonframeshift substitution" = "#BE311B",
          "Nonframeshift insertion" = "purple",
          "Frameshift substitution" = "#681B16",
          "Stoploss" = "yellow",
          "Del 9-12" = "#F6C3CB"
)
#------------------------Row annotation--------------------------#
Generole <- cbind(GenesFreq.2[,c("gene", "gene_role")])
#Reportables$Reportables <- factor(Reportables$Reportables, levels = c("Reportables", "No Reportables"))
Generole <- Generole[!duplicated(Generole$gene),]
#Generole <- with(Generole, Generole[order(gene_role, gene),])
Generole <- Generole[,"gene_role"]

Ra<-rowAnnotation("Gene role" = Generole,
                  col=list("Gene role" = c(Act = "#D90404", LoF = "#0477BF",
                                           ambiguous = "#F29422")),
                  annotation_height = unit(c(5,5,2,2,10,8), "mm"),
                  gp = gpar(col = "gray"),
                  annotation_name_side = "bottom",
                  show_legend = c(F),
                  na_col = "gray95",
                  border=T)


#------------------------Oncoprint--------------------------#
opt<-oncoPrint(myMat, 
               get_type = get_type_fun_onoco,
               alter_fun = alter_fun, col=colores,
               show_column_names = T,
               #column_names_side = "top",
               #column_names_rot = 55,
               #heatmap_legend_param = list(title ="Alteraciones"),
               column_names_gp = gpar(fontsize = 10),
               #column_order = orden, 
               #row_order = OrdenGenes,
               row_names_gp = gpar(fontsize = 10, col= "black",fontface="italic"),
               pct_gp = gpar(fontsize = 10),
               show_heatmap_legend=T,
               #row_split = Reportables,
               #column_split = rep(c("NGS", "PCR", "NGS"), c(25, 16, 82)),
               column_gap = unit(3, "mm"),
               border = TRUE,
               #left_annotation = La,
               #top_annotation = TaI#, 
               #row_split = Reportables
               #bottom_annotation = TaI#,
               right_annotation = Ra
)
opt

##############-----------Firmas-----------##############
#-----------------------Data--------------------------#
Firmas <- read.table("Firmas.txt", sep = "\t", 
                     header = T, quote = "")

Firmas.1 <- as.matrix(Firmas[rowSums(Firmas[]) > 0,])
#------------------------Plot--------------------------#
col_fun <- colorRamp2(c(0, 1), c("gray99", "red"))

#RAF<-rowAnnotation(FirmasTotal=anno_barplot(firmasfrec))
#CAF<-HeatmapAnnotation(SampleSig = anno_barplot(Muestrasfrec, width = unit(2, "cm")))


ClustersFirmas <- data.frame(cutree(fit, k = 5))
ColoresClusters <- ClustersFirmas[,1]
ColoresClusters <- gsub(1, "red", ColoresClusters)
ColoresClusters <- gsub(2, "blue", ColoresClusters)
ColoresClusters <- gsub(3, "springgreen4", ColoresClusters)
ColoresClusters <- gsub(4, "orange", ColoresClusters)
ColoresClusters <- gsub(5, "purple", ColoresClusters)


Heatmap(Firmas.1,
        col = col_fun,
        name = "Proportion",
        cluster_rows = FALSE,
        column_dend_height = unit(2.5, "cm"),
        #cluster_columns = fit,
        show_heatmap_legend = F,
        column_names_side = "top",
        column_names_gp = gpar(fontsize=10#, 
                               #col = ColoresClusters#c(rep("red", 27), rep("blue", 23), rep("green", 30))
        ),
        column_km = 5,
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(round(Firmas.1[i, j],2), x, y)
        }
        #, right_annotation = RAF, top_annotation = CAF
)
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

