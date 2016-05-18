install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)







# OLD VERSION:
# row_EGFR = HeatmapAnnotation(LAML_preds_muts$EGFR,
#                              name="EGFR",
#                              col = list(EGFR.y = c("T638M" =  "black"), c("NaN" =  "red")),
#                              which = "row",
#                              width = unit(1.333, "cm"),
#                              show_legend = F)

#works
row_EGFR <- HeatmapAnnotation(LAML_preds_muts[,"EGFRgene",drop=F],
                              name="EGFR",
                              col = list(EGFRgene = c("T638M" =  "black",
                                                      "NaN" =  "white",
                                                      "T1141S" = "blue")),
                              which = "row",
                              width = unit(0.5, "cm"),
                              show_legend = F)

row_EGFR
h1 <- Heatmap(LAML_preds_muts_test[1:3], cluster_rows = T,show_column_dend = T,show_row_dend = FALSE,
              cluster_columns = T, show_row_names = F, show_column_names = T,
              row_title_gp = gpar(fontsize =10), combined_name_fun = NULL,
              name="Scaled\nPathway\nActivity",col=my_palette, column_title = "ACC",
              show_heatmap_legend = T)

#try with 0 or 1
row_EGFR <- HeatmapAnnotation(LAML_preds_muts_test[,c("KRAS_gene", "EGFR_gene" ,"BRAF_gene"),drop=F],
                              name="KRAS",
                              col = list(KRAS_gene = c("Mutated" =  "black", "Not Mutated" =  "deeppink"),
                                            EGFR_gene = c("Mutated" =  "black","Not Mutated" =  "steelblue1"),
                                             BRAF_gene = c("Mutated" =  "black", "Not Mutated" =  "olivedrab1")),
                               which = "row",
                              width = unit(0.8, "cm"),
                              show_legend = T)

draw(h1+row_EGFR)

head(LAML_preds_muts) 

summary(LAML_preds_muts[,"KRASgene",drop=F])



LAML_preds_muts_test=LAML_preds_muts
colnames(LAML_preds_muts_test)=c( "EGFR", "KRAS", "RAF1", "BRAF_gene", "EGFR_gene", "KRAS_gene", "RAF1_gene")
colnames(LAML_preds_muts_test)
LAML_preds_muts_test[4:7]
LAML_preds_muts_test[4:7] <- lapply(LAML_preds_muts_test[4:7], as.character)
LAML_preds_muts_test[4:7][LAML_preds_muts_test[4:7]=="NaN"]<-"Not Mutated"
LAML_preds_muts_test[4:7][LAML_preds_muts_test[4:7]!="Not Mutated"]<-"Mutated"
LAML_preds_muts_test

Done LAML: 






