# Make the Heatmaps for PANCAN24

#load dependecies
install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)

# Input Files
input_rda<- "~/Documents/PhDProjects/PANCAN24_Analysis/Data/PANCAN24data.rda"
load(input_rda)

#output Files
pdf_file<- "~/Documents/PhDProjects/PANCAN24_Analysis/Results/HeatMaps_PANCAN24_scaled.pdf"

pdf(pdf_file)

row_EGFR <- HeatmapAnnotation(pred_muts_LAML[,c("KRAS_gene", "EGFR_gene" ,"BRAF_gene"),drop=F],
                              name="KRAS",
                              col = list(KRAS_gene = c("Mutated" =  "black", "Not Mutated" =  "deeppink"),
                                         EGFR_gene = c("Mutated" =  "black","Not Mutated" =  "steelblue1"),
                                         BRAF_gene = c("Mutated" =  "black", "Not Mutated" =  "olivedrab1")),
                              which = "row",
                              width = unit(0.8, "cm"),
                              show_legend = T)

h1 <- Heatmap(pred_muts_LAML[2:4], cluster_rows = T,show_column_dend = T,show_row_dend = FALSE,
              cluster_columns = T, show_row_names = F, show_column_names = T,
              row_title_gp = gpar(fontsize =10, fontface = "bold"), combined_name_fun = NULL,
              name="Scaled\nPathway\nActivity",col=my_palette, column_title = "LAML",
              show_heatmap_legend = T)

draw(h1+row_EGFR)
head(pred_muts_LAML[1:4])

make_PANCAN_heatmaps=function(CancerFile, CancerT){
  
  row_EGFR <- HeatmapAnnotation(CancerFile[,c("KRAS_gene", "EGFR_gene" ,"BRAF_gene"),drop=F],
                                name="KRAS",
                                col = list(KRAS_gene = c("Mutated" =  "black", "Not Mutated" =  "deeppink"),
                                           EGFR_gene = c("Mutated" =  "black","Not Mutated" =  "steelblue1"),
                                           BRAF_gene = c("Mutated" =  "black", "Not Mutated" =  "olivedrab1")),
                                which = "row",
                                width = unit(0.8, "cm"),
                                show_legend = T)
  
  h1 <- Heatmap(scale(CancerFile[2:4]), cluster_rows = T,show_column_dend = T,show_row_dend = FALSE,
                cluster_columns = T, show_row_names = F, show_column_names = T,
                row_title_gp = gpar(fontsize =12, fontface = "bold"), combined_name_fun = NULL,
                name="Scaled\nPathway\nActivity",col=my_palette, column_title = CancerT,
                show_heatmap_legend = T)
 
 draw(h1+row_EGFR)
                 }

pred_muts_LAML
CancerFile
pdf_file
pdf(pdf_file)


make_PANCAN_heatmaps(pred_muts_LAML, "Acute Myeloid Leukemia (LAML)") 
make_PANCAN_heatmaps(pred_muts_ACC, "Adrenocortical carcinoma (ACC)")
make_PANCAN_heatmaps(pred_muts_BLCA, "Bladder Urothelial Carcinoma (BLCA)")
make_PANCAN_heatmaps(pred_muts_LGG, "Brain Lower Grade Glioma (LGG)")
make_PANCAN_heatmaps(pred_muts_BRCA, "Breast invasive carcinoma (BRCA)")
make_PANCAN_heatmaps(pred_muts_CESC, "Cervical squamous cell carcinoma and endocervical adenocarcinoma (CESC)")
make_PANCAN_heatmaps(pred_muts_COAD, "Colon adenocarcinoma (COAD)")
make_PANCAN_heatmaps(pred_muts_GBM, "Glioblastoma multiforme  (GBM)")
make_PANCAN_heatmaps(pred_muts_HNSC, "Head and Neck squamous cell carcinoma (HNSC)")
make_PANCAN_heatmaps(pred_muts_KICH, "Kidney Chromophobe (KICH)")
make_PANCAN_heatmaps(pred_muts_KIRC, "Kidney renal clear cell carcinoma (KIRC)")
make_PANCAN_heatmaps(pred_muts_KIRP, "Kidney renal papillary cell carcinoma(KIRP)")
make_PANCAN_heatmaps(pred_muts_LIHC, "Liver hepatocellular carcinoma(LIHC)")
make_PANCAN_heatmaps(pred_muts_LUAD, "Lung adenocarcinoma (LUAD)")
make_PANCAN_heatmaps(pred_muts_LUSC, "Lung squamous cell carcinoma (LUSC)")
make_PANCAN_heatmaps(pred_muts_DLBC, "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma (DLBC)")
make_PANCAN_heatmaps(pred_muts_OV, "Ovarian serous cystadenocarcinoma (OV)")
make_PANCAN_heatmaps(pred_muts_PRAD, "Prostate adenocarcinoma (PRAD")
make_PANCAN_heatmaps(pred_muts_READ, "Rectum adenocarcinoma (READ)")
make_PANCAN_heatmaps(pred_muts_SKCM, "Skin Cutaneous Melanoma (SKCM)")
make_PANCAN_heatmaps(pred_muts_STAD, "Stomach adenocarcinoma (STAD)")
make_PANCAN_heatmaps(pred_muts_THCA, "Thyroid carcinoma (THCA)")
make_PANCAN_heatmaps(pred_muts_UCS, "Uterine Carcinosarcoma (UCS)")
make_PANCAN_heatmaps(pred_muts_UCEC, "Uterine Corpus Endometrial Carcinom (UCEC)")

dev.off()



row_UCEC <- HeatmapAnnotation(pred_muts_UCEC[,c("KRAS_gene", "EGFR_gene" ,"BRAF_gene"),drop=F],
                              name="KRAS",
                              col = list(KRAS_gene = c("Mutated" =  "black", "Not Mutated" =  "deeppink"),
                                         EGFR_gene = c("Mutated" =  "black","Not Mutated" =  "steelblue1"),
                                         BRAF_gene = c("Mutated" =  "black", "Not Mutated" =  "olivedrab1")),
                              which = "row",
                              width = unit(0.8, "cm"),
                              show_legend = T)

h_UCEC <- Heatmap(scale(pred_muts_UCEC[2:4]), cluster_rows = T,show_column_dend = F,show_row_dend = FALSE,
              cluster_columns = T, show_row_names = F, show_column_names = T,
              row_title_gp = gpar(fontsize =12, fontface = "bold"), combined_name_fun = NULL,
              name="Scaled\nPathway\nActivity",col=my_palette, column_title = "Uterine Corpus Endometrial Carcinom (UCEC)",
              show_heatmap_legend = T)

row_UCEC <- HeatmapAnnotation(pred_muts_UCS[,c("KRAS_gene", "EGFR_gene" ,"BRAF_gene"),drop=F],
                              name="KRAS",
                              col = list(KRAS_gene = c("Mutated" =  "black", "Not Mutated" =  "deeppink"),
                                         EGFR_gene = c("Mutated" =  "black","Not Mutated" =  "steelblue1"),
                                         BRAF_gene = c("Mutated" =  "black", "Not Mutated" =  "olivedrab1")),
                              which = "row",
                              width = unit(0.8, "cm"),
                              show_legend = T)

h_UCEC <- Heatmap(scale(pred_muts_UCS[2:4]), cluster_rows = T,show_column_dend = F,show_row_dend = FALSE,
                  cluster_columns = T, show_row_names = F, show_column_names = T,
                  row_title_gp = gpar(fontsize =12, fontface = "bold"), combined_name_fun = NULL,
                  name="Scaled\nPathway\nActivity",col=my_palette, column_title = "Uterine Carcinosarcoma (UCS)",
                  show_heatmap_legend = T)


draw(h_UCEC +row_UCEC +h_UCEC + row_UCEC )

probs[order( probs$Response, probs$Probability), ]

