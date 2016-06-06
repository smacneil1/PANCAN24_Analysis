# Make the Heatmaps for PANCAN24
# Shelley Macneil
# May 18, 2016 

#load dependecies
install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)


# Input Files
input_rda<- "~/Documents/PhDProjects/PANCAN24_Analysis/Data/PANCAN24data.rda"
load(input_rda)

#output Files
pdf_file<- "~/Documents/PhDProjects/PANCAN24_Analysis/Results/HeatMaps_PANCAN24_8_thicker.pdf"

pdf(pdf_file)

# create the heatmap function



palette = colorRampPalette(c("blue1","white","red1"))(n = 299)
make_PANCAN_heatmaps=function(CancerFile, CancerT, fontSize){
  
  row_EGFR <- HeatmapAnnotation(CancerFile[,c("KRAS_gene", "EGFR_gene" ,"BRAF_gene"),drop=F],
                                name="KRAS", 
                                col = list(KRAS_gene = c("Mutated" =  "black", "Not Mutated" =  "deeppink"),
                                           EGFR_gene = c("Mutated" =  "black","Not Mutated" =  "steelblue1"),
                                           BRAF_gene = c("Mutated" =  "black", "Not Mutated" =  "olivedrab1")),
                                which = "row",
                                width = unit(1.5, "cm"), 
                                show_legend = F)

  h1 <- Heatmap(scale(CancerFile[2:4]), cluster_rows = T,show_column_dend = FALSE,show_row_dend = FALSE,
                cluster_columns = T, show_row_names = F, show_column_names = T,
                row_title_gp = gpar(fontsize =12, fontface = "bold"), combined_name_fun = NULL,
                name="Scaled\nPathway\nActivity",col=palette, column_title = CancerT, column_title_gp = gpar(fontsize = fontSize, fontface = "bold"),
                show_heatmap_legend = F, width = unit(4, "cm") )
  
  draw(h1+row_EGFR, heatmap_legend_side =  "left")
}



# call the heatmap function for each cancer 

pdf(pdf_file)
# ones we are keeping 
make_PANCAN_heatmaps(pred_muts_BLCA, "Bladder urothelial carcinoma (BLCA)", 15)   # maybe: concepect of high activity but no mutations
make_PANCAN_heatmaps(pred_muts_KIRC, "Kidney renal clear cell carcinoma (KIRC)", 15) # maybe
make_PANCAN_heatmaps(pred_muts_OV, "Ovarian serous cystadenocarcinoma (OV)", 13) # maybe looks better than LIHC
make_PANCAN_heatmaps(pred_muts_BRCA, "Breast invasive carcinoma (BRCA)", 15) # add: we know know its good, low muts/high activation

# Ras driven cancers 
make_PANCAN_heatmaps(pred_muts_READ, "Rectum adenocarcinoma (READ)", 15) # correlates with KRAS, RAF, no muts/high activation
make_PANCAN_heatmaps(pred_muts_HNSC, "Head & neck squamous cell carcinoma (HNSC)", 12) # add: 
make_PANCAN_heatmaps(pred_muts_LUAD, "Lung adenocarcinoma (LUAD)", 15) 
make_PANCAN_heatmaps(pred_muts_UCS, "Uterine carcinosarcoma (UCS)", 15) # add: KRAS correlates and shows nomut/activation 

# now do some more stats

pred_muts_BLCA

#rectum adenocarcinoma
#lung


READ_KRAS_mutated_n_pathwayactivated <- pred_muts_READ[ which(pred_muts_READ$KRAS > 0.5 & pred_muts_READ$KRAS_gene=="Mutated"), ]
READ_KRAS_mutated <- pred_muts_READ[ which(pred_muts_READ$KRAS_gene=="Mutated"), ]
READ_KRAS_pathwayactivated<-pred_muts_READ[ which(pred_muts_READ$KRAS > 0.5),  ]

dim(pred_muts_READ) #167 total 
dim(READ_KRAS_pathwayactivated) #80
(READ_KRAS_mutated) #37
dim(READ_KRAS_mutated_n_pathwayactivated) #30 

#LUAD
LUAD_KRAS_mutated_n_pathwayactivated <- pred_muts_LUAD[ which(pred_muts_LUAD$KRAS > 0.5 & pred_muts_LUAD$KRAS_gene=="Mutated"), ]
LUAD_KRAS_mutated <- pred_muts_LUAD[ which(pred_muts_LUAD$KRAS_gene=="Mutated"), ]
LUAD_KRAS_pathwayactivated<-pred_muts_LUAD[ which(pred_muts_LUAD$KRAS > 0.5),  ]

dim(pred_muts_LUAD) #541 total 
dim(LUAD_KRAS_pathwayactivated) #241
dim(LUAD_KRAS_mutated) #77
dim(LUAD_KRAS_mutated_n_pathwayactivated) #34


View(LUAD_KRAS_mutated_n_pathwayactivated)
View(LUAD_KRAS_pathwayactivated)


HNSC_KRAS_mutated_n_pathwayactivated <- pred_muts_HNSC[ which(pred_muts_HNSC$KRAS > 0.5 & pred_muts_HNSC$KRAS_gene=="Mutated"), ]
HNSC_KRAS_mutated <- pred_muts_HNSC[ which(pred_muts_HNSC$KRAS_gene=="Mutated"), ]
HNSC_KRAS_pathwayactivated<-pred_muts_HNSC[ which(pred_muts_HNSC$KRAS > 0.5),  ]

dim(pred_muts_HNSC) #504 total 
dim(HNSC_KRAS_pathwayactivated) #261
dim(HNSC_KRAS_mutated) #1
dim(HNSC_KRAS_mutated_n_pathwayactivated) #34

View(pred_muts_HNSC)


BLCA_KRAS_mutated_n_pathwayactivated <- pred_muts_BLCA[ which(pred_muts_BLCA$KRAS > 0.5 & pred_muts_BLCA$KRAS_gene=="Mutated"), ]
BLCA_KRAS_mutated <- pred_muts_BLCA[ which(pred_muts_BLCA$KRAS_gene=="Mutated"), ]
BLCA_EGFR_mutated <- pred_muts_BLCA[ which(pred_muts_BLCA$EGFR_gene=="Mutated"), ]
BLCA_RAF_mutated <- pred_muts_BLCA[ which(pred_muts_BLCA$RAF_gene=="Mutated"), ]
BLCA_KRAS_pathwayactivated<-pred_muts_BLCA[ which(pred_muts_BLCA$KRAS > 0.5),  ]
BLCA_EGFR_pathwayactivated<-pred_muts_BLCA[ which(pred_muts_BLCA$EGFR > 0.5),  ]
BLCA_RAF_pathwayactivated<-pred_muts_BLCA[ which(pred_muts_BLCA$RAF > 0.5),  ]





dim(pred_muts_BLCA) #414 total 
dim(BLCA_KRAS_pathwayactivated) #92
dim(BLCA_EGFR_pathwayactivated) #176

dim(BLCA_RAF_pathwayactivated) #161




dim(BLCA_KRAS_mutated) #1
dim(BLCA_EGFR_mutated) #1
dim(BLCA_RAF_mutated) #1
dim(BLCA_KRAS_mutated_n_pathwayactivated) #34












do_stats=function(pred_file, gene){
mutated_n_pathwayactivated <- pred_file[ which(pred_file$gene > 0.5 & pred_muts_READ$KRAS_gene=="Mutated"), ]
mutated <- pred_file[ which(pred_file$KRAS_gene=="Mutated"), ]
pathwayactivated<-pred_file[ which(pred_file$gene > 0.5),  ]
  
  print("total")
  print(dim(pred_file)) #167 total 
  print("mutated")
  print(dim(mutated)) #37
  print("activated")
  print(dim(pathwayactivated)) #80
  print("both")
  print(dim(mutated_n_pathwayactivated)) #30 
}

do_stats(pred_muts_READ, KRAS)
do_stats(pred_muts_LUAD, KRAS)

View(pred_muts_LUAD)

# 7 muts didnt have activation
# 50 had activation but nt muts

View(pathwayactivated)
merge(READ_KRAS_mutated,READ_KRAS_pathwayactivated )
merge(READ_KRAS_pathwayactivated,READ_KRAS_mutated_n_pathwayactivated )



















View(pred_muts_BLCA)



dev.off()


# maybes we took out
make_PANCAN_heatmaps(pred_muts_LIHC, "Liver hepatocellular carcinoma (LIHC)", 14) # maybe does not look too amazing 
make_PANCAN_heatmaps(pred_muts_PRAD, "Prostate adenocarcinoma (PRAD)", 15) # maybe


# Not Inlcuding 
make_PANCAN_heatmaps(pred_muts_LAML, "Acute Myeloid Leukemia (LAML)", 15) 
make_PANCAN_heatmaps(pred_muts_ACC, "Adrenocortical carcinoma (ACC)", 15) 
make_PANCAN_heatmaps(pred_muts_LGG, "Brain Lower Grade Glioma (LGG)", 15) 
make_PANCAN_heatmaps(pred_muts_CESC, "Cervical squamous cell carcinoma and endocervical adenocarcinoma (CESC)", 15) 
make_PANCAN_heatmaps(pred_muts_COAD, "Colon adenocarcinoma (COAD)", 15) 
make_PANCAN_heatmaps(pred_muts_GBM, "Glioblastoma multiforme  (GBM)", 15) 
make_PANCAN_heatmaps(pred_muts_KICH, "Kidney Chromophobe (KICH)", 15) 
make_PANCAN_heatmaps(pred_muts_KIRP, "Kidney renal papillary cell carcinoma(KIRP)", 13) 
make_PANCAN_heatmaps(pred_muts_LUSC, "Lung squamous cell carcinoma (LUSC)", 15) 
make_PANCAN_heatmaps(pred_muts_DLBC, "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma (DLBC)", 15) 
make_PANCAN_heatmaps(pred_muts_SKCM, "Skin Cutaneous Melanoma (SKCM)", 15) 
make_PANCAN_heatmaps(pred_muts_STAD, "Stomach adenocarcinoma (STAD)", 15) 
make_PANCAN_heatmaps(pred_muts_THCA, "Thyroid carcinoma (THCA)", 15) 
make_PANCAN_heatmaps(pred_muts_UCEC, "Uterine Corpus Endometrial Carcinom (UCEC)", 12) # maybe not cuz KRAS doesnot look amaizng

dev.off()

# row_UCEC <- HeatmapAnnotation(pred_muts_UCEC[,c("KRAS_gene", "EGFR_gene" ,"BRAF_gene"),drop=F],
#                               name="KRAS",
#                               col = list(KRAS_gene = c("Mutated" =  "black", "Not Mutated" =  "deeppink"),
#                                          EGFR_gene = c("Mutated" =  "black","Not Mutated" =  "steelblue1"),
#                                          BRAF_gene = c("Mutated" =  "black", "Not Mutated" =  "olivedrab1")),
#                               which = "row",
#                               width = unit(0.8, "cm"),
#                               show_legend = T)
# 
# h_UCEC <- Heatmap(scale(pred_muts_UCEC[2:4]), cluster_rows = T,show_column_dend = F,show_row_dend = FALSE,
#               cluster_columns = T, show_row_names = F, show_column_names = T,
#               row_title_gp = gpar(fontsize =12, fontface = "bold"), combined_name_fun = NULL,
#               name="Scaled\nPathway\nActivity",col=my_palette, column_title = "Uterine Corpus Endometrial Carcinom (UCEC)",
#               show_heatmap_legend = T)
# 
# row_UCEC <- HeatmapAnnotation(pred_muts_UCS[,c("KRAS_gene", "EGFR_gene" ,"BRAF_gene"),drop=F],
#                               name="KRAS",
#                               col = list(KRAS_gene = c("Mutated" =  "black", "Not Mutated" =  "deeppink"),
#                                          EGFR_gene = c("Mutated" =  "black","Not Mutated" =  "steelblue1"),
#                                          BRAF_gene = c("Mutated" =  "black", "Not Mutated" =  "olivedrab1")),
#                               which = "row",
#                               width = unit(0.8, "cm"),
#                               show_legend = T)
# 
# h_UCEC <- Heatmap(scale(pred_muts_UCS[2:4]), cluster_rows = T,show_column_dend = F,show_row_dend = FALSE,
#                   cluster_columns = T, show_row_names = F, show_column_names = T,
#                   row_title_gp = gpar(fontsize =12, fontface = "bold"), combined_name_fun = NULL,
#                   name="Scaled\nPathway\nActivity",col=my_palette, column_title = "Uterine Carcinosarcoma (UCS)",
#                   show_heatmap_legend = T)
# 
# 
# draw(h_UCEC +row_UCEC +h_UCEC + row_UCEC )
# 
# probs[order( probs$Response, probs$Probability), ]
# 
# row_EGFR <- HeatmapAnnotation(pred_muts_LAML[,c("KRAS_gene", "EGFR_gene" ,"BRAF_gene"),drop=F],
#                               name="KRAS",
#                               col = list(KRAS_gene = c("Mutated" =  "black", "Not Mutated" =  "deeppink"),
#                                          EGFR_gene = c("Mutated" =  "black","Not Mutated" =  "steelblue1"),
#                                          BRAF_gene = c("Mutated" =  "black", "Not Mutated" =  "olivedrab1")),
#                               which = "row",
#                               width = unit(0.8, "cm"),
#                               show_legend = T)
# 
