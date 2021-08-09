## single cell ATAC
library(ArchR)
library(Seurat)
library(dplyr)
library(Matrix)
library(Matrix.utils)
library(ggplot2)
library(reshape2)
library(pryr)
library(RColorBrewer)
library(cowplot)
library(readr)

inputfiles <- c(
  "pool1"="../align_atac9to14/outs/fragments.tsv.gz",
  "pool2"="../align_atac15to19/outs/fragments.tsv.gz",
  "pool3"="../align_atac29to33/outs/fragments.tsv.gz",
  "pool4"="../align_atac34to39/outs/fragments.tsv.gz"
)

ArrowFiles <- createArrowFiles(
  inputFiles = inputfiles,
  sampleNames = names(inputfiles),
  filterTSS = 4,
  filterFrags=1000,
  addTileMat=T,
  addGeneScoreMat=T
)

rad1 <- ArchRProject(
  ArrowFiles=ArrowFiles,
  outputDirectory="figures", #postCovid
  copyArrows=T
)

saveArchRProject(ArchRProj = rad1, outputDirectory = "rad1", load = FALSE)

## singlet ##
sub_cellnames <- rad1$cellNames


rad <- addCellColData(
  ArchRProj = rad1,
  data = singlets,
  name = "singlets",
  cells = sub_cellnames
)


rad2 <- addCellColData(
  ArchRProj = rad,
  data =ids,
  name = "ids",
  cells = sub_cellnames
)


select_cells <- which(rad@cellColData$singlets == "singlet")
rad2_2 <- subsetCells(ArchRProj = rad2, cellNames = sub_cellnames[which(singlets == "singlet")])

##### umap ######
rad2 <- addIterativeLSI(
  ArchRProj=proj1,
  useMatrix="TileMatrix",
  name="IterativeLSI",
  iterations=2,
  clusterParams=list(
    resolution=c(0.3),
    sampleCells=10000,
    n.start=10
  ),
  varFeatures=25000,
  dimsToUse=1:30
)

rad2 <- addHarmony(
  ArchRProj=rad2,
  reducedDims="IterativeLSI",
  name="Harmony",
  groupBy="Sample"
)

rad2 <- addClusters(
  input=rad2,
  reducedDims="IterativeLSI",
  method="Seurat",
  name="Clusters",
  resolution=1
)

rad2 <- addUMAP(
  ArchRProj=rad2,
  reducedDims="Harmony",
  name="UMAPHarmony",
  nNeighbors=30,
  minDist=0.5,
  metric="cosine"
)

saveArchRProject(ArchRProj=rad2, outputDirectory="rad2_1",load=F)

p1 <- plotEmbedding(ArchRProj = rad2, colorBy = "cellColData", name = "Sample", embedding= "UMAPHarmony")
p2 <- plotEmbedding(ArchRProj = rad2, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")

# figure 2B
plotPDF(p1,p2, name = "Plots-UMAP2-Sample-Clusters.pdf", ArchRProj = rad2,addDoc=F, width=5, height=5)

markersGS <- getMarkerFeatures(
  ArchRProj=rad2,
  useMatrix="GeneScoreMatrix",
  groupBy="Clusters",
  bias=c("TSSEnrichment","log10(nFrags)"),
  testMethod="wilcoxon"
  
)

markerList <- getMarkers(markersGS, cutOff="FDR <= 0.01 & Log2FC >= 1.25")

marker.gene <- c(
  "IL7R","CCR7",,  #IL7R,CCR7-,TCF7-             Memory CD4 T   
  "CD8A","CD8B","GZMK", #CD8A+,CD8B+,GZMK+             CD8 T    
  "NKG7","NCAM1", #NKG7+,NCAM1+,CD8A-             NK    
  "CD14","LYZ","FCGR3A", #CD14+, LYZ+, FCGR3A-            CD14+ monocytes   
  #FCGR3A+, CD14-              CD16+ monocytes   
  "CST3","CD86","HLA-DRA",  #CST3+, CD86+, HLA-DRA+, CD14-, FCGR3A-   mDC, pDC (dendritic cells)    
  "CD79A", #CD79A+, CD27-, SDC1-           B cell        
  #CD79A+, CD27+, SDC1+           plasmablast        
  "PPBP" #PPBP+                                 megakaryocyte/platelet   
)

heatmapGS <- markerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR<0.01 & Log2FC >= 1.25",
  labelMarkers=marker.gene,
  transpose=TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side="bot", annotation_legend_side="bot")
plotPDF(heatmapGS, name = "Plot-geneScores.pdf", ArchRProj = rad2,addDoc=F, width=8, height=6)

## integration with scRNA supplementary figure 3##
rna <- readRDS("../rad205.rds")
new.cluster.ids2 <- c("Classical Monocytes","CD4+ T cells","CD4+ T cells","CD4+ T cells","NK cells",
                      "CD8+ T cells","Non-classical Monocytes","NK cells","B cells","B cells","CD8+ T cells",
                      "Classical Monocytes","Classical Monocytes","Dendritic cells","Classical Monocytes","Platelet",
                      "16","17","18")


names(new.cluster.ids2) <- levels(rna@meta.data$seurat_clusters)
Idents(rna) <- rna@meta.data$seurat_clusters
rna <- RenameIdents(rna, new.cluster.ids2)

rna@meta.data$clusters1 <- Idents(rna)

rad <- rad2

rad_inter <- addGeneIntegrationMatrix(
  ArchRProj = rad, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = rna,
  addToArrow = FALSE,
  groupRNA = "clusters1",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

rad_inter2 <- addGeneIntegrationMatrix(
  ArchRProj = rad_inter, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = rna,
  addToArrow = TRUE,
  force= TRUE,
  groupRNA = "clusters1",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)

rad_inter2 <- addImputeWeights(rad_inter2)

marker.gene <- c(
  "IL7R","CCR7","TCF7",  #IL7R,CCR7-,TCF7-             Memory CD4 T   
  "CD8A","CD8B","GZMK", #CD8A+,CD8B+,GZMK+             CD8 T    
  "NKG7","NCAM1", #NKG7+,NCAM1+,CD8A-             NK    
  "CD14","LYZ","FCGR3A", #CD14+, LYZ+, FCGR3A-            CD14+ monocytes   
  #FCGR3A+, CD14-              CD16+ monocytes   
  "CST3","CD86","HLA-DRA",  #CST3+, CD86+, HLA-DRA+, CD14-, FCGR3A-   mDC, pDC (dendritic cells)    
  "IL1B","TNF","IFI6","IFITM3",  #IL1B,TNF,IL6               proinflammatory markers  # macrophages
  "IFNG","GZMB","CCL5","CST7","CLEC9A","CD1D","CXCL9","CXCL10","CXCL11",
  "CD79A", #CD79A+, CD27-, SDC1-           B cell        
  #CD79A+, CD27+, SDC1+           plasmablast        
  "PPBP", #PPBP+                                 megakaryocyte/platelet   
  "KIT","TPSAB1" #KIT+, TPSAB1+                 Mast (edited)
)

p1 <- plotEmbedding(
  ArchRProj = rad_inter2, 
  colorBy = "GeneIntegrationMatrix", 
  name = marker.gene, 
  continuousSet = "horizonExtra",
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(rad_inter2)
)

p2 <- plotEmbedding(
  ArchRProj = rad_inter2, 
  colorBy = "GeneScoreMatrix", 
  continuousSet = "horizonExtra",
  name = marker.gene, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(rad_inter2)
)

p1c <- lapply(p1, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

p2c <- lapply(p2, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})


plotPDF(plotList = p1, 
        name = "Plot-UMAP-Marker-Genes-RNA-W-Imputation.pdf", 
        ArchRProj = rad_inter2, 
        addDOC = FALSE, width = 5, height = 5)


cM <- confusionMatrix(rad_inter2$Clusters, rad_inter2$predictedGroup)
labelOld <- rownames(cM)

labelNew <- colnames(cM)[apply(cM, 1, which.max)]

rad_inter2$Clusters2 <- mapLabels(rad_inter2$Clusters, newLabels = labelNew, oldLabels = labelOld)

p1 <- plotEmbedding(rad_inter2, colorBy = "cellColData", name = "Clusters2", embedding = "UMAPHarmony")


plotPDF(p1, name = "Plots-atac-rna.pdf", ArchRProj = rad_inter2,addDoc=F, width=5, height=5)


####### marker peaks #######
rad_inter2 <- addGroupCoverages(ArchRProj = rad_inter2, groupBy = "cell_dis", force = TRUE)
pathToMacs2 <- "/home/bin/macs2"

rad_temp <- addReproduciblePeakSet(
  ArchRProj = rad_inter2, 
  groupBy = "cell_dis", 
  pathToMacs2 = pathToMacs2,
  useGroups = "mono_post",
  bgdGroups = "mono_healthy"
)


getPeakSet(rad_temp)

rad_inter2 <- addPeakMatrix(rad_temp)

marker_cm <- getMarkerFeatures(
  ArchRProj = rad5, 
  useMatrix = "PeakMatrix",
  groupBy = "cell_dis2",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Classical Monocytes_post",
  bgdGroups = "Classical Monocytes_healthy"
)

markerList_cm <- getMarkers(markerTest_cm, cutOff = "Pval <= 0.05 & Log2FC >= 0.25")
pma_cm <- markerPlot(seMarker = markerTest_cm, name = "Classical Monocytes_post", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.25", plotAs = "MA")
pma_nm <- markerPlot(seMarker = markerTest_nm, name = "Non-classical Monocytes_post", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.25", plotAs = "MA")
pma_nk <- markerPlot(seMarker = markerTest_nk, name = "NK cells_post", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.25", plotAs = "MA")
pma_t <- markerPlot(seMarker = markerTest_t, name = "CD4+ T cells_post", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.25", plotAs = "MA")
pma_b <- markerPlot(seMarker = markerTest_b, name = "B cells_post", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.25", plotAs = "MA")
pma_dc <- markerPlot(seMarker = markerTest_dc, name = "Dendritic cells_post", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.25", plotAs = "MA")


plotPDF(pma_cm, pma_nm, pma_nk, pma_t, pma_b, pma_dc,name = "post-vs-healthy-marker-motif-FDR", 
        width = 5, height = 5, ArchRProj = rad5, addDOC = FALSE)


##### heatmap #######
rad4 = rad_inter2
markersGS <- getMarkerFeatures(
  ArchRProj = rad4, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


markerList <- getMarkers(markersGS, cutOff = "FDR >= 0.25")

t_atac = as.data.frame(markerList$'CD4+ T cells'[order(markerList$'CD4+ T cells'$Log2FC, decreasing = TRUE),])
nk_atac = as.data.frame(markerList$'NK cells'[order(markerList$'NK cells'$Log2FC, decreasing = TRUE),])
cm_atac = as.data.frame(markerList$'Classical Monocytes'[order(markerList$'Classical Monocytes'$Log2FC, decreasing = TRUE),])
nm_atac = as.data.frame(markerList$'Non-classical Monocytes'[order(markerList$'Non-classical Monocytes'$Log2FC, decreasing = TRUE),])
b_atac = as.data.frame(markerList$'B cells'[order(markerList$'B cells'$Log2FC, decreasing = TRUE),])
dc_atac = as.data.frame(markerList$'Dendritic cells'[order(markerList$'Dendritic cells'$Log2FC, decreasing = TRUE),])

write_tsv(as.data.frame(markerList$'16'[order(markerList$'16'$Log2FC, decreasing = TRUE),]), "../post-covid/results/rad2/marker16_all.txt")
write_tsv(as.data.frame(markerList$'B cells'[order(markerList$'B cells'$Log2FC, decreasing = TRUE),]), "../post-covid/results/rad2/markerb_all.txt")
write_tsv(as.data.frame(markerList$'Classical Monocytes'[order(markerList$'Classical Monocytes'$Log2FC, decreasing = TRUE),]), "../post-covid/results/rad2/markercm_all.txt")
write_tsv(as.data.frame(markerList$'Non-classical Monocytes'[order(markerList$'Non-classical Monocytes'$Log2FC, decreasing = TRUE),]), "../projects/post-covid/results/rad2/markernm_all.txt")
write_tsv(as.data.frame(markerList$'CD4+ T cells'[order(markerList$'CD4+ T cells'$Log2FC, decreasing = TRUE),]), "../post-covid/results/rad2/markert_all.txt")
write_tsv(as.data.frame(markerList$'NK cells'[order(markerList$'NK cells'$Log2FC, decreasing = TRUE),]), "../post-covid/results/rad2/markernk_all.txt")
write_tsv(as.data.frame(markerList$'Dendritic cells'[order(markerList$'Dendritic cells'$Log2FC, decreasing = TRUE),]), "../post-covid/results/rad2/markerdc_all.txt")

t <- read_tsv("../rad2/rad_cell_dis2/markert_all.txt")
nk <- read_tsv("../rad2/rad_cell_dis2/markernk_all.txt")
cm <- read_tsv("../rad2/rad_cell_dis2/markercm_all.txt")
nm <- read_tsv("../rad2/rad_cell_dis2/markernm_all.txt")
b <- read_tsv("../rad2/rad_cell_dis2/markerb_all.txt")
dc <- read_tsv("../rad2/rad_cell_dis2/markerdc_all.txt")

t_rna <- read.csv("../rad2/rad_cell_dis2/DEG_no/CD4DEG.csv")
t_rna2 <- t_rna[order(t_rna$avg_logFC, decreasing = TRUE), ]
t_rna2 <- t_rna2 %>% filter(avg_logFC > 0)
t_rna2 = t_rna2$X
t_int = intersect(t$name, t_rna2)

nk_rna <- read.csv("../rad2/rad_cell_dis2/DEG_no/NKDEG.csv")
nk_rna2 <- nk_rna[order(nk_rna$avg_logFC, decreasing = TRUE), ]
nk_rna2 <- nk_rna2 %>% filter(avg_logFC > 0)
nk_rna2 = nk_rna2$X
nk_int = intersect(nk$name, nk_rna2)

cm_rna <- read.csv("../rad2/rad_cell_dis2/DEG_no/CD14DEG.csv")
cm_rna2 <- cm_rna[order(cm_rna$avg_logFC, decreasing = TRUE), ]
cm_rna2 <- cm_rna2 %>% filter(avg_logFC > 0)
cm_rna2 = cm_rna2$X
cm_int = intersect(cm$name, cm_rna2)

nm_rna <- read.csv("../rad2/rad_cell_dis2/DEG_no/CD16DEG.csv")
nm_rna2 <- nm_rna[order(nm_rna$avg_logFC, decreasing = TRUE), ]
nm_rna2 <- nm_rna2 %>% filter(avg_logFC > 0)
nm_rna2 = nm_rna2$X
nm_int = intersect(nm$name, nm_rna2)

b_rna <- read.csv("../rad2/rad_cell_dis2/DEG_no/BDEG.csv")
b_rna2 <- b_rna[order(b_rna$avg_logFC, decreasing = TRUE), ]
b_rna2 <- b_rna2 %>% filter(avg_logFC > 0)
b_rna2 = b_rna2$X
b_int = intersect(b$name, b_rna2)

dc_rna <- read.csv("../rad2/rad_cell_dis2/DEG_no/mDCDEG.csv")
dc_rna2 <- dc_rna[order(dc_rna$avg_logFC, decreasing = TRUE), ]
dc_rna2 <- dc_rna2 %>% filter(avg_logFC > 0)
dc_rna2 = dc_rna2$X
dc_int = intersect(dc$name, dc_rna2)

t_atac_idx = match(markers, t$name)
t_rna_idx = match(markers, t_rna$X)

nk_atac_idx = match(markers, nk$name)
nk_rna_idx = match(markers, nk_rna$X)

cm_atac_idx = match(markers, cm$name)
cm_rna_idx = match(markers, cm_rna$X)

nm_atac_idx = match(markers, nm$name)
nm_rna_idx = match(markers, nm_rna$X)

b_atac_idx = match(markers, b$name)
b_rna_idx = match(markers, b_rna$X)

dc_atac_idx = match(markers, dc$name)
dc_rna_idx = match(markers, dc_rna$X)

t_atac_sub = t[t_atac_idx, ]
t_rna_sub = t_rna[t_rna_idx, ]

nk_atac_sub = nk[nk_atac_idx, ]
nk_rna_sub = nk_rna[nk_rna_idx, ]

cm_atac_sub = cm[cm_atac_idx, ]
cm_rna_sub = cm_rna[cm_rna_idx, ]

nm_atac_sub = nm[nm_atac_idx, ]
nm_rna_sub = nm_rna[nm_rna_idx, ]

b_atac_sub = b[b_atac_idx, ]
b_rna_sub = b_rna[b_rna_idx, ]

dc_atac_sub = dc[dc_atac_idx, ]
dc_rna_sub = dc_rna[dc_rna_idx, ]


atac_exp = data.frame(t = t_atac_sub$Log2FC, 
                      nk = nk_atac_sub$Log2FC, 
                      cm = cm_atac_sub$Log2FC, 
                      nm = nm_atac_sub$Log2FC,
                      b = b_atac_sub$Log2FC, 
                      dc = dc_atac_sub$Log2FC)

rna_exp = data.frame(t = t_rna_sub$avg_logFC, 
                     nk = nk_rna_sub$avg_logFC, 
                     cm = cm_rna_sub$avg_logFC, 
                     nm = nm_rna_sub$avg_logFC,
                     b = b_rna_sub$avg_logFC, 
                     dc = dc_rna_sub$avg_logFC)

pheatmap(as.matrix(atac_exp), 
         color = colorRamp2(c(2.5, 0, min(as.matrix(atac_exp))), c("red", "white", "blue")),
         #color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(100),
         scale = "none", angle_col = "90", cluster_cols = FALSE, cluster_rows = FALSE,
         labels_row = markers, labels_col = c("CD4T", "NK", "CM", "NM", "B", "DC"), fontsize = 12, fontsize_col = 20)

pheatmap(as.matrix(rna_exp), color = colorRamp2(c(0.1, 0, min(as.matrix(rna_exp))), c("red", "white", "blue")),
         scale = "none", angle_col = "90", cluster_cols = FALSE, cluster_rows = FALSE,
         labels_row = markers, labels_col = c("CD4T", "NK", "CM", "NM", "B", "DC"),fontsize = 12, fontsize_col = 20)
##################################################################################
