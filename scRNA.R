## single cell RNA sequencing ##
library(Seurat)
library(ggplot2)
library(plyr)
library(data.table)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(SingleR)
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(org.Hs.eg.db)
library(clusterProfiler)


####### 9to14 ##########
# Reading in the cluster assignment
clusters9 <- read.table("9to14/souporcell_9to14.tsv", header = TRUE, stringsAsFactors = FALSE)
clusters9$barcode <- substr(clusters9$barcode, start = 1, stop = nchar(clusters9$barcode)-2)
rownames(clusters9) <- clusters9$barcode
# Reading in the 10X data
data9 <- Read10X("9to14/filtered_feature_bc_matrix")
doublets9 <- subset(clusters9, status == "doublet")
unassigned9 <- subset(clusters9, status == "unassigned")
clusters9 <- subset(clusters9, status == "singlet")
# We remove all cells that are not singlets
colnames(data9) <- substring(colnames(data9),1,nchar(colnames(data9))-2)
data9 <- data9[, colnames(data9) %in% clusters9$barcode]
# We calculate the percentage of MT reads
MT.index9 <- grep(pattern = "^MT-", x = rownames(data9), value = FALSE)
all.sum9 <- Matrix::colSums(data9)
percent.MT9 <- Matrix::colSums(data9[MT.index9, ])/all.sum9
data9 <- data9[-MT.index9, ] # remove MT genes
# We calculate the percentage of Ribosomal genes
RG.index9 <- grep("^RP[SL][0-9]+$", x = rownames(data9), perl=TRUE, value = FALSE)
percent.RG9 <- Matrix::colSums(data9[RG.index9, ])/all.sum9
data9 <- data9[-RG.index9, ] 
# Adding the cluster assignment 
clusters9 <- clusters9[,c(1:3)]
clusters9$sample <- as.factor(clusters9$assignment)
sampMatch9 <- c("0" = "HC04", "1" = "HC01", "2" = "PAT02", "3" = "PAT10","4" = "PAT11", "5" = "HC12")
clusters9$sample <- as.character(revalue(clusters9$sample, sampMatch9))
cluster.assignment9 <- clusters9$sample
names(cluster.assignment9) <- clusters9$barcode
cluster.assignment9 <- cluster.assignment9[order(match(names(cluster.assignment9), names(percent.MT9)))]
# Adding the status assignment  
cluster.status9 <- clusters9$status
names(cluster.status9) <- clusters9$barcode
cluster.status9 <- cluster.status9[order(match(names(cluster.status9), names(percent.MT9)))]
# Adding the meta data
data9 <- CreateSeuratObject(data9, min.cells = 5, meta.data = data.frame(percent.mt = percent.MT9, percent.rg = percent.RG9, sample = cluster.assignment9, status = cluster.status9), project = "Pool9to14")
#VlnPlot(data9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
data9 <- subset(data9, nFeature_RNA > 250 & nFeature_RNA < 5000 & percent.mt < 0.25)
rm(cluster9, doublets9, unassigned9, MT.index9, all.sum9, percent.MT9, percent.RG9, RG.index9, sampMatch9, cluster.assignment9, cluster.status9)


# Reading in the cluster assignment
clusters15 <- read.table("15to19/souporcell_15to19.tsv", header = TRUE, stringsAsFactors = FALSE)
clusters15$barcode <- substr(clusters15$barcode, start = 1, stop = nchar(clusters15$barcode)-2)
rownames(clusters15) <- clusters15$barcode
# Reading in the 10X data
data15 <- Read10X("15to19/filtered_feature_bc_matrix")
doublets15 <- subset(clusters15, status == "doublet")
unassigned15 <- subset(clusters15, status == "unassigned")
clusters15 <- subset(clusters15, status == "singlet")
# We remove all cells that are not singlets
colnames(data15) <- substring(colnames(data15),1,nchar(colnames(data15))-2)
data15 <- data15[, colnames(data15) %in% clusters15$barcode]
# We calculate the percentage of MT reads
MT.index15 <- grep(pattern = "^MT-", x = rownames(data15), value = FALSE)
all.sum15 <- Matrix::colSums(data15)
percent.MT15 <- Matrix::colSums(data15[MT.index15, ])/all.sum15
data15 <- data15[-MT.index15, ] # remove MT genes
# We calculate the percentage of Ribosomal genes
RG.index15 <- grep("^RP[SL][0-9]+$", x = rownames(data15), perl=TRUE, value = FALSE)
percent.RG15 <- Matrix::colSums(data15[RG.index15, ])/all.sum15
data15 <- data15[-RG.index15, ] 
# Adding the cluster assignment 
clusters15 <- clusters15[,c(1:3)]
clusters15$sample <- as.factor(clusters15$assignment)
sampMatch15 <- c("0" = "PAT09", "1" = "HC08", "2" = "PAT08", "3" = "HC02","4" = "PAT03")
clusters15$sample <- as.character(revalue(clusters15$sample, sampMatch15))
cluster.assignment15 <- clusters15$sample
names(cluster.assignment15) <- clusters15$barcode
cluster.assignment15 <- cluster.assignment15[order(match(names(cluster.assignment15), names(percent.MT15)))]
# Adding the status assignment  
cluster.status15 <- clusters15$status
names(cluster.status15) <- clusters15$barcode
cluster.status15 <- cluster.status15[order(match(names(cluster.status15), names(percent.MT15)))]
# Adding the meta data
data15 <- CreateSeuratObject(data15, min.cells = 5, meta.data = data.frame(percent.mt = percent.MT15, percent.rg = percent.RG15, sample = cluster.assignment15, status = cluster.status15), project = "Pool15to19")
#VlnPlot(data15, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
data15 <- subset(data15, nFeature_RNA > 250 & nFeature_RNA < 5000 & percent.mt < 0.25)
rm(clusters15, doublets15, unassigned15, MT.index15, all.sum15, percent.MT15, percent.RG15, RG.index15, sampMatch15, cluster.assignment15, cluster.status15)

############## 29to33 ##############
# Reading in the cluster assignment
clusters29 <- read.table("29to33/souporcell_29to33.tsv", header = TRUE, stringsAsFactors = FALSE)
clusters29$barcode <- substr(clusters29$barcode, start = 1, stop = nchar(clusters29$barcode)-2)
rownames(clusters29) <- clusters29$barcode
# Reading in the 10X data
data29 <- Read10X("29to33/filtered_feature_bc_matrix")
doublets29 <- subset(clusters29, status == "doublet")
unassigned29 <- subset(clusters29, status == "unassigned")
clusters29 <- subset(clusters29, status == "singlet")
# We remove all cells that are not singlets
colnames(data29) <- substring(colnames(data29),1,nchar(colnames(data29))-2)
data29 <- data29[, colnames(data29) %in% clusters29$barcode]
# We calculate the percentage of MT reads
MT.index29 <- grep(pattern = "^MT-", x = rownames(data29), value = FALSE)
all.sum29 <- Matrix::colSums(data29)
percent.MT29 <- Matrix::colSums(data29[MT.index29, ])/all.sum29
data29 <- data29[-MT.index29, ] # remove MT genes
# We calculate the percentage of Ribosomal genes
RG.index29 <- grep("^RP[SL][0-9]+$", x = rownames(data29), perl=TRUE, value = FALSE)
percent.RG29 <- Matrix::colSums(data29[RG.index29, ])/all.sum29
data29 <- data29[-RG.index29, ] 
# Adding the cluster assignment 
clusters29 <- clusters29[,c(1:3)]
clusters29$sample <- as.factor(clusters29$assignment)
sampMatch29 <- c("0" = "HC05", "1" = "HC07", "2" = "PAT07", "3" = "HC03", "4" = "PAT04")
clusters29$sample <- as.character(revalue(clusters29$sample, sampMatch29))
cluster.assignment29 <- clusters29$sample
names(cluster.assignment29) <- clusters29$barcode
cluster.assignment29 <- cluster.assignment29[order(match(names(cluster.assignment29), names(percent.MT29)))]
# Adding the status assignment  
cluster.status29 <- clusters29$status
names(cluster.status29) <- clusters29$barcode
cluster.status29 <- cluster.status29[order(match(names(cluster.status29), names(percent.MT29)))]
# Adding the meta data
data29 <- CreateSeuratObject(data29, min.cells = 5, meta.data = data.frame(percent.mt = percent.MT29, percent.rg = percent.RG29, sample = cluster.assignment29, status = cluster.status29), project = "Pool29to33")
#VlnPlot(data24, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
data29 <- subset(data29, nFeature_RNA > 250 & nFeature_RNA < 5000 & percent.mt < 0.25)
rm(clusters29, doublets29, unassigned29, MT.index29, all.sum29, percent.MT29, percent.RG29, RG.index29, sampMatch29, cluster.assignment29, cluster.status29)



clusters34 <- read.table("34to39/souporcell_34to39.tsv", header = TRUE, stringsAsFactors = FALSE)
clusters34$barcode <- substr(clusters34$barcode, start = 1, stop = nchar(clusters34$barcode)-2)
rownames(clusters34) <- clusters34$barcode
# Reading in the 10X data
data34 <- Read10X("34to39/filtered_feature_bc_matrix")
doublets34 <- subset(clusters34, status == "doublet")
unassigned34 <- subset(clusters34, status == "unassigned")
clusters34 <- subset(clusters34, status == "singlet")
# We remove all cells that are not singlets
colnames(data34) <- substring(colnames(data34),1,nchar(colnames(data34))-2)
data34 <- data34[, colnames(data34) %in% clusters34$barcode]
# We calculate the percentage of MT reads
MT.index34 <- grep(pattern = "^MT-", x = rownames(data34), value = FALSE)
all.sum34 <- Matrix::colSums(data34)
percent.MT34 <- Matrix::colSums(data34[MT.index34, ])/all.sum34
data34 <- data34[-MT.index34, ] # remove MT genes
# We calculate the percentage of Ribosomal genes
RG.index34 <- grep("^RP[SL][0-9]+$", x = rownames(data34), perl=TRUE, value = FALSE)
percent.RG34 <- Matrix::colSums(data34[RG.index34, ])/all.sum34
data34 <- data34[-RG.index34, ] 
# Adding the cluster assignment 
clusters34 <- clusters34[,c(1:3)]
clusters34$sample <- as.factor(clusters34$assignment)
sampMatch34 <- c("0" = "HC06", "1" = "PAT05", "2" = "PAT13", "3" = "PAT12", "4" = "HC11", "5" ="PAT14")
clusters34$sample <- as.character(revalue(clusters34$sample, sampMatch34))
cluster.assignment34 <- clusters34$sample
names(cluster.assignment34) <- clusters34$barcode
cluster.assignment34 <- cluster.assignment34[order(match(names(cluster.assignment34), names(percent.MT34)))]
# Adding the status assignment  
cluster.status34 <- clusters34$status
names(cluster.status34) <- clusters34$barcode
cluster.status34 <- cluster.status34[order(match(names(cluster.status34), names(percent.MT34)))]
# Adding the meta data
data34 <- CreateSeuratObject(data34, min.cells = 5, meta.data = data.frame(percent.mt = percent.MT34, percent.rg = percent.RG34, sample = cluster.assignment34, status = cluster.status34), project = "Pool34to39")
#VlnPlot(data34, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
data34 <- subset(data34, nFeature_RNA > 250 & nFeature_RNA < 5000 & percent.mt < 0.25)
rm(clusters34, doublets34, unassigned34, MT.index34, all.sum34, percent.MT34, percent.RG34, RG.index34, sampMatch34, cluster.assignment34, cluster.status34)

data <- merge(data9, y = c(data15, data29, data34), add.cell.ids = c("Pool9to14", "Pool15to19", "Pool29to33", "Pool34to39"))
data$pool <- data@meta.data$orig.ident
data.ori <- data
rm(data9, data15, data29, data34)

Male <- subset(data, sample %in% c("HC12","PAT10","HC02","PAT08","HC06", "PAT05","PAT13", "HC03", "PAT04"))
Male$sex <- "Male"
Female <- subset(data, sample %in% c("HC01","HC04","PAT02","PAT11","HC08","PAT03","PAT09","HC11","PAT12","PAT14","HC05","PAT07"))
Female$sex <- "Female"
data <- merge(x = Male, y = Female)
rm(Male, Female)

postcovid19 <- subset(data, sample %in% c("PAT02","PAT10","PAT11","PAT03","PAT08","PAT09","PAT05", "PAT12", "PAT13", "PAT14","PAT04", "PAT07"))
postcovid19$condition <- "postcovid19"
healthy<- subset(data, sample %in% c("HC01","HC04","HC12","HC02","HC08","HC06","HC11","HC03", "HC05"))  # delete HC07
healthy$condition <- "healthycontrol"
data <- merge(x = postcovid19 , y = healthy)
rm(postcovid19,healthy)

data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- ScaleData(data, vars.to.regress = "percent.mt")
data <- RunPCA(data, features = VariableFeatures(object = data))
ElbowPlot(data)
data <- FindNeighbors(data, dims = 1:20)
data <- FindClusters(data, resolution = 0.5)
data <- RunUMAP(data, dims = 1:20)
DimPlot(data, reduction = "umap", split.by = "condition", label = T)

saveRDS(data, file = "Rad.rds")
data <- readRDS("Rad.rds")

####### single R for annotation #############
input <- GetAssayData(object = data, slot = "data", assay = "RNA")
humap.se <- HumanPrimaryCellAtlasData()  #Obtain the Humap data
blueprint.se <- BlueprintEncodeData(rm.NA="rows") #Obtain human bulk RNA-seq data from Blueprint and ENCODE
monaco.se <- MonacoImmuneData() #Obtain bulk RNA-seq data of sorted human immune cells 
immune.se  <- DatabaseImmuneCellExpressionData() #Obtain human bulk RNA-seq data from DICE 
hemato.se  <- NovershternHematopoieticData() #sorted hematopoietic cell populations 

singleR_colors <- c("T cells" = "#776fb2",
                    "T_cells" = "#776fb2",
                    "CD4+ T-cells (naive)" = "#cecce2",
                    "CD4+ T-cells" = "#cecce2",
                    "CD4+ T cells" = "#cecce2",
                    "T cells, CD4+" = "#cecce2",
                    "CD4+/CD45RA+/CD25- Naive T" = "#cecce2",
                    "CD4+ T Helper2" = "#cecce2",
                    "CD4+ Tcm" = "#cecce2",
                    "CD4+/CD45RO+ Memory" = "#cecce2",
                    "CD4+ memory T-cells" = "#cecce2",
                    "CD4+ Tem" = "#cecce2",
                    
                    "CD8+ T-cells (naive)" = "#422483",
                    "CD8+/CD45RA+ Naive Cytotoxic" = "#422483",
                    "CD8+ T-cells" = "#422483",
                    "CD8+ T cells" = "#422483",
                    "T cells, CD8+" = "#422483",
                    "CD8+ Tcm" = "#422483",
                    "CD8+ Cytotoxic T" = "#422483",
                    "CD8+ Tem" = "#422483",
                    
                    "Treg cells" = "#004c9d",
                    "CD4+/CD25 T Reg" = "#004c9d",
                    "regulatory T-cells" = "#004c9d",
                    
                    "NKT cells" = "#684495",
                    "NK T cells" = "#684495",
                    "NK cells" = "#338eb0",
                    "NK_cell" = "#338eb0",
                    "CD56+ NK" = "#338eb0",
                    
                    "ILCs" = "#d9dada",
                    
                    "naive B-cells" = "#00963f",
                    "B-cells" = "#00963f",
                    "B_cell" = "#00963f",
                    "B cells" = "#00963f",
                    "CD19+ B" = "#00963f",
                    "Pre-B_cell_CD34-" = "#00961a" ,
                    "Pro-B_cell_CD34+" = "#00961a",
                    "memory B-cells" = "#32ab6d",
                    "class-switched memory B-cells" = "#7dc29e",
                    "Plasma cells" = "#d5e7dd",
                    
                    "BM" = "#b3a930",
                    "BM & Prog." = "#b3a930",
                    "Progenitors" = "#b3a930", 
                    "HSC" = "#b3a930",
                    "HSCs" = "#b3a930",
                    "HSC_-G-CSF" = "#b3a930",
                    "HSC_CD34+" = "#b3a930",
                    "CD34+ Precursors" = "#b3a930",
                    "MPP" = "#dfd200",
                    "CLP" = "#ffed00",
                    "CMP" = "#fdef6d",
                    "CMPs" = "#fdef6d",
                    "GMP" = "#faf3a8",
                    "GMPs" = "#faf3a8",
                    "MEP" = "#e7bd00",
                    "MEPs" = "#e7bd00",
                    "Megakaryocytes" = "#efd873",
                    
                    "DC" = "#ef7c00",
                    "Dendritic" = "#ef7c00",
                    "Dendritic cells" = "#ef7c00",
                    
                    "Monocyte (CD16-)" = "#e6330f",
                    "Monocyte (CD16+)" = "#ea5552",
                    "Monocyte (CD14+)" = "#f4a5a5",
                    "Monocytes" = "#f4a5a5",
                    "Monocyte" = "#f4a5a5",
                    "CD14+ Monocyte" = "#f4a5a5",
                    
                    "Pro-Myelocyte" = "#001816",
                    "Myelocyte" = "#00312C",
                    "Granulocytes" = "#006358",
                    "Eosinophils" = "#00af9d",
                    "Neutrophils" = "#87cbbe",
                    "Basophils" = "#cae6e4",
                    "Macrophages" = "#b41909",
                    "Macrophage" = "#b41909",
                    "Erythrocytes" = "#bb79b2",
                    "Erythroblast" = "#bb79b2",
                    "Erythroid cells" = "#bb79b2",
                    "Platelets" = "#2a3937",
                    
                    "Adipocytes" = "#e2a9cd",
                    "Fibroblasts" = "#be348b",
                    "Endothelial cells" = "#7d2685",
                    "Endothelial_cells" = "#7d2685",
                    "mv Endothelial cells" = "#632282",
                    "Myocytes"="#A70000",
                    "Smooth_muscle_cells"="#A70000",
                    "Chondrocytes"="#F0F7DA",
                    "Epithelial_cells"="#A67C00",
                    "Neurons"="#63CDE3"
)

singleR.list <- list()

singleR.list$humap <- SingleR(test = input, 
                             method="single",
                             fine.tune=FALSE,
                             ref = humap.se, 
                             labels = humap.se$label.main)

singleR.list$blueprint <- SingleR(test = input, 
                                  method="single",
                                  fine.tune=FALSE,
                                  ref = blueprint.se, 
                                  labels = blueprint.se$label.main)

singleR.list$monaco <- SingleR(test = input, 
                               method="single",
                               fine.tune=FALSE,
                               ref = monaco.se, 
                               labels = monaco.se$label.main)

singleR.list$immune <- SingleR(test = input, 
                               method="single",
                               fine.tune=FALSE,
                               ref = immune.se, 
                               labels = immune.se$label.main)

singleR.list$hemato <- SingleR(test = input, 
                               method="single",
                               fine.tune=FALSE,
                               ref = hemato.se, 
                               labels = hemato.se$label.main)

data$humap.labels <- singleR.list$humap$labels
data$blueprint.labels <- singleR.list$blueprint$labels
data$monaco.labels <- singleR.list$monaco$labels
data$immune.labels <- singleR.list$immune$labels
data$hemato.labels <- singleR.list$hemato$labels


p1 <- DimPlot(object = data, reduction = 'umap', label = FALSE, group.by ="hemato.labels") + 
  ggtitle("Hematopoietic labels")+scale_color_manual(values=singleR_colors)
p2 <- DimPlot(object = data, reduction = 'umap', label = FALSE, group.by ="humap.labels") + 
  ggtitle("Humap")+scale_color_manual(values=singleR_colors)
p3 <- DimPlot(object = data, reduction = 'umap', label = FALSE, group.by ="blueprint.labels") +
  ggtitle("Blueprint-Encode")+scale_color_manual(values=singleR_colors)
p4 <- DimPlot(object = data, reduction = 'umap', label = FALSE, group.by ="monaco.labels") + 
  ggtitle("Monaco")+scale_color_manual(values=singleR_colors)
p5 <- DimPlot(object = data, reduction = 'umap', label = FALSE, group.by ="immune.labels") + 
  ggtitle("Immune")+scale_color_manual(values=singleR_colors)

pdf("plot_singleR.rad.pdf",width=20, height=10)
CombinePlots(plots = list(p1, p2, p3, p4, p5), ncol=3)
dev.off()

cell.type <- data@meta.data$seurat_clusters
levels(cell.type) <- as.factor(c("Classical Mono", "CD4", "NK", "CD4", "CD4", "CD8", "Non_classical Mono", "B", "Classical Mono", "B", "CD8", "mDC","Mega","Classical Mono", "pDC","Classical Mono", "Prol.T", "undefined"))
data@meta.data$celltype <- cell.type
data$celltype <- as.factor(data$celltype)
data@meta.data$BroadCelltype <- cell.type

pdf("umap_Rad.pdf", width = 6.5, height = 5)
DimPlot(data, group.by = "condition")
dev.off()

DimPlot(data, group.by = "celltype", pt.size = 0.5, reduction = "umap", label = T) +
        ggtitle(label = "UMAP of cell types") +
        theme(text = element_text(size = 30),
              axis.text = element_text(size = 20))

## add phenotype information 
data@meta.data$age2 <- data@meta.data$sample
data@meta.data$age2 <- gsub("PAT10", 53, data@meta.data$age)
data@meta.data$age2 <- gsub("PAT02", 69, data@meta.data$age)
data@meta.data$age2 <- gsub("PAT11", 47, data@meta.data$age)
data@meta.data$age2 <- gsub("PAT08", 67, data@meta.data$age)
data@meta.data$age2 <- gsub("PAT09", 56, data@meta.data$age)
data@meta.data$age2 <- gsub("PAT03", 71, data@meta.data$age)
data@meta.data$age2 <- gsub("PAT04", 66, data@meta.data$age)
data@meta.data$age2 <- gsub("PAT07", 58, data@meta.data$age)
data@meta.data$age2 <- gsub("PAT13", 58, data@meta.data$age)
data@meta.data$age2 <- gsub("PAT12", 30, data@meta.data$age)
data@meta.data$age2 <- gsub("PAT05", 68, data@meta.data$age)
data@meta.data$age2 <- gsub("PAT14", 63, data@meta.data$age)
data@meta.data$age2 <- gsub("HC04", 55, data@meta.data$age)
data@meta.data$age2 <- gsub("HC01", 64, data@meta.data$age)
data@meta.data$age2 <- gsub("HC12", 61, data@meta.data$age)
data@meta.data$age2 <- gsub("HC02", 63, data@meta.data$age)
data@meta.data$age2 <- gsub("HC08", 39, data@meta.data$age)
data@meta.data$age2 <- gsub("HC05", 63, data@meta.data$age)
data@meta.data$age2 <- gsub("HC03", 51, data@meta.data$age)
data@meta.data$age2 <- gsub("HC06", 60, data@meta.data$age)
data@meta.data$age2 <- gsub("HC11", 60, data@meta.data$age)

data@meta.data$gender2 <- data@meta.data$sample
data@meta.data$gender2 <- gsub("PAT10", 1, data@meta.data$gender)
data@meta.data$gender2 <- gsub("PAT02", 0, data@meta.data$gender)
data@meta.data$gender2 <- gsub("PAT11", 0, data@meta.data$gender)
data@meta.data$gender2 <- gsub("PAT08", 1, data@meta.data$gender)
data@meta.data$gender2 <- gsub("PAT09", 0, data@meta.data$gender)
data@meta.data$gender2 <- gsub("PAT03", 0, data@meta.data$gender)
data@meta.data$gender2 <- gsub("PAT04", 1, data@meta.data$gender)
data@meta.data$gender2 <- gsub("PAT07", 0, data@meta.data$gender)
data@meta.data$gender2 <- gsub("PAT13", 1, data@meta.data$gender)
data@meta.data$gender2 <- gsub("PAT12", 0, data@meta.data$gender)
data@meta.data$gender2 <- gsub("PAT05", 1, data@meta.data$gender)
data@meta.data$gender2 <- gsub("PAT14", 0, data@meta.data$gender)
data@meta.data$gender2 <- gsub("HC04", 0, data@meta.data$gender)
data@meta.data$gender2 <- gsub("HC01", 0, data@meta.data$gender)
data@meta.data$gender2 <- gsub("HC12", 1, data@meta.data$gender)
data@meta.data$gender2 <- gsub("HC02", 1, data@meta.data$gender)
data@meta.data$gender2 <- gsub("HC08", 0, data@meta.data$gender)
data@meta.data$gender2 <- gsub("HC05", 0, data@meta.data$gender)
data@meta.data$gender2 <- gsub("HC03", 1, data@meta.data$gender)
data@meta.data$gender2 <- gsub("HC06", 1, data@meta.data$gender)
data@meta.data$gender2 <- gsub("HC11", 0, data@meta.data$gender)

data@meta.data$gender2 = as.numeric(data@meta.data$gender2)
data@meta.data$age2 = as.numeric(data@meta.data$age2)

unique(data@meta.data$celltype)
data$condition2 <- data$condition
data$condition2 <- gsub("postcovid19", "Post-COVID-19", data$condition2)
data$condition2 <- gsub("healthycontrol", "Control", data$condition2)
data_copy <- data
my_levels <- c("CD4","CD8","B","NK","CD14+ Mono","CD16+ Mono","mDC","pDC","Mega","Prol.T","undefined")
data_copy@meta.data$celltype <- factor(data_copy$celltype, levels= my_levels)

cell.type2 <- data_copy@meta.data$celltype
levels(cell.type2) <- as.factor(c("CD4","CD8","B","NK","Classical Mono","Non-classical Mono","mDC","pDC","Mega","Prol.T","undefined"))
data_copy@meta.data$celltype2 <- cell.type2

## figure 1A 
DimPlot(data_copy, group.by = "celltype2", split.by = "condition2")
## figure 1C ####
identity.to.compare <- "condition2"
condition.to.compare <- as.character(levels(as.factor(data@meta.data$condition2)))

celltypes <- as.character(levels(data@meta.data$BroadCelltype))
type <- "BroadCelltype"

result <- list()
DescriptionCelltype <- function(SeuratObject, identity.to.compare, condition.to.compare, celltypes){
  
  # Result object
  result <- list()
  
  # Subsetting to the relevant subgroup
  data <- FetchData(SeuratObject, vars = identity.to.compare)
  data.subgroup <- SeuratObject[, which(x = data == condition.to.compare)]
  
  # Calculating the percentages per cell type per group
  result$CellsPerType <- table(as.character(data.subgroup$celltype))
  result$PercentagesPerType <- result$CellsPerType/dim(data.subgroup)[2]
  
  # Calculating the percentages per cell type per patient
  patients <- unique(data.subgroup$sample)
  patients <- patients[order(patients)]
  print(patients)
  CellsPerTypePerPatient <- c()
  PercentagesPerTypePerPatient <- c()
  values.for.MeanVariance <- c()
  for(person in patients){
    # Not all patients have all cell types at all time points.
    # With this dummy, we can use 0 for missing values
    result.dummy <- rep(0,length(celltypes))
    names(result.dummy) <- celltypes
    result.dummy <- result.dummy[sort(names(result.dummy))]
    
    # We subset to a sample ID (one person at one time point)
    data.person <- FetchData(data.subgroup, vars = "sample")
    data.person <- data.subgroup[, which(x = data.person == person)]
    
    # How many cells does a patient have of each type?
    result.preliminary <- table(as.character(data.person[[type]][,1]))
    result.preliminary <- result.preliminary[sort(names(result.preliminary))]
    
    result.dummy[names(result.dummy) %in% names(result.preliminary)] <- result.preliminary
    CellsPerTypePerPatient <- rbind(CellsPerTypePerPatient, result.dummy)
    rownames(CellsPerTypePerPatient)[dim(CellsPerTypePerPatient)[1]] <- paste0("Patient_", person)
    
    
    # Percentages
    result.dummy <- result.dummy/sum(result.dummy)
    PercentagesPerTypePerPatient <- rbind(PercentagesPerTypePerPatient, result.dummy)
    rownames(PercentagesPerTypePerPatient)[dim(PercentagesPerTypePerPatient)[1]] <- paste0("Patient_", person)
  }
  
  result$CellsPerPatient <- CellsPerTypePerPatient
  result$PercentagesPerPatient <- PercentagesPerTypePerPatient
  
  # Calculating the Mean, Variance and Coefficient of Variance per cell type
  result$MeanPerCelltype <- colMeans(PercentagesPerTypePerPatient)
  result$VarsPerCelltype <- apply(PercentagesPerTypePerPatient, 2, var)
  result$CoeffOfVariance <- result$VarsPerCelltype / result$MeanPerCelltype
  result$Condition <- condition.to.compare
  return(result)
}

CellTypePerCondition <- lapply(condition.to.compare, DescriptionCelltype, 
                               SeuratObject = data, 
                               identity.to.compare = identity.to.compare,
                               celltypes = celltypes)
data.boxplot.combined <- c()
for(conditions in 1:length(condition.to.compare)){
  data.boxplot <- data.frame(patient = rownames(CellTypePerCondition[[conditions]]$PercentagesPerPatient), CellTypePerCondition[[conditions]]$PercentagesPerPatient, stringsAsFactors = FALSE)
  #data.boxplot <- melt(data.boxplot, id.vars = "patient")
  
  #pdf(paste0(output_path, "Condition_", condition.to.compare[conditions], "_CellTypes.pdf"))
  #print(
  # ggplot(data.boxplot, aes(x = variable, y = value)) +
  #    geom_boxplot() +
  #    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  # dev.off()
  
  # Combining the data to plot all conditions in one boxplot
  data.boxplot <- data.frame(data.boxplot, condition = condition.to.compare[conditions])
  data.boxplot.combined <- rbind(data.boxplot.combined, data.boxplot)
}

data.boxplot.combined$name <- data.boxplot.combined$patient
data.boxplot.combined <- data.boxplot.combined[,-1]
data.boxplot.combined <- data.table::melt(data.boxplot.combined, id = 12:13, measure = colnames(data.boxplot.combined)[-c(12,13)])
data.boxplot.combined <- subset(data.boxplot.combined, variable %in% c("B","CD14..Mono","CD16..Mono","CD4","CD8","mDC","NK"))

data.boxplot.combined$variable <- factor(data.boxplot.combined$variable, levels = c("CD4","CD8","B","NK","CD14..Mono","CD16..Mono", "mDC"))

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

area <- ggplot(data.boxplot.combined, aes(x = variable, y = value)) +
  geom_boxplot(aes(fill = condition)) +
  xlab("cell type") +
  ylab("Cell proportion based on sc-RNA") +
  theme_classic() +
  # geom_hline(yintercept = 0, shape = group) +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        legend.position = c(0.9,0.9))+
  scale_x_discrete(breaks=c("CD4","CD8","B","NK","CD14..Mono","CD16..Mono", "mDC"), labels = addline_format(c("CD4","CD8","B","NK","Classical Mono","Non-classical Mono", "mDC")))
area + scale_fill_manual(values = c("steelblue", "tomato")) + geom_jitter(aes(group=factor(condition)),position=position_dodge(width=0.8))

## DEG analysis ##
significance <- function(x, treshold = 0.25, pval = 0.05){
  x$sig <- "nonsig"
  x[x$p_val_adj <= pval, "sig"] <- "sig"
  x[abs(x$avg_logFC) > treshold & x$p_val_adj <= pval, "sig"] <- "very sig"
  x <- data.frame(gene = rownames(x), x, stringsAsFactors = FALSE)
  return(x)
}  

# Volcano Plot Function
volcano <- function(x, title, top = 10){
  x <- x[order(abs(x$p_val_adj), decreasing = FALSE),]
  volcano.labels <- ifelse(x$sig %in% c("sig","very sig"), as.character(x$gene), "")
  volcano.labels[(top+1):length(volcano.labels)] <- ""
  
  ggplot(x, aes(x = avg_logFC, y = -log(p_val_adj, 10), label = gene)) +
    geom_point(aes(color = sig), size = 6) +
    geom_label_repel(aes(label = volcano.labels), 
                     hjust = 1.25, vjust = 0, size = 8) +
    scale_colour_manual(name = "Significance", 
                        values = c("very sig" = "red", "sig" = "black", "nonsig" = "gray"),
                        breaks = c("very sig", "sig", "nonsig"),
                        labels = c("very sig", "sig", "nonsig")) +
    ggtitle(label = title) +
    ylab("-log10(p_val_adj)") +
    xlab("avg_logFC") +
    theme_bw()+
    geom_hline(yintercept = -log(0.05,10)) +
    geom_vline(xintercept = c(-0.25, 0.25)) +
    theme(axis.text = element_text(size = 30),
          axis.title = element_text(size = 30),
          title = element_text(size = 30),
          legend.text = element_text(size = 30))
}

Idents(data) <- "condition"
DEG.conditions <- FindMarkers(data, logfc.threshold = 0, ident.1 = ident.1 = "postcovid19", ident.2 = "healthycontrol")
volcano.data <- significance(DEG.conditions)

p1 <- volcano(x = volcano.data, title = "DEG between disease and control")
png(paste0(OUTPUT_DIR, "pipeline/DEG/Volcano_DEG_bulk.png"), width = 1920, height = 800)
print(volcano(x = volcano.data, title = "DEG between disease and control"))
dev.off()

info <- matrix(0, ncol = 3, nrow = length(unique(data@meta.data$BroadCelltype))+1)
colnames(info) <- c("celltype", "sig", "verysig")
rownames(info) <- c("condition", levels(data@meta.data$BroadCelltype))
info["condition",] <- c("condition", sum(volcano.data$sig == "sig"), sum(volcano.data$sig == "very sig"))

volcano.ls <- list()
volcano.data.ls <- list()

Idents(data) <- "condition"
for(celltype in unique(data$BroadCelltype)){
  data.celltype <- FetchData(data, vars = "BroadCelltype")
  data.celltype <- data[, which(x = data.celltype == celltype)]
  
  DEG.celltype <- FindMarkers(data.celltype, logfc.threshold = 0, ident.1 = "postcovid19", ident.2 = "healthycontrol", min.cells.feature = 1, min.cells.group = 1,test.use = "MAST", latent.vars = c("age", "gender"))
  volcano.data <- significance(DEG.celltype)
  
  # Printing the volcano plot for the DEG per cell type
  p1 <- volcano(x = volcano.data, title = paste0("DEG between disease and control for ", celltype))
  volcano.ls[[celltype]] <- p1
  png(paste0("Rad2/Rad_age_gender/Broad_Volcano_DEG_", celltype, ".png"), width = 1920, height = 800)
  print(p1)
  dev.off()
  
  # We write the DEG in a csv
  volcano.data <- subset(volcano.data, sig %in% c("sig", "very sig"))
  info[celltype,] <- c(celltype, sum(volcano.data$sig == "sig"), sum(volcano.data$sig == "very sig"))
  
  write.csv(volcano.data, paste0("Rad2/Rad_age_gender/DEG_", celltype, ".csv"), row.names = FALSE)
  volcano.data <- subset(volcano.data, sig %in% c("very sig"))
  write.csv(volcano.data, paste0("Rad2/Rad_age_gender/DEG_", celltype, "_verysig.csv"), row.names = FALSE)
  volcano.data.ls[[celltype]] <- volcano.data
}  

write.csv(info, "Rad2/Rad_age_gender/NumberOfDEG.csv", quote = FALSE)

## pathway analysis
B <- read.csv("/Users/zli20/Documents/Post-COVID-19/Post_covid_scRNA/1to4/Rad2/Rad_age_gender/DEG_B.csv")
Mono16 <- read.csv("/Users/zli20/Documents/Post-COVID-19/Post_covid_scRNA/1to4/Rad2/Rad_age_gender/DEG_CD16+ Mono.csv")
Mono14 <- read.csv("/Users/zli20/Documents/Post-COVID-19/Post_covid_scRNA/1to4/Rad2/Rad_age_gender/DEG_CD14+ Mono.csv")
NK <- read.csv("/Users/zli20/Documents/Post-COVID-19/Post_covid_scRNA/1to4/Rad2/Rad_age_gender/DEG_NK.csv")
CD8 <- read.csv("/Users/zli20/Documents/Post-COVID-19/Post_covid_scRNA/1to4/Rad2/Rad_age_gender/DEG_CD8.csv")
CD4 <- read.csv("/Users/zli20/Documents/Post-COVID-19/Post_covid_scRNA/1to4/Rad2/Rad_age_gender/DEG_CD4.csv")
mDC <- read.csv("/Users/zli20/Documents/Post-COVID-19/Post_covid_scRNA/1to4/Rad2/Rad_age_gender/DEG_mDC.csv")

B$cell <- "B"
Mono16$cell <- "Mono16"
Mono14$cell <- "Mono14"
NK$cell <- "NK"
CD8$cell <- "CD8"
CD4$cell <- "CD4"
mDC$cell <- "mDC"

Mono16.up <- subset(Mono16, avg_logFC > 0.05)
Mono16.down <- subset(Mono16, avg_logFC < -0.05)

Mono16.down <- AnnotationDbi::select(org.Hs.eg.db, keys=Mono16.down$gene, columns=c("ENTREZID","SYMBOL"), keytype="SYMBOL")
Mono16.GO.down <- enrichGO(gene = Mono16.down$ENTREZID,
                           OrgDb = org.Hs.eg.db,
                           keyType = "ENTREZID",
                           ont = "ALL",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.01,
                           qvalueCutoff = 0.05,
                           readable      = TRUE)
Mono16_down_bp <- Mono16.GO.down[which(Mono16.GO.down$ONTOLOGY == "BP"),]
Mono16_down_mf <- Mono16.GO.down[which(Mono16.GO.down$ONTOLOGY == "MF"),]
Mono16.up <- AnnotationDbi::select(org.Hs.eg.db, keys=Mono16.up$gene, columns=c("ENTREZID","SYMBOL"), keytype="SYMBOL")
Mono16.GO.up <- enrichGO(gene = Mono16.up$ENTREZID,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID",
                         ont = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.01,
                         qvalueCutoff = 0.05,
                         readable      = TRUE)
Mono16_up_bp <- Mono16.GO.up[which(Mono16.GO.up$ONTOLOGY == "BP"),]
Mono16_up_mf <- Mono16.GO.up[which(Mono16.GO.up$ONTOLOGY == "MF"),]
Mono16_up <- Mono16.GO.up[which(Mono16.GO.up$ONTOLOGY == c("MF","BP")),]
merge_result(list(Non_classical_Mono_up_BP = Mono16_up_bp, Non_classical_Mono_up__MF = Mono16_up_mf,
                  Non_classical_Mono_down_BP = Mono16_down_bp, Non_classical_Mono_down_MF = Mono16_down_mf)) %>%
  clusterProfiler::dotplot(.,showCategory=10)+
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 50) )


## figure 3A ##
MHC2 <- c(
  'HLA-DRB5',
  'HLA-DQB1',
  'HLA-DPB1',
  'HLA-DMA',
  'HLA-DPA1',
  'HLA-DQA1',
  'HLA-DMB',
  'HLA-DRB1',
  'HLA-DRA'
)

cd14 <- subset(data, celltype == "CD14+ Mono")
cd16 <- subset(data, celltype == "CD16+ Mono")
Idents(cd14) <- "condition"
Idents(cd16) <- "condition"
p1 <- DotPlot(cd14, features = MHC2, scale.by = "radius")+ggtitle("Classical monocytes")
p2 <- DotPlot(cd16, features = MHC2, scale.by = "radius")+ggtitle("Non-classical monocytes")
plot_grid(p1, p2, ncol=1)

#### figure 3C and 3D #####
mono <- subset(data, celltype == "CD14+ Mono")
mono <- NormalizeData(mono)
mono <- FindVariableFeatures(mono, selection.method = "vst", nfeatures = 2000)
mono <- ScaleData(mono)
mono <- RunPCA(mono, features = VariableFeatures(object = mono))
mono <- FindNeighbors(mono, dims = 1:10)
mono <- FindClusters(mono, resolution = 0.2)
mono <- RunUMAP(mono, dims = 1:10)

cell.type <- data@meta.data$seurat_clusters
levels(cell.type) <- as.factor(c("c0","c1","c2","c3","c4","c5"))
data@meta.data$celltype3 <- cell.type

DimPlot(data, reduction = "umap", split.by = "condition2", group.by = "celltype3")

### supplementary figure 4A 
B <- read.csv("../Post-COVID-19/Post_covid_scRNA/1to4/Rad2/Rad_age_gender/DEG_B.csv")
Mono16 <- read.csv("../Post-COVID-19/Post_covid_scRNA/1to4/Rad2/Rad_age_gender/DEG_CD16+ Mono.csv")
Mono14 <- read.csv("../Post-COVID-19/Post_covid_scRNA/1to4/Rad2/Rad_age_gender/DEG_CD14+ Mono.csv")
NK <- read.csv("../Post-COVID-19/Post_covid_scRNA/1to4/Rad2/Rad_age_gender/DEG_NK.csv")
CD8 <- read.csv("../Post-COVID-19/Post_covid_scRNA/1to4/Rad2/Rad_age_gender/DEG_CD8.csv")
CD4 <- read.csv("../Post-COVID-19/Post_covid_scRNA/1to4/Rad2/Rad_age_gender/DEG_CD4.csv")
mDC <- read.csv("../Post-COVID-19/Post_covid_scRNA/1to4/Rad2/Rad_age_gender/DEG_mDC.csv")

B.up <- subset(B, avg_logFC > 0.05) #13
B.down <- subset(B, avg_logFC < -0.05) # 17
Mono16.up <- subset(Mono16, avg_logFC > 0.05) # 35
Mono16.down <- subset(Mono16, avg_logFC < -0.05) # 27
Mono14.up <- subset(Mono14, avg_logFC > 0.05) # 217
Mono14.down <- subset(Mono14, avg_logFC < -0.05) # 362
NK.up <- subset(NK, avg_logFC > 0.05) # 11
NK.down <- subset(NK, avg_logFC < -0.05) # 17
CD8.up <- subset(CD8, avg_logFC > 0.05) # 10
CD8.down <- subset(CD8, avg_logFC < -0.05) # 16
CD4.up <- subset(CD4, avg_logFC > 0.05) # 77
CD4.down <- subset(CD4, avg_logFC < -0.05) # 95
mDC.up <- subset(mDC, avg_logFC > 0.05) # 2
mDC.down <- subset(mDC, avg_logFC < -0.05) # 8

x <- rep(c("CD4","CD8","B","NK","Classical Mono","Non-classical Mono","mDC"), each = 2)
y <- rep(c("up","down"), times= 7)
z <- c(77,95,10,16,13,17,11,17,217,362,35,27,2,8)
#z <- c(29,21,42,33,53,47,83,67)
df1 <- data.frame(x=x, y=y,z=z)
df1$x = factor(x,levels = c("CD4","CD8","B","NK","Classical Mono","Non-classical Mono","mDC"))

ggplot(data = df1, mapping = aes(x = x, y = z, fill = y))+geom_bar(stat = "identity", position = "dodge")+scale_x_discrete(breaks = c("CD4","CD8","B","NK","Classical Mono","Non-classical Mono","mDC"), labels = addline_format(c("CD4","CD8","B","NK","Classical Mono","Non-classical Mono","mDC")))+scale_fill_manual(values = c("steelblue","tomato"),name="Post-COVID-19 vs Control")+theme_classic()+xlab("Celltype")+ylab("Number of DEGs")+
  ylim(0,400)+
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        legend.position = c(0.85,0.9))


## delete HC04 and recluster
data2 <- subset(data, sample != "HCO4")
data <- data2
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- ScaleData(data, vars.to.regress = "percent.mt")
data <- RunPCA(data, features = VariableFeatures(object = data))
ElbowPlot(data)

data <- FindNeighbors(data, dims = 1:20)
data <- FindClusters(data, resolution = 0.5)
data <- RunUMAP(data, dims = 1:20)
DimPlot(data, reduction = "umap", split.by = "condition", label = T)


res <- list(DEG.c0 = DEG.c0, 
            DEG.c1 = DEG.c1, 
            DEG.c2 = DEG.c2,
            DEG.c3 = DEG.c3,
            DEG.c4 = DEG.c4)

write.xlsx(x = res, file = 'mono/DE.xlsx')

###figure 4F###
merge_result(list(c0_down = c0_down_bp, c2_down = c2_down_bp)) %>%
  clusterProfiler::dotplot(.,showCategory=10)


merge_result(list(c0_up = c0_up_bp, c1_up= c1_up_bp)) %>%
  clusterProfiler::dotplot(.,showCategory=10)

###############################################
