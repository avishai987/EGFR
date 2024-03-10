#library("Seurat")
#install.packages("Seurat", lib="/vol/sci/bio/data/yotam.drier/Rlib/3.5")
library(Seurat, lib.loc =  "/vol/sci/bio/data/yotam.drier/Rlib/3.5")
library("dplyr")
library("ggplot2")


#lung <- readRDS("D:/yotamd/drier_lab/projects/lung_ss/raz/lung_seurat")



tpm <- function(counts, lengths) {
  
  rpk <- counts / lengths
  coef <- sum(rpk) / 1e6
  rpk/coef
}

setwd("/vol/sci/bio/data/yotam.drier/raz.sher/yotam/")
lung.data <- read.delim("fc.txt",skip=0,header=T,sep="\t",stringsAsFactors=F)
rownames(lung.data)<- make.unique(lung.data$gene_name)
lengths <- lung.data[,6] 
omitgenes <- startsWith(rownames(lung.data),"MT-")|startsWith(rownames(lung.data),"ERCC-")
lung.data <- lung.data[,8:ncol(lung.data)]
# lung.tpm <- read.table("lung_tpm.txt") 

# changing the cell labels
cell.labels <- gsub(".gene_counts.tsv","",colnames(lung.data))
cell.labels <- gsub(".sort.bam","",cell.labels)
cell.labels <- gsub("_",".",cell.labels)
cell.labels <- gsub(".S[0-9]*","",cell.labels)
t <- regexpr("[A-H][0-9]+$", cell.labels)
cell.labels[t<0] <- paste0(regmatches(cell.labels[t<0], regexpr("[^ACGT]*", cell.labels[t<0])),regmatches(colnames(lung.data)[t<0], regexpr("_S[0-9]+", colnames(lung.data)[t<0])))
# cell.labels[t<0]
well <- regmatches(cell.labels, regexpr("[A-H,S][0-9]+$", cell.labels))
plate <- substr(cell.labels,1,nchar(cell.labels)-nchar(well)-1)
# length(cell.labels)
# length(well)
# length(plate)
colnames(lung.data) <- paste(plate,well,sep="_")
lung.tpm <- apply(lung.data[!omitgenes,1:dim(lung.data)[2]], 2, function(x) tpm(x, lengths[!omitgenes]) )
colnames(lung.tpm) <- paste(plate,well,sep="_")
write.table(as.matrix(lung.tpm), file = paste("lung_tpm.txt", sep = ""), sep = "\t", row.names = TRUE)


# for using the updated Seurat library in luster


# converting the raw data to seurat object in luster
#lung.data <- read.delim("/mnt/lustre/hms-01/fs01/yotam.drier/raz.sher/140719/new_fc.txt",skip=0,header=T,sep="\t",stringsAsFactors=F)
lung <- CreateSeuratObject(counts = lung.data, project = "lung", min.cells = 3, min.features = 200)
save(lung, file = "lung_seurat.Robj")
lung.tpm <- read.delim("lung_tpm.txt",skip=0,header=T,sep="\t",stringsAsFactors=F)
lung_tpm <- CreateSeuratObject(counts = lung.tpm, project = "lung_tpm", min.cells = 3, min.features = 200)
save(lung_tpm, file = "lung_tpm_seurat.Robj")
prefix="all"

dim(lung)

id <- FetchData(lung, "orig.ident") # gives the name of the sample for each cell 
cell_number <- as.data.frame(table(id)) # creates df of number of cells sequenced for each sample

percent.mt <- PercentageFeatureSet(lung, pattern = "^MT-")
lung[["percent.mt"]] <- percent.mt

nFeature_RNA_plot <- VlnPlot(lung , features = "nFeature_RNA") + theme(legend.position="none", axis.text.x = element_text(size=8))
nCount_RNA_plot <- VlnPlot(lung , features = "nCount_RNA") + theme(legend.position="none", axis.text.x = element_text(size=8))
percent.mt_plot <- VlnPlot(lung , features = "percent.mt") + theme(legend.position="none", axis.text.x = element_text(size=8))

plot1 <- FeatureScatter(lung, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  theme(legend.position="none", axis.text.x = element_text(size=8)) + 
  geom_point(color='darkblue')
plot2 <- FeatureScatter(lung, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  theme(legend.position="none", axis.text.x = element_text(size=8)) + 
  geom_point(color='darkblue')
plot3 <- FeatureScatter(lung, feature1 = "nFeature_RNA", feature2 = "percent.mt") + 
  theme(legend.position="none", axis.text.x = element_text(size=8)) + 
  geom_point(color='darkblue')

pdf("qc.pdf")
plot3
plot1
plot2
dev.off()

lung <- CreateSeuratObject(counts = log2(lung.tpm+1), project = "lung", min.cells = 3, min.features = 200)
lung[["percent.mt"]] <- percent.mt

# subseting the data according to the QC
lung <- subset(lung, subset = nFeature_RNA > 2000 & percent.mt < 35)
dim(lung)
# normalize the data
#lung <- NormalizeData(lung)

# calculate how many cells in each sample after the QC
idQC <- FetchData(lung, "orig.ident")
cell_numberQC <- as.data.frame(table(idQC))
cell_number[["PassedQC"]] <- cell_numberQC$Freq
cell_number[["percentagePQC"]] <- cell_number$PassedQC / cell_number$Freq
write.csv(cell_number, file = "LungCellQuality.csv")

# Identification of highly variable features
lung <- FindVariableFeatures(lung, selection.method = "vst", nfeatures = 2000)

# Scaling the data
all.genes <- rownames(lung)
lung <- ScaleData(lung, features = all.genes)

#lung <- readRDS("lung_cancercells")
#prefix <- "cancercells"

# Perform linear dimensional reduction (PCA)
lung <- RunPCA(lung, features = VariableFeatures(object = lung))

pcaplot <- PCAPlot(lung, label = FALSE)

elbowplot <- ElbowPlot(lung, ndims = 50) # checking the dimensionality 

# cluster the cells
lung <- FindNeighbors(lung, dims = 1:25)
lung <- FindClusters(lung, resolution = .5)

# Run non-linear dimensional reduction (UMAP)
lung <- RunUMAP(lung, dims = 1:25)

lung_umap <- DimPlot(object = lung, reduction = "umap", pt.size = 0.1, label = TRUE)


pdf(paste(prefix,"Lung_UMAP_Plot.pdf",sep="_"))
pcaplot
elbowplot
lung_umap
dev.off()

#clusters_markers <- FindAllMarkers(lung, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#clusters_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)


pdf(paste(prefix,"white blood cells markers.pdf",sep="_"))
FeaturePlot(lung, 
            features =  c("CD68","PTPRC","CD163","CD14","CD53","CD68","CD52","C3AR1","LAPTM5","IGSF6","SRGN"), 
            pt.size = 0.1, order = TRUE)
dev.off()

pdf(paste(prefix,"endothelial cells.pdf",sep="_"))
FeaturePlot(lung, features =  c("VWF","ERG","PTPRB","TIE1"), pt.size = 0.1, order = TRUE) # endothelial cells
dev.off()

pdf(paste(prefix,"lung cells markers.pdf",sep="_"))
FeaturePlot(lung, features =  c("EGFR", "MET", "EPCAM", "NKX2-1", "NAPSA"), pt.size = 0.1, order = TRUE) # lung cells
dev.off()

pdf(paste(prefix,"platelet.pdf",sep="_"))
FeaturePlot(lung, features =  c("GP9","ITGB3","GP1BA","GP1BB","GP5"), pt.size = 0.1, order = TRUE)
dev.off()

pdf(paste(prefix,"mesothelial.pdf",sep="_"))
FeaturePlot(lung, features =  c("CALB2","UPK3B","UPK1B","HAS1","WT1","SAA2","SAA1"), pt.size = 0.1, order = TRUE)
dev.off()

pdf(paste(prefix,"markers.pdf",sep="_"))
# cancer cells
VlnPlot(lung, features = c("EGFR"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("MET"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("EPCAM"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("NKX2-1"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("NAPSA"), slot = "counts", log = TRUE)

# white blood cells "CD68","PTPRC","CD163","CD14","CD53","CD52","C3AR1","LAPTM5","IGSF6","SRGN"
VlnPlot(lung, features = c("PTPRC"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("CD163"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("CD14"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("CD53"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("CD52"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("C3AR1"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("LAPTM5"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("IGSF6"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("SRGN"), slot = "counts", log = TRUE)

# endothelial cells
VlnPlot(lung, features = c("VWF"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("ERG"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("PTPRB"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("TIE1"), slot = "counts", log = TRUE)

# platlet cells "GP9","ITGB3","GP1BA","GP1BB","GP5"
VlnPlot(lung, features = c("ITGB3"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("GP1BA"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("GP1BB"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("GP5"), slot = "counts", log = TRUE)

# mesothelial cells "CALB2","UPK3B","UPK1B","HAS1","WT1","SAA2","SAA1"
VlnPlot(lung, features = c("CALB2"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("UPK3B"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("UPK1B"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("HAS1"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("WT1"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("SAA2"), slot = "counts", log = TRUE)
VlnPlot(lung, features = c("SAA1"), slot = "counts", log = TRUE)
dev.off()



# adding patient symbol to the seurat meta.data
unique(x = lung@meta.data$orig.ident)
patient.ident <- gsub("[[:punct:]]\\w*$","", lung@meta.data$orig.ident)
patient.ident <- gsub("[[:punct:]]\\w*$","", patient.ident)
patient.ident <- gsub("[[:punct:]]\\w*$","", patient.ident)
patient.ident <- gsub("MGH1068","X1068", patient.ident)
patient.ident <- gsub("MGH1144","X1144", patient.ident)
patient.ident <- gsub("MGH1155","X1155", patient.ident)
patient.ident <- gsub("MGH1071","X1071", patient.ident)
patient.ident <- gsub("MGH1167","X1167", patient.ident)

levels(as.factor(patient.ident))
lung <- AddMetaData(object = lung, metadata = as.factor(patient.ident), col.name = "patient.ident")

# umap plot, color by patient
umap_patients <- DimPlot(lung, reduction = "umap", group.by = "patient.ident",
                         label = TRUE, pt.size = 0.1)

pdf(paste(prefix,"UmapByPatient.pdf",sep="_"))
umap_patients
dev.off()

# making table of the samples cells disributaion between the clusters
for (i in 1:length(levels(lung@active.ident))) {
  cluster_of_sample_v <- as.vector(table(lung@meta.data$orig.ident[lung@meta.data$seurat_clusters == i -1]))
  if (i == 1) {
    cluster_of_samples_df <- as.data.frame(cluster_of_sample_v)
  }
  else {
    cluster_of_samples_df <- cbind(cluster_of_samples_df, cluster_of_sample_v)
  }
}

colnames(cluster_of_samples_df) <- c(0:(length(levels(lung@active.ident))-1))
row.names(cluster_of_samples_df) <- levels(lung@meta.data$orig.ident)

# making table of the patient's cells disributaion between the clusters
for (i in 1:length(levels(lung@active.ident))) {
  cluster_of_patient_v <- as.vector(table(lung@meta.data$patient.ident[lung@meta.data$seurat_clusters == i-1]))
  if (i == 1) {
    cluster_of_patient_df <- as.data.frame(cluster_of_patient_v)
  }
  else {
    cluster_of_patient_df <- cbind(cluster_of_patient_df, cluster_of_patient_v)
  }
}

colnames(cluster_of_patient_df) <- c(0:(length(levels(lung@active.ident))-1))
row.names(cluster_of_patient_df) <- levels(lung@meta.data$patient.ident)


# assigning the type of the cell of the clusters
new.cluster.ids <- c(" 0- cancer",
                     " 1- cancer",
                     " 2- cancer",
                     " 3- WBC",
                     " 4- cancer",
                     " 5- cancer",
                     " 6- mesothelial/endothelial",
                     " 7- cancer",
                     " 8- cancer",
                     " 9- cancer",
                     "10- WBC",
                     "11- cancer",
                     "12- cancer",
                     "13- cancer",
                     "14- cancer")                
names(new.cluster.ids) <- levels(lung)
lung <- RenameIdents(lung, new.cluster.ids)
umap_labled <- DimPlot(lung, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()

pdf("LungLabledUmap.pdf")
umap_labled
dev.off()

#  How many cells of each sample are called as cancer cells
# cancer clusters = 0,1,2,3,4,6,7,8,12,13,15,16,18,19
# assigning the type of the cell of the clusters
cluster_types <- substring(new.cluster.ids,5)
cell.type <- lung@meta.data$seurat_clusters
levels(cell.type) <- cluster_types
lung <- AddMetaData(object = lung, metadata = as.factor(cell.type), col.name = "cell.type")

cell_types <- FetchData(lung, "cell.type")
cell_types_sum <- as.data.frame(table(cell_types))

cell_types_table <- as.data.frame.matrix(table(lung@meta.data$orig.ident, lung@meta.data$cell.type))


samples_summary <- cbind(cell_number, "No. of cancer cells" = cell_types_table$cancer)
CancerPer <- samples_summary$`No. of cancer cells` / samples_summary$PassedQC
samples_summary <- cbind(samples_summary, "Percentage of Cancer cells" = CancerPer)
write.csv(samples_summary, file = "PlateQuality.csv")

# new samples- 1055new, X760, X793, X802, X1099*, X1144*, X1155-9_p4, X1167-1_p2

#new_samples <- grepl(pattern = "^X1055new|X760|X793|X802|X1099|X1144|X1155|X1167", x = as.character(cell_number$id))

# write.csv(samples_summary[new_samples,], file = "C:\\Users\\razsh\\OneDrive\\Documents\\lab\\280719\\LungSamplesQuality_NewSamples.csv")
# write.csv(cluster_of_patient_df, file = "C:\\Users\\razsh\\OneDrive\\Documents\\lab\\280719\\Clusters Of Patients.csv")
# write.csv(cluster_of_samples_df, file = "C:\\Users\\razsh\\OneDrive\\Documents\\lab\\280719\\Clusters Of Samples.csv")
# write.csv(cell_types_sum, file = "C:\\Users\\razsh\\OneDrive\\Documents\\lab\\280719\\Cell Types Sum.csv")

#saveRDS(lung, "C:\\Users\\razsh\\OneDrive\\Documents\\lab\\280719\\lung_seurat")
saveRDS(lung, "lung_allcells.rds")

lung_CancerOnly <- subset(lung, subset = cell.type == "cancer")
saveRDS(lung_CancerOnly, "lung_cancercells.rds")






                         
# # adding Endothelial cluster
# select.cells <- CellSelector(plot = umap_labled) # choosing cell from plot
# Idents(lung, cells = select.cells) <- "Endothelial"
# 
# lung@active.ident == "Endothelial"
# lung@meta.data$cell.type[lung@active.ident == "Endothelial"] <- "Endothelial"
# levels(lung@meta.data$cell.type) <-  c("cancer", "WBC", "mesothelial", "Endothelial")
# lung <- RenameIdents(lung, "Endothelial" = "16- Endothelial")
# Idents(lung)


table(lung@active.ident[lung@meta.data$patient.ident == "X760"])
table(lung@active.ident[lung@meta.data$patient.ident == "X793"])
table(lung@active.ident[lung@meta.data$patient.ident == "X802"])

