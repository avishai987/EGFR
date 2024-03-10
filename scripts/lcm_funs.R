
# remove irrelevant genes and columns --------------------------------------

remove_lines <- function(data,omit_mt) {
  rownames(data)=make.unique(data$gene_name)
  if (omit_mt == TRUE){
    omitgenes=startsWith(rownames(data),"MT-")|startsWith(rownames(data),"ERCC-") #remove MT and ERCC
  }
  
  if (omit_mt == FALSE){
    omitgenes=startsWith(rownames(data),"ERCC-")
  }
  
  data= data[!omitgenes,]
  data= data[,7:ncol(data)] #remove first 7 columns
  return(data)
}




# change row names --------------------------------------------------------
change_cell_names <- function(data) {
  cell.labels <- gsub(".sort.bam","",colnames(data))
  cell.labels <- gsub("\\.","_",cell.labels)
  cell.labels <- gsub("run.*MGH.*MGH","",cell.labels)
  cell.labels = paste0("MGH",cell.labels)
  cell.labels <- gsub("ontreatment","on-treatment",cell.labels)
  cell.labels <- gsub("pretreatment","pre-treatment",cell.labels)
  cell.labels <- gsub("pre_treatment","pre-treatment",cell.labels)
  cell.labels <- gsub("on_treatment","on-treatment",cell.labels)
  
  colnames(data)=cell.labels #apply
  return(data)
}



# create Seurat object ----------------------------------------------------

create_seurat <- function(data,normalize,logTransform) {
  if (logTransform == T){
    egfr <- CreateSeuratObject(counts = log2(data+1), project = "EGFR", min.cells = 3, min.features = 2000)
  }
  
  if (logTransform == F){
    egfr <- CreateSeuratObject(counts = data, project = "EGFR", min.cells = 3, min.features = 2000)
  }
  
  percent.mt <- PercentageFeatureSet(egfr, pattern = "^MT-") #add MT genes percentages
  egfr[["percent.mt"]] <- percent.mt
  
  if (normalize == T){
    egfr <- NormalizeData(egfr) 
  }
  
  
  egfr <- AddMetaData(object = egfr, metadata = colnames(egfr), col.name = "orig.ident")
  return (egfr)
}



# Scaling -----------------------------------------------------------------

scale_data <- function(egfr) {
  
  # Identification of highly variable features
  egfr <- FindVariableFeatures(egfr, selection.method = "vst", nfeatures = 15000)
  
  # Scaling the data
  egfr <- ScaleData(egfr, vars.to.regress = c("nCount_RNA","percent.mt"))
  return(egfr)
}




# Dimensionality reductions -----------------------------------------------
Dim_reduction <- function(egfr) {
  clus_res=1
  pc2use=1:10;
  
  # Perform linear dimensional reduction (PCA)
  egfr <- RunPCA(egfr, features = VariableFeatures(object = egfr))
  pcaplot <- PCAPlot(egfr, label = FALSE)
  elbowplot <- ElbowPlot(egfr, ndims = 50) # checking the dimensionality 
  
  # cluster the cells
  egfr <- FindNeighbors(egfr, dims = pc2use)
  egfr <- FindClusters(egfr, resolution = clus_res)
  
  # Run non-linear dimensional reduction (UMAP)
  egfr <- RunUMAP(egfr, dims = pc2use)
  print(DimPlot(object = egfr, reduction = "umap", pt.size = 3, label = TRUE))
  return(egfr)
  
}


# #Add scores to dataset --------------------------------------------------

add_score <- function(genes_vector,score_name,dataset) {
  var_features=dataset@assays$RNA@var.features
  
  new_score=apply(dataset@assays$RNA@scale.data[intersect(genes_vector,var_features),],2,mean)
  dataset=AddMetaData(dataset,new_score,score_name)
}


# Add patient identity and time point to metadata ----------------------------------------

patient.iden_toMetadata <- function(egfr) {
  patient.ident = colnames(egfr)
  patient.ident = substring(patient.ident, regexpr("MGH", patient.ident),regexpr("_(?:pre-|on-|progression)", patient.ident)-1)
  egfr <- AddMetaData(object = egfr, metadata = as.factor(patient.ident), col.name = "patient.ident")
  
  egfr = SetIdent(egfr, value = egfr@meta.data$patient.ident)
  
  # DimPlot(object = egfr, reduction = "umap", pt.size = 3, label = TRUE, group.by ="patient.ident")
  return(egfr)
  
}


time.point_toMetadata <- function(egfr) {
  time.point = colnames(egfr)
  time.point = substring(time.point, regexpr("_(?:pre-|on-|progression)", time.point)+1,regexpr("_EGF816", time.point)-1)
  egfr <- AddMetaData(object = egfr, metadata = as.factor(time.point), col.name = "time.point")
  
  # DimPlot(object = egfr, reduction = "umap", pt.size = 3, label = TRUE, group.by ="time.point" )
  return(egfr)
  
}

# Add cell identity (tumor/stroma) --------------------------------------------------------
add_tumor_stroma <- function(egfr) {
  #Get sample ID from file
  cell_type = read.csv(file = "./Data/LCM/tumor or stroma.csv",header = T)
  names(cell_type)[1] = "Sample_ID"
  cell_type$Sample_ID <- gsub("EGF816_","",cell_type$Sample_ID)
  cell_type$test = paste0(cell_type$Patient,"_", cell_type$Sample_ID)
  
  #make orig.ident as the tumor/stroma file
  orig.ident = colnames(egfr)
  orig.ident2 = substring(orig.ident, regexpr("_(?:pre-|on-|progression)", orig.ident),regexpr("_EGF816", orig.ident)+6)
  orig.ident3 <- sapply(1:length(orig.ident2), function(x) gsub(orig.ident2[x], "", orig.ident[x]))
  orig.ident3 = as.data.frame(orig.ident3)
  
  #add type to orig.ident
  orig.ident3$type <- sapply(1:nrow(orig.ident3), function(x) cell_type$Tumor.Non.tumor[match(orig.ident3[x,],cell_type$test)])
  
  #save to metadata
  egfr <- AddMetaData(object = egfr, metadata = as.factor(orig.ident3$type), col.name = "cell.origin")

  return (egfr)
}




# find outliered cells  -----------------------------------------------------
library(stringi)
find_outliers <- function(data, save_plots = F,prefix = "") {
  patient_vector = unique(data@meta.data[["patient.ident"]])
  patient_vector = patient_vector[!patient_vector %in% c("MGH10124"," MGH1184")] #exclude these patients because they don't have enough cells
  
  cells_annotation = data.frame() #all cells clusters
  if(save_plots == T) {pdf(file = prefix %s+% "outliers cells heatmap.pdf",width = 12.25, height = 8.89 )}
  
  for (patient in patient_vector) {
    patient_data = subset(x = data, subset = patient.ident == patient)
    plt = plot_corr_patient_tp(as.data.frame(patient_data@assays[["RNA"]]@data), main  = patient,
                               cluster = T,return_plot = T, silent = T) #create correlation plot of all cells
    
    myannotation = as.data.frame(cutree(plt$tree_row,k = 2)) #divide into 2 clusters
    names(myannotation)[1] = "cluster"
    cluster_1_count = sum(myannotation$cluster == 1)
    cluster_2_count = sum(myannotation$cluster == 2)
    
    #keep the bigger cluster of cells
    if( cluster_1_count > cluster_2_count){
      Keep = 1
      Filter = 2
    }
    
    if( cluster_1_count < cluster_2_count){
      Keep = 2
      Filter = 1
    }
    
    myannotation[ myannotation == Keep] = "Normal"
    myannotation[ myannotation == Filter] = "Outlier"
    cells_annotation = rbind(cells_annotation,myannotation)
    
    #plot again with cluster annotation
    p = plot_corr_patient_tp(as.data.frame(patient_data@assays[["RNA"]]@data), main  = patient,
                             cluster = T,return_plot = F,annotation = myannotation, show_colnames = T) 
    print(p)
  }
  if(save_plots == T) {dev.off()}
  
  filtered_out = row.names(cells_annotation)[cells_annotation$cluster == "Outlier"] #outliers list
  tryCatch({
    patient_data = subset(x = data, subset = patient.ident == "MGH10124") #add the cell from MGH10124
    cells_annotation = rbind(cells_annotation,data.frame(cluster = "Normal", row.names = colnames(patient_data)))
    
  }, error=function(e){})
  
  
  data = AddMetaData(object = data, metadata = cells_annotation)
  return(data)
}


