annotate_hsapien_rowname <- function(data){
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  attr.string = c('ensembl_gene_id', 'hgnc_symbol')
  gene.annotation = getBM(attributes=attr.string, 
                          filters =  'ensembl_gene_id', 
                          values = rownames(data), 
                          mart = ensembl)
  gene.annotation[duplicated(gene.annotation$hgnc_symbol)==TRUE,]$hgnc_symbol <- gene.annotation[duplicated(gene.annotation$hgnc_symbol)==TRUE,]$ensembl_gene_id
  names(gene.annotation) <- c("hgnc_symbol","ensembl_id")
  gene.annotation$hgnc_symbol <- sub('[.]', '_', make.names(gene.annotation$hgnc_symbol, unique=TRUE))
  print(gene.annotation)
}

annotate_hsapien_col <- function(data,col){
  library(biomaRt)
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  attr.string = c('ensembl_gene_id', 'hgnc_symbol')
  gene.annotation = getBM(attributes=attr.string, 
                          filters =  'ensembl_gene_id', 
                          values = data[,col], 
                          mart = ensembl)
  gene.annotation[duplicated(gene.annotation$hgnc_symbol)==TRUE,]$hgnc_symbol <- gene.annotation[duplicated(gene.annotation$hgnc_symbol)==TRUE,]$ensembl_gene_id
  names(gene.annotation) <- c("ensembl_id","hgnc_symbol")
  gene.annotation$hgnc_symbol <- sub('[.]', '_', make.names(gene.annotation$hgnc_symbol, unique=TRUE))
  print(gene.annotation)
}

#10x#

output_dir = "/Users/Gautam Kumar/Desktop/Output/"

seurat10x <- function(data_10x_dir, min_cells, min_features){
  data_dir = data_10x_dir
  data_10x = Read10X(data.dir = data_dir)
  seurat_object = CreateSeuratObject(counts = data_10x, min.cells = as.integer(min_cells), min.features = as.integer(min_features))
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
  seurat_object <- subset(seurat_object)
  seurat_object <- subset(seurat_object, features=grep("^mt-",rownames(seurat_object),invert = T,val=T))
  seurat_object <- NormalizeData(seurat_object)
  all.genes <- rownames(seurat_object)
  seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst")
  seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
}

specific_dim_res <- function(seurat_object, dim_num, res_num){
  dim_num <- as.numeric(dim_num)
  res_num <- as.numeric(res_num)
  
  seurat_object <- FindNeighbors(seurat_object, dims = 1:dim_num)
  seurat_object <- FindClusters(seurat_object, resolution = res_num)
  seurat_object <- RunUMAP(seurat_object, dims = 1:dim_num)
  
  DimPlot(object = seurat_object, reduction = 'umap', label = TRUE) + labs(title = paste0("dim_",as.character(dim_num),"_res_",as.character(res_num)))
  ggsave(paste0(output_dir,"UMAP_","dim_",as.character(dim_num),"_res_",as.character(res_num),".jpg"))
  
  all.markers <- FindAllMarkers(seurat_object)
  all.markers_filter <- all.markers[all.markers$p_val_adj < 0.05,]
  
  write.csv(all.markers_filter, paste0(output_dir,"UMAP_","dim_",as.character(dim_num),"_res_",as.character(res_num),".csv"))
}

specific_dim_res_plot <- function(seurat_object, dim_num, res_num){
  dim_num <- as.numeric(dim_num)
  res_num <- as.numeric(res_num)
  
  seurat_object <- FindNeighbors(seurat_object, dims = 1:dim_num)
  seurat_object <- FindClusters(seurat_object, resolution = res_num)
  seurat_object <- RunUMAP(seurat_object, dims = 1:dim_num)
  
  DimPlot(object = seurat_object, reduction = 'umap', label = TRUE) + labs(title = paste0("dim_",as.character(dim_num),"_res_",as.character(res_num)))
}

specific_dim_res_feature <- function(seurat_object, dim_num, res_num, feature){
  dim_num <- as.numeric(dim_num)
  res_num <- as.numeric(res_num)
  
  seurat_object <- FindNeighbors(seurat_object, dims = 1:dim_num)
  seurat_object <- FindClusters(seurat_object, resolution = res_num)
  seurat_object <- RunUMAP(seurat_object, dims = 1:dim_num)
  
  FeaturePlot(seurat_object, features = deparse(substitute(feature)))
}

specific_dim_res_feature2 <- function(seurat_object, dim_num, res_num, feature, pt_size){
  dim_num <- as.numeric(dim_num)
  res_num <- as.numeric(res_num)
  
  seurat_object <- FindNeighbors(seurat_object, dims = 1:dim_num)
  seurat_object <- FindClusters(seurat_object, resolution = res_num)
  seurat_object <- RunUMAP(seurat_object, dims = 1:dim_num)
  
  FeaturePlot(seurat_object, pt.size = pt_size, features = deparse(substitute(feature)))
}

bulk_dim_res <- function(seurat_object, dim_min, dim_max, dim_seq, res_min, res_max, res_seq){
  
  dim_range <- seq(as.numeric(dim_min), as.numeric(dim_max), by = as.numeric(dim_seq))
  res_range <- seq(as.numeric(res_min), as.numeric(res_max), by = as.numeric(res_seq))
  dim_res_combo<-expand.grid(dim_range,res_range)
  
  for(row in 1:nrow(dim_res_combo)) {
    dim_num <- dim_res_combo[row,"Var1"]
    res_num <- dim_res_combo[row,"Var2"]
    
    seurat_object <- FindNeighbors(object = seurat_object, dims = 1:as.numeric(dim_num))
    
    dir.create(paste0(output_dir,as.character(dim_num), " Dimensions"))
    dimension_directory <- paste0(output_dir, as.character(dim_num), " Dimensions","/")
    
    seurat_object <- FindClusters(object = seurat_object, resolution = as.numeric(res_num))
    seurat_object <- RunUMAP(object = seurat_object, dims = 1:as.numeric(dim_num))
    
    DimPlot(object = seurat_object, reduction = 'umap', label = TRUE) + labs(title = paste0("dim_",as.character(dim_num),"_res_",as.character(res_num)))
    ggsave(paste0(dimension_directory,"UMAP_","dim_",as.character(dim_num),"_res_",as.character(res_num),".jpg"))
    
    all.markers <- FindAllMarkers(seurat_object)
    write.csv(all.markers, file = paste0(dimension_directory,"UMAP_","dim_",as.character(dim_num),"_res_",as.character(res_num),".csv"))
  }
}

##############################################################

bulk_edger <- function(count_matrix,meta_data,interest_var){
  library(edgeR)
  dge_object <<- DGEList(counts = count_matrix, genes = rownames(count_matrix))
  
  countsPerMillion <- cpm(dge_object)
  countCheck <- countsPerMillion > 1
  keep <- which(rowSums(countCheck) >= 2)
  
  dge_object <- dge_object[keep,]
  dge_object <- calcNormFactors(dge_object)

  meta_data[,c(deparse(substitute(interest_var)))] <- as.factor(meta_data[,c(deparse(substitute(interest_var)))])
  dge_object_mat <<- model.matrix(~0+meta_data[,c(deparse(substitute(interest_var)))])
  colnames(dge_object_mat) <<- gsub("\\[.+?]", "", colnames(dge_object_mat))
  colnames(dge_object_mat) <<- gsub("meta_data", "", colnames(dge_object_mat))
  
  dge_object <- estimateGLMCommonDisp(dge_object, design=dge_object_mat)
  dge_object <- estimateGLMTrendedDisp(dge_object, design=dge_object_mat)
  dge_object <- estimateGLMTagwiseDisp(dge_object, design=dge_object_mat)
  
  fit_dge_object <<- glmFit(dge_object, dge_object_mat)
}

bulk_diff <- function(fit_dge_object, dge_object_mat, group1, group2){
  x <- paste0(substitute(group1),"-",substitute(group2))
  y <- makeContrasts(contrasts = x, levels = dge_object_mat)
  assign(paste0(deparse(substitute(group1)),"_v_",deparse(substitute(group2))),glmLRT(fit_dge_object, contrast = y), envir = parent.frame())
}