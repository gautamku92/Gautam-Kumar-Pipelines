library(Seurat)  

input_10x_specific<-function(dir_10x){
  suppressPackageStartupMessages(library(Seurat))
  data <- Read10X(data.dir = dir_10x)
  seurat_object <<- CreateSeuratObject(counts = data, min.cells=0, min.features=0)
  seurat_object <<- subset(seurat_object, features=grep("^mt-",rownames(seurat_object),invert=T,val=T))
  seurat_object <<- NormalizeData(seurat_object)
  seurat_object <<- FindVariableFeatures(seurat_object, selection.method = "vst")
  seurat_object <<- ScaleData(seurat_object, features = rownames(seurat_object))
  seurat_object <<- RunPCA(seurat_object)
  dim_sel <<- as.numeric(readline(prompt = "Number of dimensions: "))
  res_sel <<- as.numeric(readline(prompt = "Resolution: "))
  seurat_object <<- FindNeighbors(object = seurat_object, dims = 1:dim_sel)
  seurat_object <<- FindClusters(object = seurat_object, resolution = res_sel)
  seurat_object <<- RunUMAP(object = seurat_object, dims = 1:dim_sel)
  DimPlot(object = seurat_object, reduction = "umap")
}

input_10x_range<-function(dir_10x){
  suppressPackageStartupMessages(library(Seurat))
  data <- Read10X(data.dir = dir_10x)
  seurat_object <<- CreateSeuratObject(counts = data, min.cells=0, min.features=0)
  seurat_object <<- subset(seurat_object, features=grep("^mt-",rownames(seurat_object),invert=T,val=T))
  seurat_object <<- NormalizeData(seurat_object)
  seurat_object <<- FindVariableFeatures(seurat_object, selection.method = "vst")
  seurat_object <<- ScaleData(seurat_object, features = rownames(seurat_object))
  seurat_object <<- RunPCA(seurat_object)

  output_dir <<- as.character(readline(prompt = "Select output directory: "))
  
  dim_low <<- as.numeric(readline(prompt = "Lowest number of dimensions: "))
  dim_hi <<- as.numeric(readline(prompt = "Highest number of dimensions: "))
  dim_seq <<- as.numeric(readline(prompt = "Sequence of dimensions to select: "))
  res_low <<- as.numeric(readline(prompt = "Lowest resolution: "))
  res_hi <<- as.numeric(readline(prompt = "Highest resolution: "))
  res_seq <<- as.numeric(readline(prompt = "Sequence of resolution to select: "))
  
  dim_range <- seq(as.numeric(dim_low), as.numeric(dim_hi), by = as.numeric(dim_seq))
  res_range <- seq(as.numeric(res_low), as.numeric(res_hi), by = as.numeric(res_seq))
  dim_res_combo<-expand.grid(dim_range,res_range)
  
  for(row in 1:nrow(dim_res_combo)) {
    dim_num <- dim_res_combo[row,"Var1"]
    res_num <- dim_res_combo[row,"Var2"]
    
    seurat_object <- FindNeighbors(object = seurat_object, dims = 1:as.numeric(dim_num))
    
    dir.create(paste0(output_dir,"/",as.character(dim_num), " Dimensions"))
    dimension_directory <- paste0(output_dir,"/", as.character(dim_num), " Dimensions","/")
    
    seurat_object <- FindClusters(object = seurat_object, resolution = as.numeric(res_num))
    seurat_object <- RunUMAP(object = seurat_object, dims = 1:as.numeric(dim_num))
    
    DimPlot(object = seurat_object, reduction = 'umap', label = TRUE) + labs(title = paste0("dim_",as.character(dim_num),"_res_",as.character(res_num)))
    ggsave(paste0(dimension_directory,"UMAP_","dim_",as.character(dim_num),"_res_",as.character(res_num),".jpg"))
    
    all.markers <- FindAllMarkers(seurat_object)
    write.csv(all.markers, file = paste0(dimension_directory,"UMAP_","dim_",as.character(dim_num),"_res_",as.character(res_num),".csv"))
  }
}

##################################################
#####################BULK#########################
##################################################

bulk_data_check<-function(count_matrix, meta_data){
  print("BEFORE RUNNING, MAKE SURE COUNT MATRIX HAS GENE NAMES AS ROWNAMES AND SAMPLE/CELL NAMES AS COLNAMES")
  meta_names <- readline(prompt="Input the column name in the metadata that contains the names of the subjects/samples: ")
  meta_names <- meta_data[,meta_names]
  print(colnames(count_matrix) %in% meta_names)
  print("IF THE ABOVE COMES OUT AS ENTIRELY/MOSTLY FALSE, CHECK NAMING SYSTEM IN BOTH DATASETS")
  meta_info <- readline(prompt="Input the column name in the metadata that contains the variable you are interested in observing: ")
  meta_info <- meta_data[,meta_info]
  print(is.na(meta_info))
  print(grepl(" ",as.character(meta_info)))
  print("IF THERE ARE ANY TRUES ABOVE, CHECK FOR EITHER MISSING DATA OR EMPTY CELLS")
  print(length(meta_info))
  print(length(colnames(pvc)))
  print("IF THE NUMBERS DO NOT MATCH ABOVE, YOU NEED TO SUBSET THE COUNT MATRIX TO INCLUDE ONLY THE CELLS THAT HAVE THE METADATA OF INTEREST")
}

bulk_data_prep<-function(count_matrix, meta_data){
  print("BEFORE RUNNING, MAKE SURE THE NAMING SYSTEM AND MISSING DATA ARE TAKEN CARE OF")
  meta_subset <<- meta_data[,readline(prompt="Input the column name in the metadata that contains the variable of interest: ")]
  meta_clean<<-meta_data[grepl(".",meta_subset) == TRUE,]
  meta_clean<<-meta_clean[grepl(" ",meta_subset) == FALSE,]
  meta_clean<<-meta_clean[is.na(meta_subset) == FALSE,]
  selected<<-meta_clean[,readline(prompt="Input the column name in the metadata that contains the names of the subjects/samples: ")]
  count_clean<<-count_matrix[selected %in% colnames(count_matrix),]
  meta_clean<<-meta_clean[selected %in% colnames(count_matrix),]
  print(dim(meta_clean))
  print(dim(count_clean))
  print("THE NUMBER OF ROWS IN THE METADATA SHOULD MATCH THE NUMBER OF COLUMNS IN THE COUNT DATA")
  print(rownames(meta_clean) %in% colnames(count_clean))
  print("IF THE ABOVE COMES OUT AS ENTIRELY FALSE, MAKE SURE THE SAMPLE NAMES ARE IN THE ROWS OF THE METADATA AND IN THE COLUMNS OF THE COUNT DATA")
}

bulk_edger <- function(count_matrix,meta_data,interest_var){
  print("MAKE SURE THAT THE VARIABLE OF INTEREST IS IN CHARACTER FORMAT")
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