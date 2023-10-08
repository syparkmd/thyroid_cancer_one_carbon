#make cell type signature matrix for CIBERSORTx

if(T){
  library(tidyverse)
  library(Seurat)
  library(ggsci)
  library(scales)
  library(ggsignif)
  library(corrplot)
  library(cowplot)
  library(ggrepel)
}


#output path
out_dir = '/path/to/out_dir'
system(paste0('mkdir -p ',out_dir))

#design
theme_syp = theme_bw() + theme(axis.text = element_text(size=15), 
                               axis.title = element_text(size=18), 
                               panel.border = element_blank(), 
                               axis.line = element_line(), 
                               plot.title = element_text(size=20, face="bold"))
cancer_type_pal<- pal_d3("category10")(10)[c(1,3,2)]
names(cancer_type_pal) <- c('Normal','PTC','ATC')

#source functions
source("/path/to/sypark_scRNAseq_plot_utils.R")

#load Seurat object
if(T){
  f4.tc.combined <- readRDS("/path/to/tc_combined_pc50_f4.rds")
}


#check cell type
f4.tc.combined@meta.data %>% colnames()
draw_umap_plot(f4.tc.combined, "cell_type2")

#cell type3 assign for publication figure
tmp_meta <- f4.tc.combined@meta.data %>% as.data.frame() %>% rownames_to_column('cell_id')
tmp_meta <- tmp_meta %>%
  mutate(cell_type3 = case_when(cell_type2 == "TPO_high_thyroid_cell" ~ "Thyroid cell",
                               cell_type2 == "thyroid_cell" ~ "Thyroid cell",
                               cell_type2 == "fibroblast" ~ "Fibroblast",
                               cell_type2 == "myofibroblast" ~ "Fibroblast",
                               cell_type2 == "macrophage" ~ "Macrophage",
                               cell_type2 == "CCL17_macrophage" ~ "Macrophage",
                               cell_type2 == "endothelial_cell" ~ "Endothelial cell",
                               cell_type2 == "dendritic_cell" ~"Dendritic cell",
                               cell_type2 == "plasma_cell" ~ "Plasma cell",
                               cell_type2 == "B_cell" ~ "B cell",
                               cell_type2 == "mast_cell" ~ "Mast cell",
                               cell_type2 == "Treg_cell" ~ "Regulatory T cell",
                               cell_type2 == "CD4_T_cell" ~ "CD4 T cell",
                               cell_type2 == "CD8_T_cell" ~ "CD8 T cell",
                               cell_type2 == "ISG_T_cell" ~ "IFN-stimulated T cell",
                               cell_type2 == "KIT_high_T_cell" ~ "CD8 T cell",
                               cell_type2 == "NK_cell" ~ "NK cell"))

f4.tc.combined@meta.data <- tmp_meta %>% as.data.frame() %>% column_to_rownames('cell_id')
draw_umap_plot(f4.tc.combined, 'cell_type3')


#randomly select 30 cells per cell type to make CIBERSORTx input matrix
tmp <- f4.tc.combined@meta.data %>% as.data.frame() %>% rownames_to_column("cell_id") %>%
  as_tibble() %>% select(cell_id, cell_type3)
set.seed(1);r_tmp <-  tmp %>% group_by(cell_type3) %>% sample_n(30) 
#make unique cell type name by appending number to cell type
r_tmp <- r_tmp %>% mutate(cell_type_unique=make.unique(cell_type3))

# re-normalize the seurat object to get CPM value
f4.tc.combined.rc <- NormalizeData(f4.tc.combined, normalization.method = "RC", 
                                   scale.factor=1e6)
#get CPM matrix per cell per selected gene
if(T){
  cpm_mx <- f4.tc.combined.rc@assays$RNA@data[selected_genes,r_tmp$cell_id] %>% as.matrix()
  dim(cpm_mx) #341, 390
  ComplexHeatmap::Heatmap(cpm_mx, show_column_names = F, show_row_names = F)
}

#get CPM matrix per cell in whole genes
cpm_mx <- f4.tc.combined.rc@assays$RNA@data[,r_tmp$cell_id] %>% as.matrix()
dim(cpm_mx) #28745, 390

#replace cell name to cell type name
tmp1 <- cpm_mx %>% t() %>% as.data.frame() %>% rownames_to_column("cell_id") %>%
  as_tibble()
tmp1 <- left_join(tmp1, r_tmp %>% select(cell_id, cell_type_unique))
tmp1 <- tmp1 %>% select(-cell_id, -cell_type3) %>% as.data.frame() %>% column_to_rownames("cell_type_unique") %>%
  t() %>% as.data.frame() %>% rownames_to_column("GeneSymbol") %>% as_tibble()
tmp1 <- tmp1 %>% mutate_at(vars(-GeneSymbol), ~ as.numeric(as.character(.)))
tmp1 %>% write_tsv(paste0(out_dir,"/single_cell_390cells.tsv"))

#make mixture file format using bulk RNAseq dataset
#load bulk RNAseq data
bulk_primary_ec_bc_path="/path/to/thy_primary_tumor_ec_bc_mx.rds"
pri_dt <- readRDS(bulk_primary_ec_bc_path)
#convert readcount to CPM data and save tsv
pri_cpm_dt <- apply(pri_dt,2, function(x) (x/sum(x))*1000000) 
pri_cpm_dt %>% as.data.frame() %>% rownames_to_column("Gene") %>% as_tibble() %>% 
  write_tsv(paste0(out_dir,"/bulk_369_cpm.tsv"))
