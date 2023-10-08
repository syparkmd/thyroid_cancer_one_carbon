#Generate_RNA_velocity_plot

if(T){
  library(tidyverse)
  library(ggsci)
  library(scales)
  library(ggsignif)
  library(corrplot)
  library(cowplot)
  library(ggrepel)
  library(velocyto.R) 
  library(mgsub)
}


#output path
out_dir = '/path/to/out_dir'
system(paste0('mkdir -p ',out_dir))

#design
theme_syp = theme_bw() + theme(axis.text = element_text(size=15), axis.title = element_text(size=18), panel.border = element_blank(), axis.line = element_line(), plot.title = element_text(size=20, face="bold"))
cancer_type_pal<- pal_d3("category10")(10)[c(1,3,2)]
names(cancer_type_pal) <- c('NT','PTC','ATC')

#load Seurat object -----

f4.tc.combined <- readRDS("/path/to/Seurat_object.rds")

f4.tc.combined@meta.data$cell_type2 %>% unique()
Idents(f4.tc.combined) <- "cell_type2"
thy.combined <- subset(f4.tc.combined, idents = c('thyroid_cell','TPO_high_thyroid_cell'))
Idents(thy.combined) <- "cancer_type"
f.thy.combined <- subset(thy.combined, idents = c('ptc','atc'))
target_cell_ids <- rownames(f.thy.combined@meta.data)
length(target_cell_ids)

if(file.exists(paste0(out_dir,'/velocyto_splice_data.rds'))==F){
  #load merged loom
  ldat <- read.loom.matrices("/path/to/thyroid_12s_merged.loom")
  #loading takes times
  orig.colnames<- ldat$spliced %>% colnames()
  as.character(lapply(orig.colnames, function(x) unlist(strsplit(x,':'))[1])) %>% unique()
  all(colnames(ldat$spliced)==colnames(ldat$unspliced)) #TRUE
  all(colnames(ldat$spliced)==colnames(ldat$ambiguous)) #TRUE
  new_colnames <- gsub("x$","-1",gsub("_GEX:","_",colnames(ldat$spliced)))
  new_colnames <- mgsub(new_colnames, c("Thy09","Thy13","Thy16","Thy17","Thy20"),c("AT9","AT13","AT16","AT17","AT20"))
  as.character(lapply(new_colnames, function(x) unlist(strsplit(x,'_'))[1])) %>% unique()
  
  m.ldat <- ldat
  colnames(m.ldat$spliced) <- new_colnames
  colnames(m.ldat$unspliced) <- new_colnames
  colnames(m.ldat$ambiguous) <- new_colnames
  target_cell_ids <- rownames(f.thy.combined@meta.data)
  setdiff(target_cell_ids, colnames(m.ldat$spliced)) #should be character(0)
  
  f.m.ldat <- m.ldat
  f.m.ldat$spliced <- m.ldat$spliced[,target_cell_ids]
  f.m.ldat$unspliced <- m.ldat$unspliced[,target_cell_ids]
  f.m.ldat$ambiguous <- m.ldat$ambiguous[,target_cell_ids]
  f.m.ldat %>% saveRDS(paste0(out_dir,'/velocyto_splice_data.rds'))  
}else{
  f.m.ldat <- readRDS(paste0(out_dir,'/velocyto_splice_data.rds'))  
  }

#extract matrix
emat <- f.m.ldat$spliced %>% as.matrix()
nmat <- f.m.ldat$unspliced %>% as.matrix()


#take embedding from PCA after Combat
if(T){
  pc <- readRDS(paste0('/path/to/thyroid_combat_pc.rds'))
  f_targets <- intersect(target_cell_ids, rownames(pc$x))
  pc_dt <- pc$x[f_targets,] %>% as.data.frame() 
  emb <- pc_dt[,c("PC1","PC2")] %>% as.matrix() #this should be matrix
  
  cell.dist <- as.dist(1-armaCor(t(pc_dt)))
  out_pdf=paste0(out_dir,'/velocity_pc_r100.pdf')
}
#take cluster labels from PCA and filter matrix
if(T){
  km_res <- kmeans(pc_dt[,1:2], 8)
  cluster.label <- km_res$cluster
  emb_dt <- emb %>% as.data.frame() %>% rownames_to_column("cell_id") %>% as_tibble()
  emb_dt <- left_join(emb_dt, tibble(cell_id = names(cluster.label), cluster=cluster.label))
  ggplot(emb_dt, aes(x=PC1, y=PC2))+
    geom_point(aes(color=as.factor(cluster)), alpha=0.5)
  
  group_dt <- read_tsv("/path/to/tds_sm_group.tsv")
  group_dt$tds_sm_group %>% unique()
  group_dt
  emb_dt <- left_join(emb_dt,group_dt)
  ggplot(emb_dt, aes(x=PC1, y=PC2))+
    geom_point(aes(color=tds_sm_group, alpha=0.5))
  emb_dt <- emb_dt %>% mutate(color=case_when(tds_sm_group == "T_hi_SM_lo" ~ "#1F77B4FF",
                                    tds_sm_group == "T_lo_SM_hi" ~  "#FF7F0EFF",
                                    TRUE ~ "grey"))
  cell.colors <- emb_dt$color
  names(cell.colors) <- emb_dt$cell_id
  emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
  nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
  length(intersect(rownames(emat),rownames(nmat))) #1502
}

fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,
                                            deltaT=1,kCells=20,
                                            cell.dist=cell.dist,
                                            fit.quantile=fit.quantile)

#plotting velocity
pdf(out_pdf, width=8, height=8)
show.velocity.on.embedding.cor(emb,rvel.cd,n=300,scale='sqrt',cex=1.5,
                               arrow.scale=10,show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5,grid.n=20,
                               arrow.lwd=1,do.par=F,cell.border.alpha = 0.1,
                               cell.colors="grey")
                               
dev.off()

