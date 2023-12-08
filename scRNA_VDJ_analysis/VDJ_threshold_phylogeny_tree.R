library(alakazam)
library(shazam)
library(tigger)
library(igraph)
library(dplyr)
genotype_dir <- ('Analysis/genotyping')
threshold_dir <- ('Analysis/threshold_estimation')
invisible( if (!dir.exists(genotype_dir)) {dir.create(genotype_dir, recursive=T)} )
invisible( if (!dir.exists(threshold_dir)) {dir.create(threshold_dir, recursive=T)} )
#######################################################
#### FIND NOVEL VDJ SEQ ALLELES AND INFER GENOTYPE ####
#######################################################

# load V-segment germline sequences
ighv <- readIgFasta( 'data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGHV.fasta')

# import ChangeO-formatted sequence database files (heavy chain)
seqdb <- NULL
for (i in c("NP-gB","NP-gH")){
  tmp <- readChangeoDb(paste0( i,"_vdj_b/",i,'_heavy_parse-select.tab'))
  tmp$SEQUENCE_ID <- paste0(tmp$SEQUENCE_ID, '-', i)
  tmp$MOUSE_NR <- i
  seqdb <- rbind(seqdb, tmp)
}

# infer genotype (performed separately for each mouse)
mice <- unique(seqdb$MOUSE_NR)
#mice <- mice[order(as.numeric(gsub('M','',mice)))]
default_threshold=0.1

predicted_thresholds <- data.frame(mouse=mice, threshold=rep(as.numeric(default_threshold), length(mice)))
for (m in mice) {
  
  cat('\n\nPROCESSING MOUSE', m, '\n\n')
  seqdb_mouse <- seqdb[seqdb$MOUSE_NR %in% m, ]
  colnames(seqdb_mouse)=tolower(colnames(seqdb_mouse))
  # # find novel alleles (if any)
  # cat('\tSearching for novel alleles ...\n')
  # novel_rows <- NULL
  # try (nv <- findNovelAlleles(seqdb_mouse, ighv,seq = "sequence_imgt"))
  # try (novel_rows <- selectNovel(nv))
  # 
  # # Extract and view the rows that contain successful novel allele calls
  # if (!is.null(novel_rows) && (nrow(novel_rows) > 0)) {
  #   png(filename=paste0(genotype_dir, '/novel_alleles_', m, '.png'), units='mm', height=250, width=180, res=300)
  #   plotNovel(seqdb_mouse, novel_rows[1, ])  # only plot first novel allele
  #   invisible(dev.off())
  # }
  
  # infer mouse genotype
  gt_mouse <- inferGenotype(seqdb_mouse, germline_db=ighv)
  png(filename=paste0(genotype_dir, '/IGHV_genotype_plot_', m, '.png'), units='mm', height=250, width=100, res=300)
  plotGenotype(gt_mouse, gene_sort="position", text_size=8) #, facet_by='ALLELES')
  invisible(dev.off())
  
  # convert genotype table to vector of nucleotide sequences
  gtseq_mouse <- genotypeFasta(gt_mouse, germline_db=ighv)
  writeFasta(gtseq_mouse, paste0(genotype_dir, '/IGHV_genotype_', m, '.fasta'))
  
  # correct allele calls based on the personalized genotype
  seqdb_mouse <- reassignAlleles(seqdb_mouse, gtseq_mouse)
  #---------
  
  
  ###################################
  ### CALCULATE NEAREST NEIGHBORS ###
  ###################################
  
  # calculate distances to nearest neighbors
  dist_ham <- distToNearest(seqdb_mouse, vCallColumn="v_call_genotyped", model="ham", normalize="len", nproc=1)
  
  # plot distance distribution
  png(filename=paste0(threshold_dir, '/distToNearestNeighbor_', m, '.png'), units='mm', height=150, width=180, res=300)
  print(ggplot(subset(dist_ham, !is.na(dist_nearest)), aes(x=dist_nearest)) +
          theme_bw() +
          xlab("Hamming distance") +
          ylab("Count") +
          scale_x_continuous(breaks=seq(0, 1, 0.1)) +
          geom_histogram(color="white", binwidth=0.02))
  invisible(dev.off())
  
  # perform automatic threshold estimation
  dist_args <- trimws(unlist(strsplit(casefold('gmm'), ',')))
  if (dist_args[1] == 'none') {
    
    cat('\n\tSkipping automatic threshold estimation. Threshold value set to default:', default_threshold, '\n')
    
  } else if (dist_args[1] %in% c('density','gmm')) {
    
    # Find threshold using the density or gmm (mixture model) methods
    if (dist_args[1] == 'density') {
      output <- findThreshold(dist_ham$dist_nearest, method=dist_args[1])
    } else {
      if (is.na(dist_args[2])) { dist_args[2] <- 'gamma-gamma' }
      output <- findThreshold(dist_ham$dist_nearest, method=dist_args[1], model=dist_args[2])  
    }
    
    if (is.null(output) || is.na(output@threshold)) {  
      cat('\n\tThreshold estimation failed. Reverting to default threshold:', default_threshold, '\n')
    } else {
      cat('\n\tEstimated hamming distance threshold:', round(output@threshold, 5), '\n')
      predicted_thresholds$threshold[predicted_thresholds$mouse == m] <- output@threshold
      png(paste0(threshold_dir, '/distToNearestNeighbor_', paste(dist_args, collapse='_'), '_fit_', m, '.png'), units='mm', height=120, width=150, res=300)
      plot(output, binwidth=0.02, title=paste0('Threshold Prediction [threshold = ', round(output@threshold, 3), ']'))
      invisible(dev.off())
    }
    
  } else {
    stop(paste0('Invalid density_method: "', dist_args[1], '". Valid options are "density", "gmm", or "none".'))
  }
  #---------
  
  
  ###########################
  ### EXPORT DATA TO FILE ###
  ###########################
  writeChangeoDb(seqdb_mouse, paste0(genotype_dir, '/IGHV-genotyped_', m, '.tab'))
  write.csv(predicted_thresholds, file=paste0(threshold_dir, '/predicted_thresholds.csv'), quote=F, row.names=F)
  
}
#---------

#### Lineage phylogeny tree ####

# Need to download and install PHYLIP from http://evolution.genetics.washington.edu/phylip/getme-new1.html
# path to dnapars executable within PHYLIP package
phylip_exec <- "/home/rs1/1-software/PHYLIP/phylip-3.697/exe/dnapars"
# specify directory containing genotyped VDJ files
proj_dir <- '/home/project7/hyzhang_2023/e007_mm_10x_multi/'
geno_dir <- paste0(proj_dir, 'Analysis/genotyping/')

# get list of available files
geno_files <- dir(geno_dir, 'germ-pass[.]tab', full.names=T)
ann1=readRDS("annotation_meta.RDS")
ann1=as.data.frame(ann1)
ann1$barcode=rownames(ann1)
####### NP-gB #########
g_file <- geno_files[1]
# load the Change-O database file with germline sequence information (*_germ-pass.tab file)
db <- readChangeoDb(g_file)
colnames(db)=tolower(colnames(db))
db$barcode=unlist(lapply(db$sequence_id,function(x){strsplit(x,'[_]')[[1]][1] }))
db$barcode=paste0(db$barcode,"_1")
db2=merge(db,ann1,by = "barcode")
db2$ann1=as.character(db2$ann1)
# select a desired clone
head(sort(table(db2$clone), decreasing=T))  # #252  62 115 176 183  67
sub_db <- subset(db2, clone == 252)
write.csv(sub_db,"clone_67_NP-gB.csv")
# create ChangeOclone object for clone
clone <- makeChangeoClone(sub_db, text_fields=c('c_call','ann1'), clone = "clone",
                          num_field='umicount',seq = "sequence_imgt",germ = "germline_imgt")
# Run PHYLIP and parse output
graph1 <- buildPhylipLineage(clone, phylip_exec, rm_temp=T, verbose=T)
# Modify graph and plot attributes
V(graph1)$color[V(graph1)$name == "Germline"] <- "black"
V(graph1)$color[grepl("Inferred", V(graph1)$name)] <- "white"
V(graph1)$color[grepl("Activated_B", V(graph1)$ann1)] <- "#FFADAD"
V(graph1)$color[grepl("Naive_B", V(graph1)$ann1)] <- "#0F7B6D"
V(graph1)$color[grepl("GC-LZ_B", V(graph1)$ann1)] <- "#6940A6"
V(graph1)$color[grepl("GC-DZ_B", V(graph1)$ann1)] <- "#DFAB00"
V(graph1)$color[grepl("PB", V(graph1)$ann1)] <- "#A0C4FF"
V(graph1)$shape <- "circle"
V(graph1)$shape[grepl("IGHG2C", V(graph1)$c_call)] <- "square"
V(graph1)$shape[grepl("IGHA", V(graph1)$c_call)] <- "sphere"
V(graph1)$shape[grepl("IGHG2B", V(graph1)$c_call)] <- "rectangle"
V(graph1)$label <- V(graph1)$c_call
E(graph1)$label <- ""

# Remove large default margins
par(mar=c(0, 0, 0, 0) + 0.1)
# Plot graph
plot(graph1, layout=layout_as_tree, edge.arrow.mode=0, vertex.frame.color="black",
     vertex.label.color="black", vertex.size=5, vertex.label.family='sans', vertex.label.cex=0.65,
     vertex.label.cex=0.7,vertex.label.dist=1.2)
# Add legend
legend("topleft", c("Germline", "Inferred", "Activated_B","Naive_B","GC-LZ_B","GC-DZ_B","PB"), 
       fill=c("black", "white", "#FFADAD","#0F7B6D", "#6940A6","#DFAB00","#A0C4FF"), cex=0.75)



####### NP-gH #########
g_file <- geno_files[2]
# load the Change-O database file with germline sequence information (*_germ-pass.tab file)
db <- readChangeoDb(g_file)
colnames(db)=tolower(colnames(db))
db$barcode=unlist(lapply(db$sequence_id,function(x){strsplit(x,'[_]')[[1]][1] }))
db$barcode=paste0(db$barcode,"_2")
db2=merge(db,ann1,by = "barcode")
db2$ann1=as.character(db2$ann1)
# select a desired clone
head(sort(table(db2$clone), decreasing=T))  # 15 135 204 232 272 302 
sub_db <- subset(db2, clone == 302)
write.csv(sub_db,"clone_302_NP-gH.csv")
# create ChangeOclone object for clone
clone <- makeChangeoClone(sub_db, text_fields=c('c_call','ann1'), clone = "clone",
                          num_field='umicount',seq = "sequence_imgt",germ = "germline_imgt")
# Run PHYLIP and parse output
graph1 <- buildPhylipLineage(clone, phylip_exec, rm_temp=T, verbose=T)
# Modify graph and plot attributes
V(graph1)$color[V(graph1)$name == "Germline"] <- "black"
V(graph1)$color[grepl("Inferred", V(graph1)$name)] <- "white"
V(graph1)$color[grepl("Activated_B", V(graph1)$ann1)] <- "#FFADAD"
V(graph1)$color[grepl("Naive_B", V(graph1)$ann1)] <- "#0F7B6D"
V(graph1)$color[grepl("GC-LZ_B", V(graph1)$ann1)] <- "#6940A6"
V(graph1)$color[grepl("GC-DZ_B", V(graph1)$ann1)] <- "#DFAB00"
V(graph1)$color[grepl("PB", V(graph1)$ann1)] <- "#A0C4FF"
V(graph1)$shape <- "circle"
V(graph1)$shape[grepl("IGHG2C", V(graph1)$c_call)] <- "square"
V(graph1)$shape[grepl("IGHA", V(graph1)$c_call)] <- "sphere"
V(graph1)$shape[grepl("IGHG2B", V(graph1)$c_call)] <- "rectangle"
V(graph1)$label <- V(graph1)$c_call
E(graph1)$label <- ""

# Remove large default margins
par(mar=c(0, 0, 0, 0) + 0.1)
# Plot graph
plot(graph1, layout=layout_as_tree, edge.arrow.mode=0, vertex.frame.color="black",
     vertex.label.color="black", vertex.size=5, vertex.label.family='sans', vertex.label.cex=0.5,
     vertex.label.cex=0.7,vertex.label.dist=1)
# Add legend
legend("topleft", c("Germline", "Inferred", "Activated_B","Naive_B","GC-LZ_B","GC-DZ_B","PB"), 
       fill=c("black", "white", "#FFADAD","#0F7B6D", "#6940A6","#DFAB00","#A0C4FF"), cex=0.75)
