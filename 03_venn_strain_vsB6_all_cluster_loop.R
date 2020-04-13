# making Venn diagram to show the intersection of the genes from strain effect: WSB, PWK, CAST vs B6

library(tidyverse)
library(ggrepel)

output <- "./03_venn_strain_vsB6/"

# make overlaping gene with GWAS genes
# the pool of GWAS gene is from paper "GWAS on family history of Alzheimer’s disease" https://www.nature.com/articles/s41398-018-0150-6#Sec8 
# table S5 MAGMA gene-based associations for the UK Biobank and IGAP meta-analysis summary results
# downloaded to "03_venn_strain_vsB6/41398_2018_150_MOESM2_ESM.xlsx"

library(readxl)
GWAS <- read_excel(path = "03_venn_strain_vsB6/41398_2018_150_MOESM2_ESM.xlsx", sheet = "Table S5", skip = 2)
names(GWAS) <- gsub(pattern = " ", replacement = "", names(GWAS))
names(GWAS) <- gsub(pattern = "\\+", replacement = "_", names(GWAS))
GWAS <- GWAS %>% mutate(Symbol=str_to_title(SYMBOL)) %>% filter(PUKB_IGAP<0.05)
GWAS_gene <- GWAS  %>% select(Symbol) %>% unlist() 

clusters <- c("H", "6", "7", "8", "9", "10", "11", "12")
DE_files <- paste("strain_vsB6_", clusters, ".txt", sep = "")
DE_files <- file.path("./02_combine/DE_vsB6", DE_files)


for (i in seq_along(clusters)){
  cluster <- clusters[[i]]
  df <- read_delim(DE_files[[i]], delim = "\t")
  
  split_strain <- function(strain, data){
    df_subset <- data %>% select(ENSMUSG, Orig_Symbol, contains(strain)) 
    names(df_subset)[3:4] <- c("logFC", "FDR")
    df_subset <- df_subset %>% filter(FDR<0.05, abs(logFC)>=1)
    return(df_subset$Orig_Symbol %>% unlist())
  }
  
  DE_CAST <- split_strain("CAST", df)
  DE_PWK <- split_strain("PWK", df)
  DE_WSB <- split_strain("WSB", df)
  
  library(VennDiagram)
  
  # Prepare a palette of 3 colors with R colorbrewer:
  strain_colors <- c("#888888", "#00AA00", "#FF0000", "#9900EE") %>% as.list()
  names(strain_colors) <- c("B6", "CAST", "PWK", "WSB")
  strains <- c("CAST", "PWK", "WSB")
  myCol <- c(strain_colors[[strains[1]]], strain_colors[[strains[2]]], strain_colors[[strains[3]]])
  
  # Chart
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  venn.diagram(
    x = list(DE_CAST, DE_PWK, DE_WSB),
    category.names = paste(strains, "vs", "B6", sep=" "),
    filename = paste(output, "cluster_", cluster, "_venn_vsB6", ".png", sep=""),
    #output=TRUE,
    
    # Output features
    imagetype="png" ,
    height = 540 , 
    width = 540 , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = myCol,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "Arial",
    
    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "Arial",
    rotation = 1
  )
  
  
  
  
  GWAS_color <- "#3792cb"
  myCol <- c(strain_colors[[strains[1]]], strain_colors[[strains[2]]], strain_colors[[strains[3]]], GWAS_color)
  
  # Chart
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  venn.diagram(
    x = list(DE_CAST, DE_PWK, DE_WSB, GWAS_gene),
    category.names = c(paste(strains, "vs", "B6", sep=" "), "GWAS"),
    filename = paste(output, "cluster_", cluster, "_venn_vsB6_GWAS", ".png", sep=""),
    #output=TRUE,
    
    # Output features
    imagetype="png" ,
    height = 520 , 
    width = 620 , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = myCol,
    
    # Numbers
    cex = .4,
    fontface = "bold",
    fontfamily = "Arial",
    
    # Set names
    cat.cex = 0.4,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    #cat.pos = c(-27, 27, 135),
    #cat.dist = c(0.055, 0.055, 0.085),
    #cat.fontfamily = "Arial",
    #rotation = 1
  )
  
  ## get the genelist of each intersection category
  
  library(gplots)
  input = list(DE_CAST, DE_PWK, DE_WSB, GWAS_gene)
  x <- venn(input, names=c(paste(strains, "vs", "B6", sep=" "), "GWAS"))
  # the venn function in gplots is more intuitive in terms of what interaction to look at
  
  intersections <- attr(x, "intersections")
  intersections_length <- map(intersections, length)
  intersections_sum <- map_df(intersections, length) %>% t() %>% as.data.frame()
  intersections_sum <- rownames_to_column(intersections_sum)
  colnames(intersections_sum) <- c("Intersection" ,"number_of_DE_genes")
  
  intersections_name <- names(intersections) %>% as.list()
  ## convert the list into a character vector
  Gene_ID <- character()
  for (i in seq_along(intersections)){
    for(j in seq_along(intersections[[i]])){
      Gene_ID <- c(Gene_ID, intersections[[i]][j])
    }
  }
  
  intersect_label <- map2(intersections_name, intersections_length, rep)
  intersect_label_vecter <- character()
  for (i in seq_along(intersect_label)){
    for(j in seq_along(intersect_label[[i]])){
      intersect_label_vecter <- c(intersect_label_vecter, intersect_label[[i]][j])
    }
  }
  gene_intersect <- data.frame(Gene_ID, intersect_label_vecter) 
  colnames(gene_intersect) <- c("Orig_Symbol", "Intersections")
  
  ### grab all the q and FC from all comparisons
  ## df is the all in one file, just full join GWAS list with P-value
  df <- df %>% full_join(GWAS, by=c("Orig_Symbol"="Symbol"))
  
  gene_intersect <- left_join(gene_intersect, df, by=c("Orig_Symbol"))
  write_delim(gene_intersect, paste(output, "cluster_", cluster, "_vennlist_vsB6", ".txt", sep=""), delim = "\t")
  
  gene_intersect_GWAS <-  gene_intersect %>% filter(str_detect(Intersections, "\\:GWAS"))
  write_delim(gene_intersect_GWAS, paste(output, "cluster_", cluster, "_vennlist_vsB6_GWAS", ".txt", sep=""), delim = "\t")
  
}

