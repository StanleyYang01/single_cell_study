# making Venn diagram to show the intersection of the genes from strain effect: WSB, PWK, CAST vs B6 and how they overlap with GWAS

library(tidyverse)
library(ggrepel)

output <- "./04_volcano_strain_vsB6/"

clusters <- c("H", "6", "7", "8", "9", "10", "11", "12")
vennlist <- paste("cluster_", clusters[1], "_vennlist_vsB6", ".txt", sep="")
vennlist <- file.path("./03_venn_strain_vsB6", vennlist)

strain_colors <- c("#888888", "#00AA00", "#FF0000", "#9900EE") %>% as.list()
names(strain_colors) <- c("B6", "CAST", "PWK", "WSB")
strains <- c("CAST", "PWK", "WSB")



df <- read_delim(vennlist[1], delim = "\t", 
                 col_types = "ccccdcddddddcddddddddcd")

#c = character, i = integer, n = number, d = double, l = logical, f = factor, D = date, T = date time, t = time, ? = guess, or _/- to skip the column

# set the limits for the plot for a given cluster
x_min <-df %>% select(contains("logFC")) %>% map(min, na.rm=TRUE) %>% unlist() %>% min()
x_max <-df %>% select(contains("logFC")) %>% map(max, na.rm=TRUE) %>% unlist() %>% max()
y_max <-df %>% select(contains("FDR")) %>% mutate_all(~-log10(.)) %>% map(max, na.rm=TRUE) %>% unlist() %>% max()

x_GWAS_min <-df %>% drop_na() %>% select(contains("logFC")) %>% map(min, na.rm=TRUE) %>% unlist() %>% min()
x_GWAS_max <-df %>% drop_na() %>% select(contains("logFC")) %>% map(max, na.rm=TRUE) %>% unlist() %>% max()
y_GWAS_max <-df %>% drop_na() %>% select(contains("FDR")) %>% mutate_all(~-log10(.)) %>% map(max, na.rm=TRUE) %>% unlist() %>% max()

# do volcano plot for CAST vs B6
df_subset <- df %>% select(Orig_Symbol, contains(strains[2]),PUKB_IGAP, Intersections) 
names(df_subset)[2:3] <- c("logFC", "FDR")



df_subset <- df_subset %>% 
  mutate(logFDR= -log10(FDR))
         #mark_Strain = ifelse(logFDR>12 & abs(logFC)>10, "Sig", "NotSig"),
         #ark_GWAS = ifelse(logFDR>5 & abs(logFC)>3 & PUKB_IGAP<0.02, "Sig", "NotSig"))

top_gene_n <- 25

# plot gene of strain difference
ggplot(df_subset, aes(x=logFC, y=logFDR))+
  geom_point(alpha=0.5, size=0.5, color="grey")+
  geom_point(data = df_subset %>% filter(logFDR>10) %>% top_n(top_gene_n, abs(logFC)),
             size=1, color=strain_colors[[strains[2]]])+
  scale_x_continuous(limits = c(x_min, x_max))+
  scale_y_continuous(limits = c(NA, y_max))+
  geom_text_repel(data = df_subset %>% filter(logFDR>10) %>% top_n(top_gene_n, abs(logFC)), 
                  aes(label = Orig_Symbol),
                  size = 3,
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.1, "lines"))+
  coord_fixed()+
  labs(x="log2FC", y="-log10(FDR)")+
  ggtitle(paste("cluster", clusters[1], strains[2], "vs", "B6", sep = " "))+
  theme_bw()
ggsave(paste(output, "cluster_", clusters[1], "_volcano_", strains[2], "vsB6.png", sep=""), 
       width = 4.5, height = 4, units = "in", dpi = 300)


# plot gene of strain difference overlap with GWAS
ggplot(df_subset, aes(x=logFC, y=logFDR))+
  geom_point(alpha=0.5, size=0.5, color="grey")+
  geom_point(data = df_subset %>% 
               filter(logFDR>5 & PUKB_IGAP<0.02) %>% top_n(top_gene_n, abs(logFC)),
             size=1, color=strain_colors[[strains[2]]])+
  scale_x_continuous(limits = c(x_GWAS_min, x_GWAS_max))+
  scale_y_continuous(limits = c(NA, y_GWAS_max))+
  geom_text_repel(data = df_subset %>% 
                    filter(logFDR>5 & PUKB_IGAP<0.02) %>% top_n(top_gene_n, abs(logFC)), 
                  aes(label = Orig_Symbol),
                  size = 3,
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.1, "lines"))+
  coord_fixed()+
  labs(x="log2FC", y="-log10(FDR)")+
  ggtitle(paste("cluster", clusters[1], strains[2], "vs", "B6", "GWAS", sep = " "))+
  theme_bw()
ggsave(paste(output, "cluster_", clusters[1], "_volcano_", strains[2], "vsB6_GWAS.png", sep=""), 
       width = 4.5, height = 4, units = "in", dpi = 300)


