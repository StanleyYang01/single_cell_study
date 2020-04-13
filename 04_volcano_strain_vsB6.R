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

# do volcano plot for CAST vs B6
df_subset <- df %>% select(Orig_Symbol, contains(strains[1]),PUKB_IGAP, Intersections) 
names(df_subset)[2:3] <- c("logFC", "FDR")

df_subset <- df_subset %>% 
  mutate(logFDR= -log10(FDR),
         mark = ifelse(logFDR>5 & abs(logFC)>3 & PUKB_IGAP<0.02, "Sig", "NotSig"))

# df_subset$mark %>% table() # check the number of significant genes

ggplot(df_subset, aes(x=logFC, y=logFDR))+
  geom_point(alpha=0.5, size=0.5)+
  geom_point(data = df_subset %>% filter(mark=="Sig"), size=1, color=strain_colors[[strains[1]]])+
  geom_text_repel(data = df_subset %>% filter(mark=="Sig"), 
                  aes(label = Orig_Symbol),
                  size = 3,
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.1, "lines"),
                  color="white")+
  coord_fixed()+
  labs(x="log2FC", y="-log10(FDR)")+
  ggtitle(paste("cluster", clusters[1], strains[1], "vs", "B6", sep = " "))+
  theme_dark()

ggsave(paste(output, "cluster_", clusters[1], "_volcano_", strains[1], "vsB6.png", sep=""), 
       width = 4.5, height = 4, units = "in", dpi = 300)
