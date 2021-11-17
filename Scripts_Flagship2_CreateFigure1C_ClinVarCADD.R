library(dplyr)
library(ggplot2)

setwd("Path/")
ClinVar <- read.csv("ClinVarHg38_AllClasses_IntersectedWithPhyloP_AllChromosomes.bed", sep="\t", header = T)

#####Create new phyloP for collapsing data < -10
ClinVar$New_PhyloP <- ClinVar$PhyloP
ClinVar$New_PhyloP[ClinVar$New_PhyloP < -10] <- -10


CADD <- read.csv("Variants_WithCADDGr20_WithHeader.bed", sep="\t", header = T)   ####We have to add a header to the Variants_WithCADDGr20.bed file (already generated) --> "chr\tstart\tend\tPhyloP\tgene_name\tClass\tOrigin"
    
#####Create new phyloP for collapsing data < -10
CADD$New_PhyloP <- CADD$PhyloP
CADD$New_PhyloP[CADD$New_PhyloP < -10] <- -10

data <- rbind(ClinVar, CADD)

data <- data %>% 
  mutate(
    "Category" = case_when(
      Class == "Pathogenic" | Class == "Likely_pathogenic" ~ "ClinVar Pathogenic",
      Class == "Benign" | Class == "Likely_benign" ~ "ClinVar Benign",
      Class == "cadd>=20" ~ "CADD>=20"
    )
  )

data$Category <- factor(data$Category, levels=c("ClinVar Benign", "ClinVar Pathogenic", "CADD>=20"))

pdf("Figure_1C.pdf")
ggplot(data, aes(x=Category, y=New_PhyloP, fill=Category)) + 
  geom_violin() + scale_fill_manual(values=c("dodgerblue4", "firebrick4","orangered2")) +
  ylab("phyloP score") +
  theme_classic() +
  xlab("") +
  theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), legend.position = "none") +
  scale_y_continuous(breaks = seq(-10, 10, by = 2)) +
  scale_x_discrete(breaks=c("ClinVar Benign","ClinVar Pathogenic","CADD>=20"), labels=c("ClinVar Benign" = "ClinVar benign", "ClinVar Pathogenic" = "ClinVar pathogenic", "CADD>=20" = expression("CADD">="20")))

  dev.off()


