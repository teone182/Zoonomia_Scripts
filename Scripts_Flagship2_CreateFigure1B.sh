#!/bin/bash -l

#####TopMed files with all SNVs were available for every chromosome (.Rdata). Preliminary steps involved intersection with phyloP scores####
####header -> chr end AC AF phyloPam###5 tab separated columns with chromosome, position, allele count, allele frequency, phyloP score####

WorkFolder='Path/To/Working/Folder'
inFolder='Path/To/TopMed/Starting/Files'
PCspacesFolder='Path/To/Genomic/Annotations'

cd $inFolder

R --vanilla <<EOF

library("dplyr")
library("ggplot2")

file_list <- list.files("$inFolder", pattern=glob2rx("TopMedv8_chr*.RData"), full.names = T)
file_list
total <- do.call(rbind, lapply(file_list, readRDS)

total <- total %>%
  mutate(
    Category = case_when(
      AC == 1 ~ "AC=1",
      AC == 2 ~ "AC=2",
      AC >= 3 & AC <= 5 ~ "AC3-5",
      AC >= 6 & AC <= 10 ~ "AC6-10",
      AC >=11 & AF < 0.005 ~ "AC>=11AF<0.005",
      AF >= 0.005 ~ "AF>=0.005"
    )
 )

data_2export <- total %>% select(chr, end, phyloPam, Category)
write.table(data_2export, "$WorkFolder/TopMed_Vars_ALLChromosomes_WithCategories_4Plot_4Paper", quote=F, sep="\t", row.names=F, col.names=F)

EOF


cd $WorkFolder

#####create bed file for intersection#####
awk -v OFS="\t" '{print $1,$2-1,$2,$3,$4,$5}' TopMed_Vars_ALLChromosomes_WithCategories_4Plot_4Paper > TopMed_Vars_ALLChromosomes_WithCategories_4Plot_4Paper.bed

####overlap with cds (pc)######
bedtools intersect -a TopMed_Vars_ALLChromosomes_WithCategories_4Plot_4Paper.bed -b $PCspacesFolder/pc_pickOne_cds -wo | awk -v OFS="\t" '{print $1,$3,$4,$5,$6}' > TopMed_Vars_ALLChromosomes_IntersectedWith_PC_cds    ##We refer to the pc_pickOne_cds annotation#####   

awk '$4 == "AC=1" {print}' TopMed_Vars_ALLChromosomes_IntersectedWith_PC_cds | sort -n -k3,3 | awk -f Calculate_Median.awk | awk -v OFS="" '{print "a","\t",$1,"\t","cds","\t",1}' > TopMed_Vars_ALLChromosomes_PC_cds_Category_AC1_MedianPhyloP
awk '$4 == "AC=2" {print}' TopMed_Vars_ALLChromosomes_IntersectedWith_PC_cds | sort -n -k3,3 | awk -f Calculate_Median.awk | awk -v OFS="" '{print "b","\t",$1,"\t","cds","\t",2}' > TopMed_Vars_ALLChromosomes_PC_cds_Category_AC2_MedianPhyloP
awk '$4 == "AC3-5" {print}' TopMed_Vars_ALLChromosomes_IntersectedWith_PC_cds | sort -n -k3,3 | awk -f Calculate_Median.awk | awk -v OFS="" '{print "c","\t",$1,"\t","cds","\t",3}' > TopMed_Vars_ALLChromosomes_PC_cds_Category_AC3-5_MedianPhyloP
awk '$4 == "AC6-10" {print}' TopMed_Vars_ALLChromosomes_IntersectedWith_PC_cds | sort -n -k3,3 | awk -f Calculate_Median.awk | awk -v OFS="" '{print "d","\t",$1,"\t","cds","\t",4}' > TopMed_Vars_ALLChromosomes_PC_cds_Category_AC6-10_MedianPhyloP
awk '$4 == "AC>=11AF<0.005" {print}' TopMed_Vars_ALLChromosomes_IntersectedWith_PC_cds | sort -n -k3,3 | awk -f Calculate_Median.awk | awk -v OFS="" '{print "e","\t",$1,"\t","cds","\t",5}' > TopMed_Vars_ALLChromosomes_PC_cds_Category_AC11AF0005_MedianPhyloP
awk '$4 == "AF>=0.005" {print}' TopMed_Vars_ALLChromosomes_IntersectedWith_PC_cds | sort -n -k3,3 | awk -f Calculate_Median.awk | awk -v OFS="" '{print "f","\t",$1,"\t","cds","\t",6}' > TopMed_Vars_ALLChromosomes_PC_cds_Category_AF0005_MedianPhyloP

#####overlap with non-coding ######get all the vars overlapping with cds and grep out them from the total set of variants######
grep -vFf TopMed_Vars_ALLChromosomes_IntersectedWith_PC_cds TopMed_Vars_ALLChromosomes_WithCategories_4Plot_4Paper > TopMed_Vars_ALLChromosomes_IntersectedWith_PC_non_coding

awk '$4 == "AC=1" {print}' TopMed_Vars_ALLChromosomes_IntersectedWith_PC_non_coding | sort -n -k3,3 | awk -f Calculate_Median.awk | awk -v OFS="" '{print "a","\t",$1,"\t","non_cds","\t",1}' > TopMed_Vars_ALLChromosomes_PC_non_coding_Category_AC1_MedianPhyloP
awk '$4 == "AC=2" {print}' TopMed_Vars_ALLChromosomes_IntersectedWith_PC_non_coding | sort -n -k3,3 | awk -f Calculate_Median.awk | awk -v OFS="" '{print "b","\t",$1,"\t","non_cds","\t",2}' > TopMed_Vars_ALLChromosomes_PC_non_coding_Category_AC2_MedianPhyloP
awk '$4 == "AC3-5" {print}' TopMed_Vars_ALLChromosomes_IntersectedWith_PC_non_coding | sort -n -k3,3 | awk -f Calculate_Median.awk | awk -v OFS="" '{print "c","\t",$1,"\t","non_cds","\t",3}' > TopMed_Vars_ALLChromosomes_PC_non_coding_Category_AC3-5_MedianPhyloP
awk '$4 == "AC6-10" {print}' TopMed_Vars_ALLChromosomes_IntersectedWith_PC_non_coding | sort -n -k3,3 | awk -f Calculate_Median.awk | awk -v OFS="" '{print "d","\t",$1,"\t","non_cds","\t",4}' > TopMed_Vars_ALLChromosomes_PC_non_coding_Category_AC6-10_MedianPhyloP
awk '$4 == "AC>=11AF<0.005" {print}' TopMed_Vars_ALLChromosomes_IntersectedWith_PC_non_coding | sort -n -k3,3 | awk -f Calculate_Median.awk | awk -v OFS="" '{print "e","\t",$1,"\t","non_cds","\t",5}' > TopMed_Vars_ALLChromosomes_PC_non_coding_Category_AC11AF0005_MedianPhyloP
awk '$4 == "AF>=0.005" {print}' TopMed_Vars_ALLChromosomes_IntersectedWith_PC_non_coding | sort -n -k3,3 | awk -f Calculate_Median.awk | awk -v OFS="" '{print "f","\t",$1,"\t","non_cds","\t",6}' > TopMed_Vars_ALLChromosomes_PC_non_coding_Category_AF0005_MedianPhyloP

cat TopMed_Vars_ALLChromosomes_PC_cds_Category_AC1_MedianPhyloP TopMed_Vars_ALLChromosomes_PC_non_coding_Category_AC1_MedianPhyloP TopMed_Vars_ALLChromosomes_PC_cds_Category_AC2_MedianPhyloP TopMed_Vars_ALLChromosomes_PC_non_coding_Category_AC2_MedianPhyloP TopMed_Vars_ALLChromosomes_PC_cds_Category_AC3-5_MedianPhyloP TopMed_Vars_ALLChromosomes_PC_non_coding_Category_AC3-5_MedianPhyloP TopMed_Vars_ALLChromosomes_PC_cds_Category_AC6-10_MedianPhyloP TopMed_Vars_ALLChromosomes_PC_non_coding_Category_AC6-10_MedianPhyloP TopMed_Vars_ALLChromosomes_PC_cds_Category_AC11AF0005_MedianPhyloP TopMed_Vars_ALLChromosomes_PC_non_coding_Category_AC11AF0005_MedianPhyloP TopMed_Vars_ALLChromosomes_PC_cds_Category_AF0005_MedianPhyloP TopMed_Vars_ALLChromosomes_PC_non_coding_Category_AF0005_MedianPhyloP > TopMed_Vars_ALLChromosomes_PC_ANNOTATION_FILE_MEDIANS_minimal

####create a file with phyloP scores and AC/AF classes amenable for plotting (essential info)
awk -v OFS="\t" '{print $3,$4}' TopMed_Vars_ALLChromosomes_WithCategories_4Plot_4Paper | sed 's/AC=1/a/g' | sed 's/AC=2/b/g' | sed 's/AC3-5/c/g' | sed 's/AC6-10/d/g' | sed 's/AC>=11AF<0\.005/e/g' | sed 's/AF>=0\.005/f/g' > TopMed_Vars_ALLChromosomes_WithCategories_4Plot_4Paper_minimal

R --vanilla <<EOF

library(dplyr)
library(ggplot2)

data <- read.csv("TopMed_Vars_ALLChromosomes_WithCategories_4Plot_4Paper_minimal", head=F, sep="\t")
data\$V2 <- factor(data\$V2, levels=c("a","b","c","d","e","f"))

data <- data[complete.cases(data), ]

anno <- read.csv("TopMed_Vars_ALLChromosomes_PC_ANNOTATION_FILE_MEDIANS_minimal", head=F, sep="\t")
anno\$V1 <- factor(anno\$V1, levels=c("a","b","c","d","e","f"))

pdf("Figure_1B_BiggerFontSize.pdf")
ggplot(data = data, aes(x=V2, y=V1)) +
  geom_boxplot(outlier.shape = NA, color="grey", fill="grey", alpha=0.4) + coord_cartesian(ylim=c(-3,4)) +
  ylab("phyloP score") +
  xlab("") +
  theme_classic() +
  scale_y_continuous(breaks = seq(-20, 10, by = 2)) +
  scale_x_discrete(breaks=c("a","b","c","d","e","f"),
                   labels=c("a" = "AC=1", "b" = "AC=2", "c" = "AC 3-5", "d" = "AC 6-10",
                            "e" = expression("AC">="11 AF<0.005"), "f" = expression("AF">="0.005"))) +
  theme(axis.text.x = element_text(angle=75, hjust=1, size=14), axis.title.y = element_text(size=16), axis.text.y = element_text(size=16),
       legend.title=element_text(size=14), legend.text = element_text(size=14)) +
  geom_point(data = anno, aes(x=V1, y=V2, colour=factor(V3)), size=3) +
  geom_line(data = anno, aes(x=V4, y=V2, colour=factor(V3)), size=0.8) +
  scale_color_manual(name="Category",
                     labels=c("Median phyloP cds","Median phyloP non-cds"),
                     values=c("firebrick4","orange"))
dev.off()
EOF

