#!/bin/bash -l

InFolder='Path/To/Start/file'
WorkFolder='Path'
AnnoFolder='Path/To/All/Gencode/Annotations/and/Cadd/values'  ###this dir contains Gencode and other regulatory annotation bed files#####
OutFolder='Path/To/Summary/files/and/plots'

##########
cd $WorkFolder

while read n
do

cd $WorkFolder

######All PhyloP >= 2.27 only##########
awk -v var="${n}" -v FS="\t" -v OFS="\t" 'FNR == 1 {next} $3 >= 2.27 && $4 < var {print $1,$2-1,$2,$3,$4}' $InFolder/MERGED_PHYLOP_CADD_ALL_CHROMOSOMES.txt > $WorkFolder/Variants_WithPhyloPFDR05_CADDLw${n}.bed &&

######All CADD >= 10-20##########
awk -v var="${n}" -v FS="\t" -v OFS="\t" 'FNR == 1 {next} $4 >= var && $3 < 2.27 {print $1,$2-1,$2,$3,$4}' $InFolder/MERGED_PHYLOP_CADD_ALL_CHROMOSOMES.txt > $WorkFolder/Variants_WithCADDGr${n}_WithPhyloPLwFDR05.bed &&

######All PhyloP >= 2.27 and CADD >=10-20##########
awk -v var="${n}" -v FS="\t" -v OFS="\t" 'FNR == 1 {next} $3 >= 2.27 && $4 >= var {print $1,$2-1,$2,$3,$4}' $InFolder/MERGED_PHYLOP_CADD_ALL_CHROMOSOMES.txt > $WorkFolder/Variants_WithPhyloPFDR05_and_CADDGr${n}.bed

done<$AnnoFolder/List_CADD_values.list  #####list of values from 10 to 20#####



cd $WorkFolder

while read n
do
while read class
do

cd $WorkFolder

#####################
###PhyloP

Nr_vars_class=$(cat $AnnoFolder/${class}.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}') &&

bedtools intersect -a $AnnoFolder/${class}.bed -b $WorkFolder/Variants_WithPhyloPFDR05_CADDLw${n}.bed -wo | awk -v OFS="\t" '{print $5,$7,$8,$9,$4}' > ${class}_PhyloP_FDR05_CADDLw${n} &&

Nr_ThisClass_PhyloPyes=$(wc -l "${class}_PhyloP_FDR05_CADDLw${n}" | cut -f1 -d' ') &&

if [[ $Nr_ThisClass_PhyloPyes -eq 0 ]]
then
Fraction_PhyloPyes_ThisClass=0
else
Fraction_PhyloPyes_ThisClass=$(echo "scale=3; $Nr_ThisClass_PhyloPyes/$Nr_vars_class" | bc)
fi

echo -e "${class}\t$Fraction_PhyloPyes_ThisClass\tPhyloP>=2.27_CADD<${n}" > Fraction_PhyloPYES_CADDLw${n}_${class}.txt &&



#####################
###CADD
Nr_vars_class=$(cat $AnnoFolder/${class}.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}') &&

bedtools intersect -a $AnnoFolder/${class}.bed -b $WorkFolder/Variants_WithCADDGr${n}_WithPhyloPLwFDR05.bed -wo | awk -v OFS="\t" '{print $5,$7,$8,$9,$4}' > ${class}_CADD_Gr${n}_PhyloPLwFDR05 &&

Nr_ThisClass_CADDyes=$(wc -l "${class}_CADD_Gr${n}_PhyloPLwFDR05" | cut -f1 -d' ') &&

if [[ $Nr_ThisClass_CADDyes -eq 0 ]]
then
Fraction_CADDyes_ThisClass=0
else
Fraction_CADDyes_ThisClass=$(echo "scale=3; $Nr_ThisClass_CADDyes/$Nr_vars_class" | bc)
fi

echo -e "${class}\t$Fraction_CADDyes_ThisClass\tCADD>=${n}_PhyloP<2.27" > Fraction_CADDyes_CADD${n}_${class}.txt &&



#####################
###Both PhyloP and CADD

Nr_vars_class=$(cat $AnnoFolder/${class}.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}') &&

bedtools intersect -a $AnnoFolder/${class}.bed -b $WorkFolder/Variants_WithPhyloPFDR05_and_CADDGr${n}.bed -wo | awk -v OFS="\t" '{print $5,$7,$8,$9,$4}' > ${class}_PhyloP_FDR05_CADD_Gr${n} &&

Nr_ThisClass_PhyloPyes_CADDyes=$(wc -l "${class}_PhyloP_FDR05_CADD_Gr${n}" | cut -f1 -d' ') &&

if [[ $Nr_ThisClass_PhyloPyes_CADDyes -eq 0 ]]
then
Fraction_PhyloPyes_CADDyes_ThisClass=0
else
Fraction_PhyloPyes_CADDyes_ThisClass=$(echo "scale=3; $Nr_ThisClass_PhyloPyes_CADDyes/$Nr_vars_class" | bc)
fi

echo -e "${class}\t$Fraction_PhyloPyes_CADDyes_ThisClass\tPhyloP>=2.27_CADD>=${n}" > Fraction_PhyloPyes_CADDyes_CADD${n}_${class}.txt

done<$AnnoFolder/List_ALL_Genomic_Categories.list   #####list of annotation categories (e.g pc_pickOne_cds, pc_pickAll_intron, etc.#####   
done<$AnnoFolder/List_CADD_values.list      #####list of values from 10 to 20#####   


##########
cd $WorkFolder

while read n
do
while read class
do

cd $WorkFolder

sed "s/$/\t${n}\tOnly PhyloP constrained variant positions/g" Fraction_PhyloPYES_CADDLw${n}_${class}.txt > $OutFolder/Fraction_PhyloPYES_CADDLw${n}_${class}.txt &&
sed "s/$/\t${n}\tPhyloP constrained and CADD score variant positions/g" Fraction_PhyloPyes_CADDyes_CADD${n}_${class}.txt > $OutFolder/Fraction_PhyloPyes_CADDyes_CADD${n}_${class}.txt &&
sed "s/$/\t${n}\tOnly CADD score variant positions/g" Fraction_CADDyes_CADD${n}_${class}.txt > $OutFolder/Fraction_CADDyes_CADD${n}_${class}.txt

done<$AnnoFolder/List_ALL_Genomic_Categories.list    #####list of annotation categories (e.g pc_pickOne_cds, pc_pickAll_intron, etc.#####
done<$AnnoFolder/List_CADD_values.list    #####list of values from 10 to 20#####  

#####summarise####
cd $OutFolder

while read class
do

cd $OutFolder

cat Fraction_*_${class}.txt | sed "1iAnnotation\tProportion\tMeasure\tValue\tCategory" > FRACTION_ALL_CATEGORIES_${class}.txt &&
rm Fraction_*_${class}.txt

done<$AnnoFolder/List_ALL_Genomic_Categories.list    #####list of annotation categories (e.g pc_pickOne_cds, pc_pickAll_intron, etc.#####



#####Plot FigureS5_C#######
####Create a dataframe with number of variants for only phyloP, only CADD and PhyloP-CADD####
cd $WorkFolder

while read n
do

cd $WorkFolder

#####phyloP
NrVars_OnlyPhyloP=$(wc -l Variants_WithPhyloPFDR05_CADDLw${n}.bed | cut -f1 -d' ')
echo -e "${n}\tPhyloP>=2.27\t$NrVars_OnlyPhyloP\tOnly PhyloP constrained variant positions" > Tmp_OnlyPhyloPFDR05_CADDLw${n}

####cadd####
NrVars_OnlyCADD=$(wc -l Variants_WithCADDGr${n}_WithPhyloPLwFDR05.bed | cut -f1 -d' ')
echo -e "${n}\tCADD>=${n}\t$NrVars_OnlyCADD\tOnly CADD score variant positions" > Tmp_OnlyCADDGr${n}_PhyloPLwFDR05

#####phyloP and cadd######
NrVars_PhyloP_and_CADD=$(wc -l Variants_WithPhyloPFDR05_and_CADDGr${n}.bed | cut -f1 -d' ')
echo -e "${n}\tPhyloP>=2.27_CADD>=${n}\t$NrVars_PhyloP_and_CADD\tPhyloP constrained and CADD score variant positions" > Tmp_PhyloPFDR05_CADDGr${n}

done<$AnnoFolder/List_CADD_values.list      #####list of values from 10 to 20#####

cat Tmp_OnlyPhyloPFDR05_CADDLw* Tmp_OnlyCADDGr*_PhyloPLwFDR05 Tmp_PhyloPFDR05_CADDGr* | sed '1iCutoff\tMeasure\tVariants\tCategory' > $OutFolder/Table_NrVariants_UNIQUE.txt
rm Tmp_*


R --vanilla <<EOF

library("dplyr")
library("ggplot2")

data <- read.csv("$OutFolder/Table_NrVariants_UNIQUE.txt", head=T, sep="\t")
pdf("$OutFolder/Figure5S-C.pdf")
data$Category <- factor(data$Category, levels = c("Only CADD score variant positions","PhyloP constrained and CADD score variant positions","Only PhyloP constrained variant positions"))
ggplot(data, aes(x=Cutoff, y=Variants, fill=Category)) +
  geom_bar(stat="identity") +
  scale_x_continuous(breaks = seq(10, 20, by = 1), labels=c("10" = expression("">="10"), "11" = expression("" >="11"), "12" = expression("" >="12"), "13" = expression("" >="13"), "14" = expression("" >="14"), "15" = expression("" >="15"),
                                                              "16" = expression("" >="16"),"17" = expression("" >="17"),"18" = expression("" >="18"),"19" = expression("" >="19"),"20" = expression("" >="20"))) +
  ylab("Number of variant positions (million)") +
  xlab("CADD score threshold") +
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=75, hjust=1)) +
  scale_fill_manual(values = c("firebrick4","#619CFF","darkblue"))
dev.off()

EOF




#################
######Merge the genomic annotation of interest and then plot (FigureS5_D)###########

R --vanilla <<EOF

library("dplyr")
library("ggplot2")

wanted_anno_pc_PickOne <- list.files("$OutFolder", pattern=glob2rx("FRACTION_ALL_CATEGORIES_pc_pickOne_*.txt"), full.names = T) ###here pc_pickOne
wanted_anno_reg <- list.files("$OutFolder", pattern=glob2rx("FRACTION_ALL_CATEGORIES_reg_*.txt"), full.names = T) ###here reg

file_list <- c(wanted_anno_pc_PickOne,wanted_anno_reg)
file_list

data <- do.call(rbind, lapply(file_list, read.csv, header = T, sep="\t"))

data$Proportion_PerCent <- data$Proportion * 100
data$Category <- factor(data$Category, levels = c("Only CADD score variant positions","PhyloP constrained and CADD score variant positions","Only PhyloP constrained variant positions"))
pdf("$OutFolder/Figure_S5-D.pdf", width = 15 , height = 15)
ggplot(data, aes(x=Value, y=Proportion_PerCent)) + geom_point(aes(color=Category)) + geom_line(aes(color=Category)) +
  ylab("Percentage of category in annotation (%)") +
  xlab("CADD score threshold") +
  scale_x_continuous(breaks=c(10,11,12,13,14,15,16,17,18,19,20),
                     labels=c("10" = expression("">="10"), "11" = expression("">="11"), "12" = expression("">="12"), "13" = expression("">="13"), "14" = expression("">="14"), "15" = expression("">="15"),
                              "16" = expression("">="16"),"17" = expression("">="17"),"18" = expression("">="18"),"19" = expression("">="19"),"20" = expression("">="20"))) +
  #   scale_x_continuous(breaks=c(10,11,12,13,14,15,16,17,18,19,20)) +
  facet_wrap(~Annotation) +
  theme(strip.text = element_text(size=10), axis.text.x = element_text(size = 10, angle=75, hjust=1),
        axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.title=element_text(size=16), legend.text=element_text(size=16)) +
  scale_color_manual(values = c("firebrick4","#619CFF","darkblue"))
dev.off()

EOF
