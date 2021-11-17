#!/bin/bash -l

InFolder='Path/To/Start/file'
WorkFolder='Path'
AnnoFolder='Path/To/All/Gencode/Annotations/and/Cadd/values'  ###this dir contains Gencode and other regulatory annotation bed files#####
OutFolder='Path/To/Summary/files/and/plots'

##########
cd $WorkFolder

######All PhyloP >= 2.27##########
awk -v FS="\t" -v OFS="\t" 'FNR == 1 {next} $3 >= 2.27 {print $1,$2-1,$2,$3,$4}' $InFolder/MERGED_PHYLOP_CADD_ALL_CHROMOSOMES.txt > $WorkFolder/Variants_WithPhyloPFDR05.bed &&

cd $WorkFolder
while read n
do

cd $WorkFolder

######All CADD >= 10-20##########
awk -v var="${n}" -v FS="\t" -v OFS="\t" 'FNR == 1 {next} $4 >= var {print $1,$2-1,$2,$3,$4}' $InFolder/MERGED_PHYLOP_CADD_ALL_CHROMOSOMES.txt > $WorkFolder/Variants_WithCADDGr${n}.bed

done<$AnnoFolder/List_CADD_values.list   #####list of values from 10 to 20#####

while read n
do
while read class
do

cd $WorkFolder

#####################
###PhyloP

Nr_vars_class=$(cat $AnnoFolder/${class}.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}') &&

bedtools intersect -a $AnnoFolder/${class}.bed -b Variants_WithPhyloPFDR05.bed -wo | awk -v OFS="\t" '{print $5,$7,$8,$9,$4}' > ${class}_PhyloP_FDR05 &&

Nr_ThisClass_PhyloPyes=$(wc -l "${class}_PhyloP_FDR05" | cut -f1 -d' ') &&

if [[ $Nr_ThisClass_PhyloPyes -eq 0 ]]
then
Fraction_PhyloPyes_ThisClass=0
else
Fraction_PhyloPyes_ThisClass=$(echo "scale=3; $Nr_ThisClass_PhyloPyes/$Nr_vars_class" | bc)
fi

echo -e "${class}\t$Fraction_PhyloPyes_ThisClass\tPhyloP>=2.27" > Fraction_PhyloPYES_Nr${n}_${class}.txt &&


#####################
###CADD
Nr_vars_class=$(cat $AnnoFolder/${class}.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}') &&

bedtools intersect -a $AnnoFolder/${class}.bed -b Variants_WithCADDGr${n}.bed -wo | awk -v OFS="\t" '{print $5,$7,$8,$9,$4}' > ${class}_CADD_Gr${n} &&

Nr_ThisClass_CADDyes=$(wc -l "${class}_CADD_Gr${n}" | cut -f1 -d' ') &&

if [[ $Nr_ThisClass_CADDyes -eq 0 ]]
then
Fraction_CADDyes_ThisClass=0
else
Fraction_CADDyes_ThisClass=$(echo "scale=3; $Nr_ThisClass_CADDyes/$Nr_vars_class" | bc)
fi

echo -e "${class}\t$Fraction_CADDyes_ThisClass\tCADD>=${n}" > Fraction_CADDyes_CADD${n}_${class}.txt


####
cd $WorkFolder

sed "s/$/\t${n}\tPhyloP constrained variant positions/g" Fraction_PhyloPYES_Nr${n}_${class}.txt > $OutFolder/Fraction_PhyloPYES_Nr${n}_${class}.txt &&
sed "s/$/\t${n}\tCADD score variant positions/g" Fraction_CADDyes_CADD${n}_${class}.txt > $OutFolder/Fraction_CADDyes_CADD${n}_${class}.txt

done<$AnnoFolder/List_ALL_Genomic_Categories.list  #####list of annotation categories (e.g pc_pickOne_cds, pc_pickAll_intron, etc.#####
done<$AnnoFolder/List_CADD_values.list      #####list of values from 10 to 20#####


####summarise

cd $OutFolder

while read class
do

cd $OutFolder

cat Fraction_*_${class}.txt | sed "1iAnnotation\tProportion\tMeasure\tValue\tCategory" > FRACTION_ALL_CATEGORIES_${class}.txt

done<$AnnoFolder/List_ALL_Genomic_Categories.list #####list of annotation categories (e.g pc_pickOne_cds, pc_pickAll_intron, etc.#####


#####Plot FigureS5_1_A#######
####Create a dataframe with number of variants for phyloP constraint and every CADD class####
cd $WorkFolder

while read n
do

cd $WorkFolder

#####phyloP
NrVars_PhyloPConstraint=$(wc -l Variants_WithPhyloPFDR05.bed | cut -f1 -d' ')
echo -e "${n}\tPhyloP>=2.27\t$NrVars_PhyloPConstraint\tPhyloP constrained variant positions" > Tmp_NrPhyloPConstraint_Nr${n}

####cadd####
NrVars_CADD=$(wc -l Variants_WithCADDGr${n}.bed | cut -f1 -d' ')
echo -e "${n}\tCADD>=${n}\t$NrVars_CADD\tCADD score variant positions" > Tmp_NrCADDGr${n}

done<$AnnoFolder/List_CADD_values.list      #####list of values from 10 to 20#####

cat Tmp_NrPhyloPConstraint_Nr* Tmp_NrCADDGr* | sed '1iCutoff\tMeasure\tVariants\tCategory' > $OutFolder/Table_NrVariants_TOTAL.txt
rm Tmp_*

#####Plot###
R --vanilla <<EOF
library("dplyr")
library("ggplot2")

data <- read.csv("$OutFolder/Table_NrVariants_TOTAL.txt", head=T, sep="\t")
head(data)
pdf("$OutFolder/FigureS5_1-A.pdf")
ggplot(data, aes(x=Cutoff, y=Variants, fill=Category)) +
  geom_bar(stat="identity") +
  scale_x_continuous(breaks = seq(10, 20, by = 1), labels=c("10" = expression("">="10"), "11" = expression("" >="11"), "12" = expression("" >="12"), "13" = expression("" >="13"), "14" = expression("" >="14"), "15" = expression("" >="15"),
                                                            "16" = expression("" >="16"),"17" = expression("" >="17"),"18" = expression("" >="18"),"19" = expression("" >="19"),"20" = expression("" >="20"))) +
  ylab("Number of variant positions (million)") +
  xlab("CADD score threshold") +
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=75, hjust=1)) +
scale_fill_manual(values=c("firebrick4","darkblue"))
dev.off()
EOF


#################
######Merge the genomic annotation of interest and then plot (FigureS5_1_B)###########

R --vanilla <<EOF
library("dplyr")
library("ggplot2")

wanted_anno_pc_PickOne <- list.files("$OutFolder", pattern=glob2rx("FRACTION_ALL_CATEGORIES_pc_pickOne_*.txt"), full.names = T) ###here pc_pickOne
wanted_anno_reg <- list.files("$OutFolder", pattern=glob2rx("FRACTION_ALL_CATEGORIES_reg_*.txt"), full.names = T) ###here reg

file_list <- c(wanted_anno_pc_PickOne,wanted_anno_reg)
file_list

data <- do.call(rbind, lapply(file_list, read.csv, header = T, sep="\t"))

data$Proportion_PerCent <- data$Proportion * 100

pdf("$OutFolder/Figure_S5_1-B.pdf", width = 15 , height = 15)
ggplot(data, aes(x=Value, y=Proportion_PerCent)) + geom_point(aes(color=Category)) + geom_line(aes(color=Category)) +
  ylab("Percentage of category in annotation (%)") +
  xlab("CADD score threshold") +
  scale_x_continuous(breaks=c(10,11,12,13,14,15,16,17,18,19,20),
                     labels=c("10" = expression("">="10"), "11" = expression("">="11"), "12" = expression("">="12"), "13" = expression("">="13"), "14" = expression("">="14"), "15" = expression("">="15"),
                              "16" = expression("">="16"),"17" = expression("">="17"),"18" = expression("">="18"),"19" = expression("">="19"),"20" = expression("">="20"))) +
  facet_wrap(~Annotation) +
  theme(strip.text = element_text(size=10), axis.text.x = element_text(size = 10, angle=75, hjust=1),
        axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.title=element_text(size=16), legend.text=element_text(size=16)) +
  scale_color_manual(values = c("firebrick4","darkblue"))
dev.off()

EOF
