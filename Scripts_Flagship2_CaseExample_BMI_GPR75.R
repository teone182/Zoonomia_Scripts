library(dplyr)
library(ggplot2)

#####Figure S5_2A####

#####Data_ForFigure_S5_2A####
#Class	Feature	FracPhyloP
#"CDS"	"GPR75 pLOF"	0.70
#"Exon 1"	"GPR75 pLOF	1
#"Exon 2"	"GPR75 pLOF"	0.69
#"CDS"	"GPR75 all bps"	0.587
#"Exon 1"	"GPR75 all bps"	0.435
#"Exon 2"	"GPR75 all bps"	0.506
################################

data <- read.csv("Data_ForFigure_S5_2A", sep="\t", head=T)
pdf("Figure_S5_2A.pdf")
ggplot(data=data, aes(x=Feature, y=FracPhyloP, fill=Class)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values = c("#736F6E", "#98AFC7","#153E7E")) +
    ylab("Fraction Constraint (5% FDR PhyloP)") +
  xlab("") +
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=12),
      axis.title=element_text(size=14),
      legend.title=element_text(size=12), legend.text=element_text(size=12))
dev.off()



#####Figure S5_2B####

#########Data_ForFigure_S5_2B####
#Name	PhyloP	Group	Position	Class
#"-110+1G>A"	7.901000	1	1	"SNV"
#"Met1?"	6.134000	1	2	"SNV"
#"His6fs"	3.805000	1	3	"Non-SNV"
#"Ser21*"	-0.057000	1	4	"SNV"
#"Gly24fs"	0.502000	1	5	"Non-SNV"
#"His38fs"	6.766000	1	6	"Non-SNV"
#"Leu91fs"	 8.598000	1	7	"Non-SNV"
#"Cys94fs"	6.195000	1	8	"Non-SNV"
#"Gly95*"	5.794000	1	9	"SNV"
#"Ala110fs"	0.5591429	1	10	"Non-SNV"
#"Cys118fs"	3.253000	1	11	"Non-SNV"
#"Ser126*"	8.598000	1	12	"SNV"
#"Gln151*"	6.772000	1	13	"SNV"
#"Leu206fs"	4.817	1	14	"Non-SNV"
#"Tyr207fs"	1.6145	1	15	"Non-SNV"
#"Ser219fs"	 -0.27	1	16	"Non-SNV"
#"Gln227*"	6.817000	1	17	"SNV"
#"Gln234*"	6.819000	1	18	"SNV"
#"Arg236*"	-0.346000	1	19	"SNV"
#"Arg236fs"	1.89875	1	20	"Non-SNV"
#"Val241fs"	5.697000	1	21	"Non-SNV"	
#"Gln250*"	5.610000	1	22	"SNV"
#"Pro263fs"	2.6265	1	23	"Non-SNV"
#"Leu271fs"	"NA"	"NA"	24	"NA"
#"Tyr277fs"	4.483000	1	25	"Non-SNV"
#"Tyr277*"	3.815000	1	26	"SNV"
#"Gln294*"	5.543000	1	27	"SNV"
#"Arg302*"	3.088000	1	28	"SNV"
#"Ser329*"	8.598000	1	29	"SNV"
#"Gln343*"	8.598000	1	30	"SNV"
#"Tyr355*"	2.766000	1	31	"SNV"
#"Tyr355fs"	4.118385	1	32	"Non-SNV"
#"Asn372fs"	5.4025	1	33	"Non-SNV"
#"Cys400fs"	1.280000	1	34	"Non-SNV"
#"Lys404*"	6.220000	1	35	"SNV"
#"Arg408*"	4.687000	1	36	"SNV"
#"Arg419*"	3.609000	1	37	"SNV"
#"Lys459fs"	 2.131000	1	38	"Non-SNV"
#"Cys467fs"	1.6095	1	39	"Non-SNV"
#"Gln469*"	6.737000	1	40	"SNV"
#"Gln501*"	6.793000	1	41	"SNV"
#"Asn504fs"	-0.038000	1	42	"Non-SNV"
#"Thr519fs"	4.517	1	43	"Non-SNV"
#"Ter541Gluext*?"	4.665000	1	44	"SNV"
#"Ter541Serext*?"	0.757000	1	45	"SNV"
#"Ter541Tyrext*?"	1.837000	1	46	"SNV"
#############################################


data <- read.csv("Data_ForFigure_S5_2B", head=T, sep="\t")
data$Name <- factor(data$Name, levels = data$Name[order(data$Position)])
fit0 <- lm(data$PhyloP ~ data$Position)
summary(fit0)
fit1_data <- subset(data, data$Class == "SNV")
fit1 <- lm(fit1_data$PhyloP ~ fit1_data$Position)
summary(fit1)
pdf("Figure_S5_2B.pdf")
ggplot(data, aes(x=Name, y=PhyloP)) + geom_point(aes(group = 1, color=Class)) + 
  scale_color_manual(values = c("SNV" = "black", "Non-SNV" = "red")) +
  geom_line(aes(group = 1),linetype = "dashed") +
  theme(text = element_text(size=13),
        axis.text.x = element_text(angle=90, hjust=1)) +
  geom_abline(intercept = coef(fit0)[1], slope = coef(fit0)[2], col='blue') +
  geom_abline(intercept = coef(fit1)[1], slope = coef(fit1)[2], col='black') +
  ylab("PhyloP score") +
  xlab("GPR75 pLOF variants") +
  theme(legend.position = "none") +
  coord_fixed(ratio=2.5)
dev.off()




#####Figure S5_2C####

#########Data_ForFigure_S5_2B####
#Beta	PhyloP	L	U	Group
#-2.2	4.50866	-2.9	-1.4	"N-terminal variants: truncation before last 100 amino acids"
#-2.1	4.552041	-2.9	-1.3	"N-terminal variants: truncation before last intracellular domain"
#-1.4	3.586214	-2.5	-0.3	"C-terminal variants: truncation within last intracellular domain"
#-0.7	3.223167	-2.1	0.7	"C-terminal variants: truncation within last 100 amino acids"
####################################


data <- read.csv("Data_ForFigure_S5_2B", head=T, sep="\t")
pdf("Figure_S5_2C.pdf")
plot(data$Beta, data$PhyloP, xlab = "Beta (95% CI) per allele in kg/m2 units of BMI", ylab= "Mean PhyloP", cex=2, pch=20, col=seq_along(data$Group), xlim = c(-3,1))  
arrows(x0=data$L, y0=data$PhyloP, x1=data$U, y1=data$PhyloP, code=3, angle=90, length=0.05)
abline(lm(data$PhyloP~data$Beta), col="black")
legend("topright", legend=data$Group, cex=0.6, pch=20, col=seq_along(data$Group))
dev.off()


