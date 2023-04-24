############# HCV-induced HCC Gene expression profiling analysis script #############
#####################################################################################
####################################################################################

##########loading the required libraries###########

library(BiocManager)
library(hgu133plus2.db)
library(hgu133acdf)
library(limma)
library(hgu133plus2cdf)
library(hgu133a2.db)
library(GSEABase)
library(GOstats)
library(ggplot2)
library(curl)
library(RCurl)
library(affy)
library(readr)
library(hgu133a.db)
library(genefilter)
library(multtest)
library(affyPLM)
library(pheatmap)
library(affyio)
library(readr)
library(scatterplot3d)
library(rgl)
library(VennDiagram)
library(ggVennDiagram)
library(VennDiagram)
library(RDAVIDWebService)
library(GOplot)


#############################################################
##writing a function for multiple test correction:          #
correctPvalueandReturnAll<-function(tt.pval,method)         #      
{                                                           #  
  mt=mt.rawp2adjp(tt.pval,proc=method)                      #
  adjp=mt$adjp[order(mt$index),]                            #   
  return(adjp[,2])                                          #
}                                                           #   
######################################################################################################
####################################### Loading the files############################################
#####################################################################################################

#Loading the cel fileS list: 
Cel_files <- list.files("D:/Bioinformatics/NU courses/Integrative-Mohamed Hamed/project/working/GEO", pattern = "CEL", full.names = TRUE)
#Checking the chip types of the cel files:
table(sapply(Cel_files, function(x) read.celfile.header(x)$cdfName))

#Here, we have two batches, the first one is excluded due to its small number (9 samples).
#We will work only on the second batch which have 115 samples

#Reading the cel files to extract the files of the working batch( the second batch):
Batches<- split(Cel_files, sapply(Cel_files, function(x) read.celfile.header(x)$cdfName))
Batch_2 <- ReadAffy(filenames = Batches$"HG-U133A_2") # 732*732

######################################################################################################
################################# Data exploration of Batch_2##########################################
#######################################################################################################

class(Batch_2)
#Explore samples name 
sampleNames(Batch_2)
#Explore the affrmytrics chip probe (probe ids)
head(featureNames(Batch_2))
#Explore the anotation  
annotation(Batch_2)   # "hgu133a2"
#Explore the dimensions 
dim(Batch_2)
#Explore the expression of our raw data before any processing

#Visualization of raw data by using Histogram:
cols_batch2 <- seq(1:length(sampleNames(Batch_2)))  # 115
hist(Batch_2, main = "Figure 1 Histogram for raw data", col=cols_batch2, xlab= "log intensity") 
legend(15,1, sampleNames(Batch_2),col=cols_batch2,lty=1, lwd=1,cex=0.3)
boxplot(Batch_2, main = paste("Figure 2 Boxplot for raw data"), col = cols_batch2, names=paste("cel.",cols_batch2, sep=""),
        xlab="Samples" , ylab= "Sample Genes Intensity")

######################################################################################################
##############################Pre-processing of the data #############################################
#####################################################################################################

# three steps (background correction, normalization, summarization from probes to probesets)
# Normalization to remove florescence due to technical error, summarization here, remove miss match (MM) from perfect match (PM)related to probe
Batch2_set <- threestep(Batch_2,
                        background.method = "IdealMM",
                        normalize.method = "quantile",
                        summary.method = "median.polish")
hist(Batch2_set, main= "Figure 3 Histogram after Background Correction", xlab= "log intensity")
legend(12.5,0.25, sampleNames(Batch_2),col=cols_batch2,lty=1, lwd=1,cex=0.3)
boxplot(Batch2_set, main = paste("Figure 5 Boxplot after Background Correction"), col = cols_batch2, names=paste("cel.",cols_batch2, sep=""),
        xlab="Samples" , ylab= "Sample Genes Intensity")

######################################################################################################
#################################### Loading the processed data #####################################
#####################################################################################################

Batch_2_EXP_data= as.data.frame(exprs(Batch2_set))
Batch_2_EXP_data$probe_id <- row.names(Batch_2_EXP_data)
#mapping the probe ids into gene symbol
Batch_2_mapper <-as.data.frame(hgu133a2SYMBOL)
#Merging the two data frames (Batch_2_EXP_data and Batch_2_mapper) to add the gene symbols to our data
Merged_Batch_2 <- merge(Batch_2_EXP_data,Batch_2_mapper, by = "probe_id", all.y = T )
#Removing the probe_id column
Merged_Batch_2$probe_id <- NULL
#Check if we have any NAs or symbol duplication 
sum(is.na(Merged_Batch_2)) # no NA
sum(duplicated(Merged_Batch_2$symbol)) #7409 duplication
#Because we having duplication in symbol, we will do aggregation
#The symbol column will be removed, then convert all the data frame into numeric, so we can do aggregation
Merged_Batch_2_final<- Merged_Batch_2[-dim(Merged_Batch_2)[2]]
Merged_Batch_2_final <- apply(Merged_Batch_2_final, 2, as.numeric) #19808 * 115 before aggregation
#Aggregation
Aggregated_Batch_2 <- aggregate(Merged_Batch_2_final,list(Merged_Batch_2$symbol), FUN= mean)
row.names(Aggregated_Batch_2) <- Aggregated_Batch_2$Group.1
Aggregated_Batch_2$Group.1 <- NULL 
#Log transformation
#All_data_Exp = log2(Aggregated_Batch_2)   #12399 * 115 after aggregation ### I will not use log as the data is already log
All_data_Exp = Aggregated_Batch_2
#Upload sample info file
sample_info <-  read_delim("D:/Bioinformatics/NU courses/Integrative-Mohamed Hamed/project/working/E-GEOD-14323.sdrf.txt", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)


#we will rename the columns of our All_data_Exp with the scan names columns of the sample_info file
#however we need to match the ids of the two files 
All_data_Exp_names=names(All_data_Exp)  #115 names for Batch2 only
Exp_id=strsplit(All_data_Exp_names,".gz")
colnames(All_data_Exp)=c(Exp_id)
Exp_ids=colnames(All_data_Exp)
sample_info_ids= sample_info$`Scan Name`
index.files=match(Exp_ids,sample_info_ids)
colnames(All_data_Exp)=sample_info$`Comment [Sample_characteristics]`[index.files]  

#############################################################################################
############################## Differential expression analysis #############################
#############################################################################################

#Extract the target dataframes from the All_data_Exp which are:
Normal_df= All_data_Exp[,grep("Normal$", colnames(All_data_Exp))] ##### 19 samples
Cirrhosis_df = All_data_Exp[,grep("cirrhosis$", colnames(All_data_Exp))]   ##### 41 samples
HCC_df = All_data_Exp[,grep("Tissue: HCC", colnames(All_data_Exp))]   ##### 38 samples
HCV_HCC_df = All_data_Exp[,grep("cirrhosisHCC$", colnames(All_data_Exp))] ##### 17 samples

############################################################################################
#################################PCA Analysis###############################################
############################################################################################

#preparing step for PCA
# we need to get the gene expression matrix for the cases in just one matrix (we have dataframe)
all.con.data <- cbind(Normal_df, Cirrhosis_df, HCC_df, HCV_HCC_df) # 12399 *115

# we have 4 different conditions so we should have 4 different color
genotype <- colnames(all.con.data) 
condition <- c(replicate(19,"Normal"), replicate(41,"cirrhosis"),replicate(38,"HCC"), replicate(17,"Cirrhosis.HCC"))  ### 115
metadata.all.con.data <- data.frame(genotype, condition)

####2D PCA####

norm.all.con.data <- all.con.data # (we already normalized the data in previous steps)
pca <- prcomp(t(norm.all.con.data), scale= TRUE) # it has all values of pca
plot(pca$x[,1], pca$x[,2]) #it will not give fancy plotting for pca1 and pca2 only
pca.var <- pca$sdev^2
percentVar <- round(pca.var/sum(pca.var)*100,1)
barplot(percentVar, main= "Scree Plot", xlab=
          "Principal Component", ylab="Percent Variation")  # pca scree plot, we will have 114 PCAs

df_pca_data=data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], pc3 = pca$x[,3], sample = genotype, condition= condition)
pca2d<-ggplot(df_pca_data, aes(PC1,PC2, color = condition))+
  geom_point(size=3)+
  labs(x=paste0("PC1 (",percentVar[1],")"), y=paste0("PC2 (",percentVar[2],")"))
print(pca2d + ggtitle("               PCA in two dimension"))

####3D PCA###

scatterplot3d(df_pca_data[1:3], pch =20, color = "red")
plot3d(df_pca_data[1:3])
pca3d<-with(df_pca_data, plot3d(df_pca_data[1:3],type ="s",col = c("red","green","blue","yellow")))

################################################################################################
########DEGS Function To get All DEGs, Significant DEGs and its name  ##########################
#################################################################################################

get_DEGS_info <- function(all_cases, original_df){
  original_index = 1:ncol(original_df)
  All_DEGs_Res <- list()  
  All_sign_DEGs <- list()     #####sign of HCC&Normal = 6214, sign of HCC&HCV-HCC = 3225
  All_sign_DEGs_genes <- list()
  for (i in 1:length(all_cases)){
    ## calculating LFC (log fold change) 
    lfc_diff=apply(all_cases[[i]],1, function(x)  mean(x[original_index]) -mean(x[(ncol(original_df)+1) : ncol(all_cases[[i]])]))
    lfc_diff_df=as.data.frame(lfc_diff)
    ## calcualting p values: 
    f=factor( c( rep(1, length(original_index)) , rep(2, length((ncol(original_df)+1) : ncol(all_cases[[i]])) )))
    t_pval=rowttests(as.matrix(all_cases[[i]][,1:ncol(all_cases[[i]])]),f)$p.value
    t_pval_adj=correctPvalueandReturnAll(t_pval,"BH")
    ## res= results
    All_DEGs_Res[[i]] <- cbind(lfc_diff_df, t_pval, t_pval_adj)
    
    All_sign_DEGs[[i]] <-All_DEGs_Res[[i]][abs(lfc_diff_df) > log2(2)   & t_pval_adj <0.05,]  ##### identify DEGs based on both  LFC and the significance level 
    #log2(2)=1
    All_sign_DEGs_genes[[i]] <- rownames(All_sign_DEGs[[i]])
    
  }
  return(list(Res= All_DEGs_Res, DEGs_sign=All_sign_DEGs, sign_genes =  All_sign_DEGs_genes))
}


##############################################################################################
##############################################################################################


##### DEGs estimating for comparing (Normal & Cirrhosis) and (Normal & HCV-HCC): 
Normal_Cirrhosis <- cbind(Normal_df, Cirrhosis_df)   ### 60 samples
Normal_HCV_HCC <- cbind(Normal_df, HCV_HCC_df)       ### 36 samples

All_Normal_cases <- list(Normal_Cirrhosis, Normal_HCV_HCC)

Normal_DEGS_info <- get_DEGS_info(All_Normal_cases, Normal_df)
All_Normal_DEGs_Res <- Normal_DEGS_info[['Res']]  ### 1st d.f.>>>>Normal_Cirrhosis=12399, 2nd d.f>>>>Normal_HCV_HCC=12399   .... constant for different LFC 
All_Normal_Sign_DEGs <- Normal_DEGS_info[['DEGs_sign']] 
#for LFC=log2(1)=0    ### 1st d.f.>>>>Normal_Cirrhosis=6520 ,2nd d.f.>>>>>>Normal_HCV_HCC=5084
#for LFC= log2(2)= 1  ### 1st d.f.>>>>Normal_Cirrhosis=1286 ,2nd d.f.>>>>>>Normal_HCV_HCC=1047
All_Normal_Sign_DEGs_genes <- Normal_DEGS_info[['sign_genes']]
#for LFC=log2(1)=0    ### 1st d.f.>>>>Normal_Cirrhosis=6520 ,2nd d.f.>>>>>>Normal_HCV_HCC=5084
#for LFC= log2(2)= 1  ### 1st d.f.>>>>Normal_Cirrhosis=1286 ,2nd d.f.>>>>>>Normal_HCV_HCC=1047

##### DEGs estimating for comparing HCC tissue with other types tissues,
#####(HCC & Normal, (HCC &  Cirrhosis) and (HCC & HCV-HCC) : 
HCC_Cases_list <- list(Normal_df, Cirrhosis_df, HCV_HCC_df)
Final_HCC_Cases_list <- lapply(HCC_Cases_list, function(x) cbind(HCC_df,x) ) #HCC_df 38 samples
HCC_Normal <- Final_HCC_Cases_list[[1]]        ### 57 samples
HCC_Cirrhosis <- Final_HCC_Cases_list[[2]]     ### 79 samples
HCC_HCV_HCC <- Final_HCC_Cases_list[[3]]       ### 55 samples

All_HCC_cases <- list(HCC_Normal,HCC_Cirrhosis,HCC_HCV_HCC)

HCC_DEGS_info <- get_DEGS_info(All_HCC_cases, HCC_df)
All_HCC_DEGs_Res <- HCC_DEGS_info[['Res']]
### HCC_Normal=12399, HCC_Cirrhosis=12399, HCC_HCV_HCC= 12399 .... constant for different LFC
All_HCC_Sign_DEGs <- HCC_DEGS_info[['DEGs_sign']]
#for LFC=log2(1)    ### HCC_Normal=5400, HCC_Cirrhosis=5011, HCC_HCV_HCC= 1308 
#for LFC=log2(2)=1  ### HCC_Normal=853, HCC_Cirrhosis=308,   HCC_HCV_HCC= 105 
All_HCC_Sign_DEGs_genes <- HCC_DEGS_info[['sign_genes']]     
#for LFC=log2(1)    ### HCC_Normal=5400, HCC_Cirrhosis=5011, HCC_HCV_HCC= 1308 
#for LFC=log2(2)=1  ### HCC_Normal=853, HCC_Cirrhosis=308,   HCC_HCV_HCC= 105

##### DEGs estimating for comparing (Cirrhosis & HCV-HCC):
Cirrhosis_HCV_HCC <- cbind(Cirrhosis_df, HCV_HCC_df)  ### 58 samples
Cirrhosis_case <- list(Cirrhosis_HCV_HCC)
Cirrhosis_DEGS_info <- get_DEGS_info(Cirrhosis_case, Cirrhosis_df)
All_Cirrhosis_DEGs_Res <- Cirrhosis_DEGS_info[['Res']] ### Cirrhosis_HCV_HCC = 12399 ...constant for different LFC
All_Cirrhosis_Sign_DEGs <- Cirrhosis_DEGS_info[['DEGs_sign']]  
#for LFC=log2(1)   ### Cirrhosis_HCV_HCC = 1314
#for LFC=log2(2)=1 ### Cirrhosis_HCV_HCC = 47
All_Cirrhosis_Sign_DEGs_genes <- Cirrhosis_DEGS_info[['sign_genes']]   
#for LFC=1  ### Cirrhosis_HCV_HCC = 1314
#for LFC=log2(2)=1  ### Cirrhosis_HCV_HCC = 47


###################################################################
##Our target study will be only on (HCC&Normal) and (HCC& HCV-HCC)#
#for LFC=log2(2)=1 ### HCC_Normal=853,  HCC_HCV_HCC= 47           #    
###################################################################


##Extracting the unique DEGs that are involved in HCV in the case of HCV-HCC:
rows_unique <- !(rownames(All_HCC_Sign_DEGs[[3]]) %in%  rownames(All_HCC_Sign_DEGs[[1]]))
HCV_HCC_DEGs_unique <- All_HCC_Sign_DEGs[[3]][rows_unique,] 
#for LFC=log2(2)=1 ##### 69
degs=rownames(HCV_HCC_DEGs_unique)
## Saving the unique DEGs in TXT file
write.table(degs, 'D:/Bioinformatics/NU courses/Integrative-Mohamed Hamed/project/HCV_HCC_DEGs_names_unique.txt', col.names=F, row.names=F, quote = FALSE)
write.table(HCV_HCC_DEGs_unique, 'D:/Bioinformatics/NU courses/Integrative-Mohamed Hamed/project/HCV_HCC_DEGs_unique.txt', quote = FALSE)
HCV_HCC_DEGs_unique_df=as.data.frame(HCV_HCC_DEGs_unique)

##########################################################################################
################ Make a basic volcano plot#################################################
#########################################################################################

par("mar")
par(mfrow=c(1,1))
with(All_Cirrhosis_DEGs_Res[[1]], plot(lfc_diff, -log10(t_pval_adj), pch=20, main="Volcano plot for HCV with HCV_HCC_df", xlim=c(-3,3)))
# Add colored points: blue if padj<0.05, red if log2FC>2 and padj<0.05)
with(subset(All_Cirrhosis_DEGs_Res[[1]], All_Cirrhosis_DEGs_Res[[1]]$t_pval_adj<.05 ), points(lfc_diff, -log10(t_pval_adj), pch=20, col="blue"))
with(subset(All_Cirrhosis_DEGs_Res[[1]], t_pval_adj<.05 & lfc_diff> log2(2)), points(lfc_diff, -log10(t_pval_adj), pch=20, col="red"))
with(subset(All_Cirrhosis_DEGs_Res[[1]], t_pval_adj<.05 & lfc_diff < -log2(2)), points(lfc_diff, -log10(t_pval_adj), pch=20, col="green"))
################
with(All_Normal_DEGs_Res[[1]], plot(lfc_diff, -log10(t_pval_adj), pch=20, main="Volcano plot for HCC with Normal_df", xlim=c(-5,5)))
with(subset(All_Normal_DEGs_Res[[1]], t_pval_adj<.05 ), points(lfc_diff, -log10(t_pval_adj), pch=20, col="blue"))
with(subset(All_Normal_DEGs_Res[[1]], t_pval_adj<.05 & lfc_diff >log2(2)), points(lfc_diff, -log10(t_pval_adj), pch=20, col="red"))
with(subset(All_Normal_DEGs_Res[[1]], t_pval_adj<.05 & lfc_diff < -log2(2)), points(lfc_diff, -log10(t_pval_adj), pch=20, col="green"))
################
with(All_HCC_DEGs_Res[[2]], plot(lfc_diff, -log10(t_pval_adj), pch=20, main="Volcano plot for HCC with HCV_HCC_df", xlim=c(-3,3)))
with(subset(All_HCC_DEGs_Res[[2]], t_pval_adj<.05 ), points(lfc_diff, -log10(t_pval_adj), pch=20, col="blue"))
with(subset(All_HCC_DEGs_Res[[2]], t_pval_adj<.05 & lfc_diff > log2(2)), points(lfc_diff, -log10(t_pval_adj), pch=20, col="red"))
with(subset(All_HCC_DEGs_Res[[2]], t_pval_adj<.05 & lfc_diff < -log2(2)), points(lfc_diff, -log10(t_pval_adj), pch=20, col="green"))

########################################################################################
###############  ####Venn Diagram  ####  #### ##########################################
###############################################################################

##Visualization by venn diagram between HCC-Normal & HCC-HCV_HCC genes:
#using LFC=1

venn_genes1 <- list('A'= All_HCC_Sign_DEGs_genes[[1]], 'B'=All_Normal_Sign_DEGs_genes[[2]],
                    'C'=All_Normal_Sign_DEGs_genes[[1]])

#A: Normal& HCC
#B: Normal& HCV_HCC
#C: Normal& HCV
ggVennDiagram(venn_genes1, color = "grey")

venn_genes2 <- list('A'= All_HCC_Sign_DEGs_genes[[1]], 'B'=All_HCC_Sign_DEGs_genes[[2]], 
                    'C'= All_HCC_Sign_DEGs_genes[[3]])
#A: HCC& Normal 
#B: HCC& HCV 
#C: HCC& HCC-HCV
ggVennDiagram(venn_genes2, color = "grey")

venn_genes3 <- list('A'=All_Normal_Sign_DEGs_genes[[1]], 'B'=All_HCC_Sign_DEGs_genes[[2]], 
                  'C'= All_Cirrhosis_Sign_DEGs_genes[[1]])
#A: HCV& Normal
#B: HCV& HCC
#C: HCV& HCC_HCV
ggVennDiagram(venn_genes3, color = "grey")

venn_genes4 <- list('A'=All_Normal_Sign_DEGs_genes[[2]],  'B'= All_HCC_Sign_DEGs_genes[[3]], 
                    'C'= All_Cirrhosis_Sign_DEGs_genes[[1]])
#A: HCC-HCV& Normal
#B: HCC-HCV& HCC 
#C: HCC_HCV& HCV
ggVennDiagram(venn_genes4, color = "grey")

##########################################################################################
###########################  heat map  ######################################################3
##############################################################################################

exp.degs=HCC_HCV_HCC[rownames(HCC_HCV_HCC) %in% rownames(All_HCC_Sign_DEGs[[3]]),]
dsm=exp.degs

HCC_index =c(1:38)
HCC_HCV_index = c(39:dim(HCC_HCV_HCC)[2])
index1 <- c()
index2 <- c()
for (i in 1:length(HCC_index)) {
  index1[i] <- paste("HCC", sep = ".",i)
}
for (i in 1:length(HCC_HCV_index)) {
  index2[i] <- paste("HCV_HCC", sep = ".",i)
}


Colnames <- c(index1,index2)
colnames(dsm) <- Colnames

HCC.vector=rep("HCC", length(HCC_index))
HCV_HCC.vector=rep("HCV_HCC", length(HCC_HCV_index ))
Sample=c(HCC.vector, HCV_HCC.vector) 
annotation=as.data.frame(Sample)
rownames(annotation)=Colnames

# Specify colors
Sample = c("lightgreen", "navy")
names(Sample) = c("HCC", "HCV_HCC")
ann_colors = list(Sample = Sample)
m2=scale(t(dsm),center=T,scale=T)
m2=t(m2)
pheatmap(m2, annotation = annotation, annotation_colors = ann_colors, fontsize_row = 3, fontsize_col = 6, main = "Heat map for HCC& HCV_HCC")

#########################################################################################
###################### Enrichment analysis& GO Plot #####################################
########################################################################################

#preparing the output of david to be used in Go plot
david.output <-  read_delim("D:/Bioinformatics/NU courses/Integrative-Mohamed Hamed/project/chart_F962E8C77D121613507973691.txt",
                            "\t", escape_double = FALSE, trim_ws = TRUE)
sub.david.output <- david.output[c(1, 2, 6, 12 )]


first.col <- c()
second.col <- c()
third.col <- c()
for (i in 1:dim(sub.david.output)[1]) {
  first.col[i] <- strsplit(sub.david.output$Category[i], "_")[[1]][2]
  second.col[i] <-  strsplit(sub.david.output$Term[i],"~")[[1]][1]
  third.col[i] <- strsplit(sub.david.output$Term[i], "~")[[1]][2]
}


go.plot.input<- data.frame(first.col, second.col, third.col, 
                           sub.david.output$Genes,
                           sub.david.output$Benjamini)
colnames(go.plot.input) <- c("Category", "ID", "Term", "Genes", "adj_pval")
write.table(go.plot.input, "D:/Bioinformatics/NU courses/Integrative-Mohamed Hamed/project/go.plot.input.txt", quote = FALSE,
            sep = '\t',row.names = F, col.names = T)


#for uniq gene preparation
View(HCV_HCC_DEGs_unique)
uniq.genes <- data.frame(rownames(HCV_HCC_DEGs_unique), HCV_HCC_DEGs_unique$lfc_diff)
colnames(uniq.genes) <- c("ID", "logFC")
write.table(uniq.genes, "D:/Bioinformatics/NU courses/Integrative-Mohamed Hamed/project/uniq.genes.txt", quote = FALSE,
            sep = '\t',row.names = F, col.names = T)


david= read.delim("D:/Bioinformatics/NU courses/Integrative-Mohamed Hamed/project/go.plot.input.txt",header = T, sep = "\t")
head(david)

genelist = read.delim("D:/Bioinformatics/NU courses/Integrative-Mohamed Hamed/project/uniq.genes.txt",header = T, sep = "\t")
head(genelist)

# Generate the plotting object
circ <- circle_dat(david,genelist)
sub.cir <- circ[ grep("metabolic", circ$term),]

# Generate a simple barplot
GOBar(subset(circ, category == 'BP'))
# Facet the barplot according to the categories of the terms 
GOBar(circ, display = 'multiple')
# Facet the barplot, add a title and change the colour scale for the z-score
GOBar(circ, display = 'multiple', title = 'Z-score coloured barplot', order.by.zscore = T,zsc.col = c('yellow', 'black', 'cyan'))

# Generate the bubble plot with a label threshold of zero
GOBubble(circ,title = 'Bubble plot',bg.col = F, colour = c('orange', 'darkred', 'gold'),ID = T, table.legend= T,display = 'multiple',table.col=F,labels = 0)
# Reduce redundant terms with a gene overlap >= 0.75.and plot on adjusted p-value =1
reduced_circ <- reduce_overlap(circ, overlap = 0.75)
#sub.reduced_circ <-reduced_circ[ grep("metabolic", reduced_circ$term),]
sub.reduced_circ <- reduced_circ[reduced_circ$adj_pval==1,]
GOBubble(reduced_circ, labels = 0)
#GOBubble(sub.reduced_circ, labels = 0)
#Generate a circular visualization of the results of gene- annotation enrichment analysis
GOCircle(circ)
# Generate a circular visualization for 12 terms
GOCircle(circ, nsub = 12)

# Generate the matrix with a list of selected genes and the matrix with selected processes .
#chord <- chord_dat(data = circ, genes = genelist, process = unique(sub.cir$term))
chord <- chord_dat(data = circ, genes = genelist, process = sub.reduced_circ$term)
GOChord(chord, limit = c(0,1), gene.order = 'logFC')
GOChord(chord, space=0.02, gene.order='logFC',lfc.col=c('red','black','cyan'))

########################################################################################################
##################################### Heatmap of genes and terms (GOHeat)#################################
#########################################################################################################
## First, we use the chord object without logFC column to create the heatmap
#GOHeat(chord[,-8], nlfc = 1)
# Now we create the heatmap with logFC values and user-defined colour scale

GOHeat(chord, nlfc = 1, fill.col = c('red', 'yellow', 'green'))


#######################################################################################################3
###################################TFmir################################################################3
######################################################################################################

#TFmir, convert the previus txt file to  another format to be read by Tfmir
TFmirdata <- HCV_HCC_DEGs_unique
TFmirdata$flag <- HCV_HCC_DEGs_unique$lfc_diff

TFmirdata[TFmirdata$lfc_diff>0,]$flag = 1  # to check if the value of LFC>0
TFmirdata[TFmirdata$lfc_diff<0,]$flag= -1
table(TFmirdata$flag)  # I have 657 for -1 flag and I have 806 for 1 flag
TFmir<- data.frame(rownames(TFmirdata), TFmirdata$flag)

write.table(TFmir, file = "D:/Bioinformatics/NU courses/Integrative-Mohamed Hamed/project/TFmir.common.DEGs.txt", sep="\t",quote = F , row.names= F, col.names = F )

###################################################################################################
########################################### CMAP ##################################################
##################################################################################################

## Extracting the up regulated and down regulated genes to be inserted in the L1000CDS2 tool
# to estimating the disease signature and drug signature
## For up regulated genes:
upregulated_genes_index <- HCV_HCC_DEGs_unique$lfc_diff > 0
All_upregulated_genes <- rownames(HCV_HCC_DEGs_unique[upregulated_genes_index,]) ###334 genes
write.table(All_upregulated_genes, 'D:/Bioinformatics/NU courses/Integrative-Mohamed Hamed/project/All_upregulated_genes.txt', quote = FALSE,
            sep = '\t', row.names = FALSE)

## For down regulated genes:
downregulated_genes_index <- HCV_HCC_DEGs_unique$lfc_diff < 0
All_downregulated_genes <- rownames(HCV_HCC_DEGs_unique[downregulated_genes_index,]) ###305 genes
write.table(All_downregulated_genes, 'D:/Bioinformatics/NU courses/Integrative-Mohamed Hamed/project/All_downregulated_genes.txt', quote = FALSE,
            sep = '\t', row.names = FALSE)



