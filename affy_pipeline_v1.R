##analysis Human gene 2.1 ST
require(oligo)
require(limma)
require(pd.hugene.2.1.st)
require(pd.mogene.2.1.st)
require(genefilter)
require(ggplot2)
require(reshape)
require(dendextend)
require(colorspace)
require(scatterplot3d)

########## General arguments ###########
# this script should be located in the root folder of the analysis
# cel files should be in a cels directory

#no group?
noGroup=F

#select array type:
#arrayType="human"
arrayType="mouse"

#filtering parameter:ra
global_nsamples=4
group_nsamples=3

#read info and cel files
cat("Reading info file\n")
info<-read.table("sample_info2.txt", sep="\t",header=T) #should have columns (Sample.Name, Sample.Group)
cat("Info file loaded\n")


########################################
###########loading the data#############
cat("Loading data\n")
cels.files<-list.celfiles("cels/",full.names=T)
raw.data<-read.celfiles(cels.files)
cat("Data loaded\n")
#########################################


##########Process sample names################

sampleNames(raw.data)<- gsub(sampleNames(raw.data), pattern = ".CEL",replacement = "")
index_sample<-match(info$Sample.Name,sampleNames(raw.data) )


##########Selected sample present in the info file ######

raw.data<-raw.data[,index_sample]


##########Get group data ##############################

group<-levels(as.factor(info$Sample.Group))
ngroup<-length(levels(as.factor(info$Sample.Group)))


############create all the directory needed######

cat("Creating directories\n")
dir.create("controls", showWarnings = F)
dir.create("raw_data_plots", showWarnings = F)
dir.create("rma_plots", showWarnings = F)
dir.create("results/graphs", showWarnings = F,recursive = T)
dir.create("results/tables", showWarnings = F,recursive = T)


############plot boxplots#######################

cat("Starting plotting raw data\n")
cat("Ploting boxplot\n")
pdf(file="raw_data_plots/boxplot_raw_groups.pdf",width=10,height=5)
boxplot(raw.data, col=info$Sample.Group,las=2)
dev.off()


#############Plot MDS clustering (group coloring)########

cat("Ploting HC\n")
pdf(file="raw_data_plots/Hierarchical_clustering_raw_data.pdf",width=10, height=5)
plot(hclust(dist(t(exprs(raw.data)))),xlab="Distance based on expression",main="Hierarchical clustering raw data")
dev.off()

cat("Ploting MDS\n")
pdf(file="raw_data_plots/MDS_clustering_raw_data.pdf",width=10, height=5)
fit_mds_raw<-cmdscale(dist(t(exprs(raw.data))),eig = T, k=3)
pca1_raw<-fit_mds_raw$points[,1]
pca2_raw<-fit_mds_raw$points[,2]
par(mar=c(5.1,4.1,4.1,12.1),xpd=T)
plot(pca1_raw,pca2_raw,type="n")
if(!noGroup){
  text(pca1_raw,pca2_raw,labels=sampleNames(raw.data),cex=0.7,col=as.numeric(info$Sample.Group))
  legend(legend=levels(info$Sample.Group),fill=c(1:length(levels(info$Sample.Group))), "topright",inset=c(-0.2,0))
}
dev.off()


############### 3D plots ###########

cat("Plotting 3D MDS\n")
pdf(file="raw_data_plots/MDS_clustering_raw_data_3D.pdf",width=10, height=5)
scatterplot3d(fit_mds_raw$points[,c(1:3)], xlab="PCA1",ylab="PCA2", zlab="PCA3", color=as.numeric(info$Sample.Group),pch=19, type="h" )
#alternative (interactive): plot3d(fit_mds_raw$points[,c(1:3)], col=as.numeric(info$Sample.Group), siz=1, type="s", xlab="PCA1", ylab="PCA2", zlab="PCA3")
#scatterplot3d(pca$rotation[,c(1:3)], xlab="PCA1",ylab="PCA2", zlab="PCA3", color=as.numeric(info$Sample.Group),pch=19, type="h" )
legend("topleft", title="Groups", levels(as.factor(info$Sample.Group)), fill=levels(as.factor(as.numeric(info$Sample.Group))))
dev.off()

############# coloured dendogram #######

cat("Generating coloured dendogram\n")
hc<-hclust(dist(t(exprs(raw.data))))
pdf(file="raw_data_plots/Hierarchical_clustering_raw_data_coloured.pdf",width=10, height=5)
dend<-as.dendrogram(hc)
dend <- rotate(dend, 1:dim(info)[1])
dend <- color_branches(dend, k=ngroup)
dend <- set(dend, "labels_cex", 0.7)
match_order<-match(labels(dend), info$Sample.Name)
group_order<-info$Sample.Group[match_order]
label_legend<-unique(group_order)
labels_colors(dend) <-  rainbow_hcl(ngroup)[sort_levels_values(as.numeric(info$Sample.Group)[order.dendrogram(dend)])]
plot(dend, main="Hierarchical cluster of sample groups", horiz=T, nodePar=list(cex=.007))
legend("topleft", title="Groups", as.character(label_legend), fill=rainbow_hcl(ngroup))
par(mar=c(5.1,4.1,4.1,2.1),xpd=F)
dev.off()


########### plot GC content intensity dependance#########

cat("Plotting GC content intensity dependence\n")
pmSeq<-pmSequence(raw.data)
pmsLog2<-log2(pm(raw.data))
coefs<-getAffinitySplineCoefficients(pmsLog2, pmSeq)
counts<-Biostrings::alphabetFrequency(pmSeq,baseOnly=T)
GCcontent<-ordered(counts[,"G"]+counts[,"C"])
total_graph=ceiling(length(sampleNames(raw.data))/4)
graph_n=0
for(i in 1:length(sampleNames(raw.data))){
  current_n=ceiling(i/4)
  if(current_n==graph_n)
  {
    plot(y=pmsLog2[,i], x=GCcontent, main=paste("GC content ",sampleNames(raw.data)[i]), xlab="GC frequency", ylab="Intensity" )
  }else{
    if(current_n>1){
      dev.off()
    }
    graph_n=current_n
    jpeg(file=paste("raw_data_plots/GCContent_",graph_n,"_raw_data.jpg",sep=""),width = 480, height=480)
    par(mfrow=c(2,2))
    plot(y=pmsLog2[,i], x=GCcontent, main=paste("GC content ",sampleNames(raw.data)[i]), xlab="GC frequency", ylab="Intensity" )
  } 
}
dev.off()


########### Plot density plot ####################

cat("Plot density plot\n")
pdf(file="raw_data_plots/density_plot_raw_data.pdf",width=10, height=5)
par(mar=c(5.1,4.1,4.1,12.1),xpd=T)
hist(raw.data, col=1:length(levels(info$Sample.Group)))
legend(legend=levels(info$Sample.Group),fill=c(1:length(levels(info$Sample.Group))), "topright",inset=c(-0.2,0))
dev.off()
par(mar=c(5.1,4.1,4.1,2.1),xpd=F)

########## Plot samples MA plot vs agregate from the others sample ########
cat("Plotting MA data\n")
graph_n=0
for(i in 1:length(sampleNames(raw.data))){
  current_n=ceiling(i/4)
  if(current_n==graph_n)
  {
    MAplot(raw.data, which=i)
  }else{
    if(current_n>1){
      dev.off()
    }
    graph_n=current_n
    jpeg(file=paste("raw_data_plots/MAplot_",graph_n,"_raw_data.jpg",sep=""),width = 800, height=500)
    par(mfrow=c(2,2))
    MAplot(raw.data, which=i)
  } 
}
dev.off()

########## plot group HC #######
cat("Plot group HC and MDS\n")
if(!noGroup){
  for(i in 1:length(group)){
    selectedSamples<-info[info$Sample.Group==as.character(group[i]),]$Sample.Name
    pdf(file=paste("raw_data_plots/Hist_",group[i],"_raw_data.pdf",sep=""),width = 10, height=5)
    hist(raw.data[,selectedSamples],main=group[i],col=1:length(selectedSamples), lty=1:length(selectedSamples))
    legend("topright", legend=selectedSamples, col=1:length(selectedSamples),lty=1:length(selectedSamples))
    dev.off()
    
    pdf(file=paste("raw_data_plots/MDS_",group[i],"_raw_data.pdf",sep=""),width = 10, height=5)
    fit_mds_raw<-cmdscale(dist(t(exprs(raw.data[,selectedSamples]))),eig = T, k=2)
    pca1_raw<-fit_mds_raw$points[,1]
    pca2_raw<-fit_mds_raw$points[,2]
    par(mar=c(5.1,4.1,4.1,12.1),xpd=T)
    plot(pca1_raw,pca2_raw,type="n")
    dev.off()
  }
}




############# normalisation rma data ########

cat("Performing RMA normalisation\n")
rma.data.raw<-rma(raw.data, background=F, normalize=F) #to get the controls
rma.probe<-rma(raw.data, target="probe")
rma.data.core<-rma(raw.data, target="core")
rma.data.core<-rma.data.core[,as.character(info$Sample.Name)]
median_exprs<-median(exprs(rma.data.core))


############# Filter expression levels #########
cat("Filtering dataa\n")
filter_setting<-kOverA(global_nsamples,median_exprs)
filter_function<-filterfun(filter_setting)
probe_filter<-genefilter(exprs(rma.data.core),filter_function)
rma.data.core<-rma.data.core[probe_filter,]
data.annot<-getNetAffx(rma.data.core, type="transcript")
fData(rma.data.core)<-pData(data.annot)
rma.data.core<-rma.data.core[fData(rma.data.core)$category=="main",]

############# Filter for each group and generate MDS #######
cat("Plotting MDS post normalisation for each group\n")
for(i in 1:ngroup){
  current_group<-group[i]
  cat(paste0(current_group,"\n"))
  current_samples<-info[info$Sample.Group==current_group,]
  color_selected<-rainbow_hcl(length(current_samples$Sample.Name))
  selected_rma_data<-rma.data.core[,match(as.vector(current_samples[,1]),sampleNames(rma.data.core))]
  if(dim(selected_rma_data)[2]>3){
    pdf(file=paste("rma_plots/MDS_clustering_rma_data_3D_",current_group,".pdf",sep=""),width=10, height=5)
    par(xpd = T)
    fit_mds_rma<-cmdscale(dist(t(exprs(selected_rma_data))),eig = T, k=3)
    write.table(fit_mds_rma$points,file=paste("rma_plots/MDS_clustering_rma_data_table_",current_group,".txt",sep=""), sep="\t" )
    scatterplot3d(fit_mds_rma$points[,c(1:3)], xlab="MDS1",ylab="MDS2", zlab="MDS3",pch=19, type="h",color=color_selected[1:length(current_samples$Sample.Name)])
    legend("topleft", title="Samples", legend = current_samples$Sample.Name, fill=color_selected, cex = 0.7, inset=c(-0.05,-0.1))
    dev.off()
  }

  pdf(file=paste("rma_plots/MDS_clustering_rma_data_",current_group,".pdf",sep=""),width=10, height=5)
  par(xpd = T)
  fit_mds_rma<-cmdscale(dist(t(exprs(selected_rma_data))),eig = T, k=2)
  mds1<-fit_mds_rma$points[,1]
  mds2<-fit_mds_rma$points[,2]
  plot(mds1, mds2, main="MDS plot", xlab="MDS 1", ylab="MDS 2", col=rainbow_hcl(length(current_samples$Sample.Name)), pch=16)
  legend("topleft", title="Samples", legend = current_samples$Sample.Name, fill=rainbow_hcl(length(current_samples$Sample.Name)), cex=0.7, inset(0,-0.1))
  dev.off()
  pdf(file=paste("rma_plots/Hierarchical_clustering_rma_data_",current_group,".pdf",sep=""),width=10, height=5)
  plot(hclust(dist(t(exprs(selected_rma_data)))),xlab="Distance based on expression",main="Hierarchical clustering normalised data")
  dev.off()
 
} 


############## get the control probes ###########
############## load the correct array type:#######

cat(paste0("Loading array type: ", arrayType,"\n"))
if(arrayType=="mouse"){
  controls_info<-read.table("MoGene-2_1-st.qcc",sep="\t", header=T)
}else{
  controls_info<-read.table("HuGene-2_1-st.qcc",sep="\t", header=T)
}

match_control_probe<-match(controls_info$probeset_id, row.names(rma.probe))
controlsLevels<-exprs(rma.probe)[match_control_probe,]
controlsLevels<-cbind(as.data.frame(controlsLevels), controls_info$probeset_name, controls_info$group_name)
colnames(controlsLevels)[length(colnames(controlsLevels))]<-"probe_group"
colnames(controlsLevels)[length(colnames(controlsLevels))-1]<-"probe_type"


############# now get the specific controls
#############BioC_5at
cat("Bioc_5at\n")
BioControls<-controlsLevels[grep("Bio.+5_at", controlsLevels$probe_type),]
CreControls<-controlsLevels[grep("cre-5_at",controlsLevels$probe_type),]
BioControls<-rbind(BioControls,CreControls)
BioControls<-melt(BioControls[,-length(colnames(BioControls))], id="probe_type")
pdf(file="controls/control_spykes.pdf",width = 15, height=7.5)
g<-ggplot(data=BioControls, aes(colour=probe_type,y=value, group=probe_type, x=as.character(variable)))+geom_point()+geom_line()
g<-g+theme(axis.text.x=element_text(angle=90))+labs(x="Samples", y="Intensity", colour="Control")
print(g)
dev.off()

############poly A
cat("PolyA\n")
polyAControls<-controlsLevels[grep("polya", controlsLevels$probe_group),]
polyAControls<-polyAControls[grep("r2.+5.+st", polyAControls$probe_type),]
polyAControls<-melt(polyAControls[,-length(colnames(polyAControls))], id="probe_type")
pdf(file="controls/control_polyA.pdf",width = 15, height=7.5)
g<-ggplot(data=polyAControls, aes(colour=probe_type,y=value, group=probe_type, x=as.character(variable)))+geom_point()+geom_line()
g<-g+theme(axis.text.x=element_text(angle=90))+labs(x="Samples", y="Intensity", colour="Control")
print(g)
dev.off()

#############pos control/neg control
cat("Positive negative controls\n")
posControls<-controlsLevels[grep("pos_control", controlsLevels$probe_group),]
negControls<-controlsLevels[grep("neg_control", controlsLevels$probe_group),]
pos_cont_avg<-colMeans(posControls[, -c(length(colnames(posControls)),length(colnames(posControls))-1)])
pos_cont_avg2<-colMeans(posControls[, -c(length(colnames(posControls)),length(colnames(posControls))-1)])
neg_cont_avg<-colMeans(negControls[, -c(length(colnames(negControls)),length(colnames(negControls))-1)])
neg_cont_avg2<-colMeans(negControls[, -c(length(colnames(negControls)),length(colnames(negControls))-1)])
neg_cont_avg<-cbind(as.data.frame(neg_cont_avg), "neg_control")
pos_cont_avg<-cbind(as.data.frame(pos_cont_avg), "pos_control")
neg_cont_avg<-cbind(as.data.frame(neg_cont_avg),row.names(neg_cont_avg))
pos_cont_avg<-cbind(as.data.frame(pos_cont_avg), row.names(pos_cont_avg))
colnames(neg_cont_avg)<-c("Intensity","Control","SampleID")
colnames(pos_cont_avg)<-c("Intensity","Control", "SampleID")
posneg_contr<-rbind(pos_cont_avg,neg_cont_avg)
pdf(file="controls/control_positive_negative.pdf",width = 15, height=7.5)
g<-ggplot(data=posneg_contr, aes(colour=Control, y=Intensity, x=SampleID,group=Control))+geom_point()+geom_line()
g<-g+theme(axis.text.x=element_text(angle=90))+labs(x="Samples", y="Intensity", colour="Control")
print(g)
dev.off()

############# Adding  signal to noise value to rank table
signal_to_noise<-pos_cont_avg2-neg_cont_avg2
signal_to_noise_rank<-rank(-signal_to_noise)
rank_table<-as.data.frame(signal_to_noise)
rank_table<-cbind(rank_table, signal_to_noise_rank)
colnames(rank_table)[1]<-"Signal to noise"
colnames(rank_table)[2]<-"Signal to noise rank (higher better)"

#fit plm model
pmData<-fitProbeLevelModel(raw.data)

############# Plot NUSE 
cat("Plot NUSE\n")
pdf(file="controls/control_NUSE.pdf",width = 15, height=7.5)
NUSE(pmData, las=2,cex.axis=0.7)
dev.off()

#NUSE table, extract the NUSE value, it's the se divide by its median. SE is extracted when getting the PM data.
#NUSE_table<-sweep(se(pmData),1,rowMedians(se(pmData)),"/")
NUSE_table<-NUSE(pmData,type="value")
NUSE_median<-apply(NUSE_table, 2, "median",na.rm=T)
NUSE_class<-abs(1-NUSE_median)
NUSE_class_rank<-rank(NUSE_class)

########### adding NUSE to rank_table
rank_table<-cbind(rank_table, NUSE_class, NUSE_class_rank)
colnames(rank_table)[3]<-"Nuse control"
colnames(rank_table)[4]<-"Nuse rank lower better"

#############Plot RLE 
cat("Plot RLE\n")
pdf(file="controls/control_RLE.pdf",width = 15, height=7.5)
RLE(pmData,las=2,cex.axis=0.7)
dev.off()

####### adding RLE to rank_table
RLE_table<-RLE(pmData,type="value")
RLE_median<-apply(RLE_table, 2, "median",na.rm=T)
RLE_class<-abs(RLE_median)
RLE_class_rank<-rank(RLE_class)
rank_table<-cbind(rank_table, RLE_class, RLE_class_rank)
colnames(rank_table)[5]<-"RLE control"
colnames(rank_table)[6]<-"RLE rank (lower is better)"


###### adding correlation to rank_table
first_cor=TRUE
for(i in 1: length(group)){
  sample.group<-info[info$Sample.Group%in%group[i],]$Sample.Name
  rma_group<-rma.data.core[,match(sample.group,sampleNames(rma.data.core))]
  cor_group<-cor(exprs(rma_group), method="pearson")
  if(first_cor){
    cor_table<-as.data.frame(rowMeans(cor_group))
    first_cor=FALSE
  }else{
    cor_table<-rbind(cor_table, as.data.frame(rowMeans(cor_group)))
  }
}
index_cor<-match(row.names(rank_table), row.names(cor_table))
cor_table_rank<-rank(-cor_table[index_cor,])
rank_table<-cbind(rank_table, cor_table[index_cor,],cor_table_rank)
colnames(rank_table)[7]<-"Withing group correlation score"
colnames(rank_table)[8]<-"Within group correlation score (higher better)"

##### writting out the rank table ######
rank_table<-cbind(rownames(rank_table), rank_table)
colnames(rank_table)[1]<-"Samples"
write.table(file="Control_ranking_table.txt", rank_table, sep="\t", row.names = F)

#################### plot after normalisation #########
cat("Boxplot after normalisation\n")
pdf(file="rma_plots/boxplot_rma_core.pdf",width=10,height=5)
boxplot(rma.data.core, col=info$Sample.Group,las=2)
dev.off()

cat("HC after normalisation\n")
pdf(file="rma_plots/Hierarchical_clustering_norm_data.pdf",width=10, height=5)
plot(hclust(dist(t(exprs(rma.data.core)))),xlab="Distance based on expression",main="Hierarchical clustering normalised data")
dev.off()

cat("MDS after normalisation\n")
pdf(file="rma_plots/MDS_clustering_rma_data.pdf",width=10, height=5)
fit_mds_rma<-cmdscale(dist(t(exprs(rma.data.core))),eig = T, k=2)
pca1_rma<-fit_mds_rma$points[,1]
pca2_rma<-fit_mds_rma$points[,2]
par(mar=c(5.1,4.1,4.1,12.1),xpd=T)
plot(pca1_rma,pca2_rma,type="n")
text(pca1_rma,pca2_rma,labels=sampleNames(rma.data.core),cex=0.7,col=as.numeric(info$Sample.Group))
legend(legend=levels(info$Sample.Group),fill=c(1:length(levels(info$Sample.Group))), "topright",inset=c(-0.2,0))
dev.off()
par(mar=c(5.1,4.1,4.1,2.1),xpd=F)

cat("Density plot after normalisation\n")
pdf(file="rma_plots/density_plot_rma_data.pdf",width=15, height=5)
par(mar=c(5.1,4.1,4.1,12.1),xpd=T)
hist(rma.data.core, col=1:length(levels(info$Sample.Group)))
legend(legend=levels(info$Sample.Group),fill=c(1:length(levels(info$Sample.Group))), "topright",inset=c(-0.2,0))
dev.off()
par(mar=c(5.1,4.1,4.1,2.1),xpd=F)

cat("MA plot after normalisation\n")
graph_n=0
for(i in 1:length(sampleNames(rma.data.core))){
  current_n=ceiling(i/4)
  if(current_n==graph_n)
  {
    MAplot(rma.data.core, which=i)
  }else{
    if(current_n>1){
      dev.off()
    }
    graph_n=current_n
    jpeg(file=paste("rma_plots/MAplot_",graph_n,"_rma_data.jpg",sep=""),width = 480, height=480)
    par(mfrow=c(2,2))
    MAplot(rma.data.core, which=i)
  } 
}
dev.off()

#############3D MDS #####################
cat("3D MDS after normalisation\n")
fit_mds_rma<-cmdscale(dist(t(exprs(rma.data.core))),eig = T, k=3)
pdf(file="rma_plots/MDS_clustering_norm_data_3D.pdf",width=10, height=5)
scatterplot3d(fit_mds_rma$points[,c(1:3)], xlab="MDS1",ylab="MDS2", zlab="MDS3", color=as.numeric(info$Sample.Group),pch=19, type="h" )
legend("topleft", title="Groups", levels(as.factor(info$Sample.Group)), fill=levels(as.factor(as.numeric(info$Sample.Group))))
dev.off()

#alternative (interactive): plot3d(fit_mds_raw$points[,c(1:3)], col=as.numeric(info$Sample.Group), siz=1, type="s", xlab="PCA1", ylab="PCA2", zlab="PCA3")
#scatterplot3d(pca$rotation[,c(1:3)], xlab="PCA1",ylab="PCA2", zlab="PCA3", color=as.numeric(info$Sample.Group),pch=19, type="h" )

##### coloured dendogram #######
cat("Coloured dendogram after normalisation\n")
pdf(file="rma_plots/Hierarchical_clustering_norm_data_coloured.pdf",width=10, height=5)
hc<-hclust(dist(t(exprs(rma.data.core))))
dend<-as.dendrogram(hc)
dend <- rotate(dend, 1:dim(info)[1])
dend <- color_branches(dend, k=ngroup)
dend <- set(dend, "labels_cex", 0.5)
match_order<-match(labels(dend), info$Sample.Name)
group_order<-info$Sample.Group[match_order]
label_legend<-unique(group_order)
labels_colors(dend) <-  rainbow_hcl(ngroup)[sort_levels_values(as.numeric(info$Sample.Group)[order.dendrogram(dend)])]
plot(dend, main="Hierarchical cluster of sample groups", horiz=T, nodePar=list(cex=.007))
legend("topleft", title="Groups", as.character(label_legend), fill=rainbow_hcl(ngroup))
dev.off()

# Comparisons ####
cat("Performing comparisons\n")
summary_sig<-data.frame(Comparisons=character(), Significant_0.05=integer(),Significant_0.01=integer(),stringsAsFactors = F)
n_comp=0
first_comp=TRUE
for(i in 1:length(group)){
  for( l in i+1:(length(group)-1)){
    if(l<=length(group)){
      n_comp=n_comp+1
      group1=group[i]    
      group2=group[l]
      cat(paste0("Comparing groups: ",group[i]," v ",group[l],"\n"))
      index1<-row.names(info[info$Sample.Group==group1,])
      index2<-row.names(info[info$Sample.Group==group2,])
      index_all<-c(index1,index2)
      selected_samples<-(info[index_all,])
      selected_samples<-droplevels(selected_samples)
      raw.data.s<-raw.data[,as.numeric(index_all)]
      rma.data.core.s<-rma(raw.data.s, target="core")
      rma.data.core.s<-rma.data.core.s[,as.character(selected_samples$Sample.Name)]
      median_exprs<-median(exprs(rma.data.core.s))
      
      ####Filter expression levels#
      
      filter_setting<-kOverA(length(sampleNames(rma.data.core.s))/2,median_exprs)
      filter_function<-filterfun(filter_setting)
      probe_filter<-genefilter(exprs(rma.data.core.s),filter_function)
      rma.data.core.s<-rma.data.core.s[probe_filter,]
      data.annot<-getNetAffx(rma.data.core.s, type="transcript")
      fData(rma.data.core.s)<-pData(data.annot)
      rma.data.core.s<-rma.data.core.s[fData(rma.data.core.s)$category=="main",]
      design<-model.matrix(~0+selected_samples$Sample.Group)
      colnames(design)<-gsub(colnames(design), pattern = "selected_samples.Sample.Group",replacement = "")
      
      #### generate graphs during analysis ###
      for (k in 1:2){
        group_selected<-group[k]
        current_samples<-info[info$Sample.Group==current_group,]
        color_selected<-rainbow_hcl(length(current_samples$Sample.Name))
        if(dim(selected_rma_data)[2]>3){
          pdf(file=paste("results/graphs/MDS_clustering_rma_data_3D_",group1,"_v_",group2,"_comp",".pdf",sep=""),width=10, height=5)
          par(xpd = T)
          fit_mds_rma<-cmdscale(dist(t(exprs(selected_rma_data))),eig = T, k=3)
          write.table(fit_mds_rma$points,file=paste("rma_plots/MDS_clustering_rma_data_table_",group1,"_v_",group2,".txt",sep=""), sep="\t" )
          scatterplot3d(fit_mds_rma$points[,c(1:3)], xlab="MDS1",ylab="MDS2", zlab="MDS3",pch=19, type="h",color=color_selected[1:length(current_samples$Sample.Name)])
          legend("topleft", title="Samples", legend = current_samples$Sample.Name, fill=color_selected, cex = 0.7, inset=c(-0.05,-0.1))
          dev.off()
        }
        pdf(file=paste("results/graphs/MDS_clustering_rma_data_",group1,"_v_",group2,"_comp",".pdf",sep=""),width=10, height=5)
        par(xpd = T)
        fit_mds_rma<-cmdscale(dist(t(exprs(rma.data.core.s))),eig = T, k=2)
        mds1<-fit_mds_rma$points[,1]
        mds2<-fit_mds_rma$points[,2]
        plot(mds1, mds2, main="MDS plot", xlab="MDS 1", ylab="MDS 2", col=rainbow_hcl(length(current_samples$Sample.Name)), pch=16)
        legend("topleft", title="Samples", legend = current_samples$Sample.Name, fill=rainbow_hcl(length(current_samples$Sample.Name)), cex=0.7, inset(0,-0.1))
        dev.off()
        pdf(file=paste("results/graphs/Hierarchical_clustering_rma_data_",group1,"_v_",group2,"_comp",".pdf",sep=""),width=10, height=5)
        plot(hclust(dist(t(exprs(rma.data.core.s)))),xlab="Distance based on expression",main="Hierarchical clustering normalised data")
        dev.off()
      }
      
      
      
      fit<-lmFit(exprs(rma.data.core.s),design)
      contrast.matrix=makeContrasts(paste(as.character(group2),as.character(group1),sep="-"),levels=design)
      fit.contrast<-contrasts.fit(fit, contrast.matrix)
      fit.ebayes<-eBayes(fit.contrast)
      results.table<-topTable(fit.ebayes,number=dim(exprs(rma.data.core.s))[1],adjust.method="BH")
      results.ids<-row.names(results.table)
      results<-cbind(fData(rma.data.core.s[results.ids,]),results.table,exprs(rma.data.core.s)[results.ids,])
      write.table(results, file=paste("results/tables/results_",group1,"_vs_",group2,".csv",sep=""), sep=",", row.names = F)
      logFC<-results.table$logFC
      log_adj.pvalue=-log(results.table$adj.P.Val,base = 10)
      significant<-results.table$adj.P.Val<0.05
          
      n_sig_0.05<-dim(results[results.table$adj.P.Val<0.05,])[1]
      n_sig_0.01<-dim(results[results.table$adj.P.val<0.01,])[1]
      comp<-paste(group1,"_vs_",group2,sep="")
      cat(paste0(n_comp," ",comp, " ",n_sig_0.05," ",n_sig_0.01,"\n"))
      if (first_comp){
        summary_sig[n_comp,]<-c(comp,n_sig_0.05,n_sig_0.01)
        first_comp=FALSE
      }else{
        summary_sig<-rbind(summary_sig, c(comp, n_sig_0.05, n_sig_0.01))
      }
      print(summary_sig)
          
      volcano<-as.data.frame(cbind(logFC,log_adj.pvalue))
      g<-ggplot(volcano, aes(logFC,log_adj.pvalue,col=significant))+geom_point()+ggtitle(paste("Volcano_",group1,"_vs_",group2,sep=""))+ylab("significance level (-log(adj.pval))")
      pdf(file=paste("results/graphs/Volcano_",group1,"_vs_",group2,".pdf",sep=""),width = 10, height=5)
      print(g)
      dev.off()
    }
  }
}
write.table(file="results/summary_comparisons.txt",summary_sig, row.names = F, sep="\t")



  