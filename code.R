#############################m6A调节子与lncRNA共表达#################################
library(limma)
corFilter=0.4            
pvalueFilter=0.001       
setwd("D:\\TCGA")   

rt=read.table("lncRNA.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.15,]

group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
lncRNA=data[,group==0]
conNum=length(group[group==1])       
treatNum=length(group[group==0])     
sampleType=c(rep(1,conNum), rep(2,treatNum))

rt1=read.table("m6aGeneExp.txt", header=T, sep="\t", check.names=F)
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames1=list(rownames(exp1),colnames(exp1))
m6A=matrix(as.numeric(as.matrix(exp1)), nrow=nrow(exp1), dimnames=dimnames1)
m6A=avereps(m6A)
m6A=m6A[rowMeans(m6A)>0.15,]

group=sapply(strsplit(colnames(m6A),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
m6A=m6A[,group==0]

outTab=data.frame()
for(i in row.names(lncRNA)){
	if(sd(lncRNA[i,])>0.1){
		test=wilcox.test(data[i,] ~ sampleType)
		if(test$p.value<0.05){
			for(j in row.names(m6A)){
				x=as.numeric(lncRNA[i,])
				y=as.numeric(m6A[j,])
				corT=cor.test(x,y)
				cor=corT$estimate
				pvalue=corT$p.value
				if((cor>corFilter) & (pvalue<pvalueFilter)){
					outTab=rbind(outTab,cbind(m6A=j,lncRNA=i,cor,pvalue,Regulation="postive"))
				}
				if((cor< -corFilter) & (pvalue<pvalueFilter)){
					outTab=rbind(outTab,cbind(m6A=j,lncRNA=i,cor,pvalue,Regulation="negative"))
				}
			}
		}
	}
}

write.table(file="net.network.txt",outTab,sep="\t",quote=F,row.names=F)

lncNode=data.frame(Node=unique(as.vector(outTab[,"lncRNA"])), Type="lncRNA")
mrnaNode=data.frame(Node=unique(as.vector(outTab[,"m6A"])), Type="m6A")
nodeOut=rbind(lncNode, mrnaNode)
write.table(nodeOut, file="net.node.txt", sep="\t", quote=F, row.names=F)

m6aLncRNA=unique(as.vector(outTab[,"lncRNA"]))
m6aLncRNAexp=data[m6aLncRNA,]
m6aLncRNAexp=rbind(ID=colnames(m6aLncRNAexp), m6aLncRNAexp)
write.table(m6aLncRNAexp,file="m6aLncExp.txt",sep="\t",quote=F,col.names=F)

####################################三分型免疫疗效violin图#######################################
library(limma)
m6aCluFile="mirlncRNACluster.txt"        
ExpressionFile="ips_ctla4_pos_pd1_pos.txt"           
setwd("C:\\Users\\11976\\Desktop\\wenxian_G\\7ips\\1")    

m6aClu=read.table(m6aCluFile, header=T, sep="\t", check.names=F, row.names=1)
dim(m6aClu)
Expression=read.table(ExpressionFile, header=T, sep="\t", check.names=F, row.names=1)
dim(Expression)

sameSample=intersect(row.names(m6aClu), row.names(Expression))
dim(sameSample)
data=cbind(Expression[sameSample,,drop=F], m6aClu[sameSample,,drop=F])
write.table(data,file="ips_ctla4_pos_pd1_pos.txt",sep="\t",quote=F,col.names=T)
############把输出的文件列明改成"id ips cluster"

library(ggpubr)           
inputFile="ips_ctla4_pos_pd1_pos.txt"      
outFile="ips_ctla4_pos_pd1_pos.pdf"      
setwd("C:\\Users\\11976\\Desktop\\wenxian_G\\7ips\\1")    

rt=read.table(inputFile,header=T,sep="\t",check.names=F)
x=colnames(rt)[3]
y=colnames(rt)[2]
colnames(rt)=c("id","Expression","Type")

group=levels(factor(rt$Type))
rt$Type=factor(rt$Type, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

pdf(file=outFile, width=6, height=5)
ggviolin(rt, x="Type", y="Expression", fill = "Type", 
         xlab=x, ylab=y,
         legend.title=x,
         add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons)
#stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()