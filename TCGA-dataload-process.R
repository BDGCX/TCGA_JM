###生物信息云
rm(list=ls())
options(stringsAsFactors = F)

##首先在GDC网页在线下载Manifest,(Metadata)json和Cart文件,将Cart文件解压缩
##数据下载

##数据处理#####
##将下载的所有文件移动到同一个文件夹SampleFiles中
dir.create("SampleFiles")#创建工作目录
filepath <- dir(path = "./Gene_Expression_Quantification",
                full.names = T) 
for (wd in filepath) {
  files <- dir(path = wd,pattern = "gz$")#查看满足条件文件
  fromfilepath <- paste(wd,"/",files,sep = "")
  tofilepath <- paste("./SampleFiles/",files,sep = "")
  file.copy(fromfilepath,tofilepath)
  
}

####解压所有文件并删除原文件#######
setwd("./SampleFiles")
countsFiles <- dir(path = "./",pattern = "gz$")#查看满足条件文件
length(countsFiles)#查看文件数量
library(R.utils)
sapply(countsFiles, gunzip)#解压函数gunzip需要R.utils包

####处理json文件######
library(rjson)
metadata_json_File <- fromJSON(file="../metadata.cart.time.json")
View(metadata_json_File)
json_File_Info <- data.frame(fileNames=c(),TCGA_Barcode=c())
for (i in 1:length(metadata_json_File)) {
  TCGA_Barcode <- metadata_json_File[[i]][["associated_entities"]][[1]][["entity_submitted_id"]]
  file_name <- metadata_json_File[[i]][["file_name"]]
  json_File_Info <- rbind(json_File_Info,
                          data.frame(fileNames=file_name),TCGA_Barcode=TCGA_Barcode)
  
}
rownames(json_File_Info) <- json_File_Info[,1]
write.csv(json_File_Info,file = "../json_File_Info.csv")

####获取Counts矩阵######
fileName_To_TCGA_BarcodeFile <- json_File_Info[-1]
countsFileNames <- dir(pattern = "counts$")#list.files函数也可以用

allSampleRawCounts <- data.frame()
for (txtfile in countsFileNames) {
  #每一个循环读取一个文件
  SampleCounts <- read.table(txtfile,header = F)
  rownames(SampleCounts) <- SampleCounts[,1]
  SampleCounts <- SampleCounts[-1]#colnames出错的话SampleCounts定义成数据框
  ##根据fileName_To_TCGA_BarcodeFile文件中文件名称与Barcode对应关系，命名列名
  colnames(SampleCounts) <- fileName_To_TCGA_BarcodeFile[paste(txtfile,".gz",sep = ""),]
  if (dim(allSampleRawCounts)[1]==0) {
    allSampleRawCounts <- SampleCounts
  }
  else{allSampleRawCounts <- cbind(allSampleRawCounts,SampleCounts)}
}
write.csv(allSampleRawCounts,file = "../allSampleRawCounts.csv")
ensembl_id <- substr(row.names(allSampleRawCounts),1,15)
rownames(allSampleRawCounts) <- ensembl_id
#RawCounts.csv文件与allSampleRawCounts.csv的区别在于行名的ensembl去掉了版本号
write.csv(allSampleRawCounts,file = "../RawCounts.csv")

#####ID转换#####
#添加一列Ensembl_ID到RawCounts数据框中
RawCounts <- allSampleRawCounts
Ensembl_ID <- data.frame(Ensembl_ID=rownames(RawCounts))
rownames(Ensembl_ID) <- Ensembl_ID[,1]
RawCounts <- cbind(Ensembl_ID,RawCounts)


#一个函数，通过gtf文件获取Ensembl_ID与基因名称的对应关系
get_map <- function(input){
  if(is.character(input)){
    if (!file.exists(input)) stop("Bad input file.")
    message("Treat input as file")
    input=data.table::fread(input,header = F)
  } else{
        data.table::setDT(input)
  }
  
  input=input[input[[3]]=="gene",]
  
  pattern_id=".*gene_id\"([^;]+)\";.*"
  pattern_name=".*gene_name\"([^;]+)\";.*"
  
  gene_id=sub(pattern_id,"\\1",input[[9]])
  gene_name=sub(pattern_name,"\\1",input[[9]])
  
  Ensembl_ID_To_Genename <- data.frame(gene_id=gene_id,
                                       gene_name=gene_name,
                                       stringsAsFactors = F)
  return(Ensembl_ID_To_Genename)
  
}

Ensembl_ID_To_Genename <- get_map("../gencode.v33lift37.annotation.gtf")

gtf_Ensembl_ID <- substr(Ensembl_ID_To_Genename[,1],1,15)
Ensembl_ID_To_Genename <- data.frame(gtf_Ensembl_ID,Ensembl_ID_To_Genename[,2])
colnames(Ensembl_ID_To_Genename) <- c("Ensembl_ID","gene_id")
write.csv(Ensembl_ID_To_Genename,file = "../Ensembl_ID_To_Genename.csv")

#融合数据
mergeRawCounts <- merge(Ensembl_ID_To_Genename,RawCounts,by="Ensembl_ID")

#按照gene_id列进行排序
mergeRawCounts <- mergeRawCounts[order(mergeRawCounts[,"gene_id"]),]
#按照gene_id列建立索引
index <- duplicated(mergeRawCounts$gene_id)
#我们想要的那一行为FALSE,素以要取反
mergeRawCounts <- mergeRawCounts[!index,]
#利用基因名称作为行名
rownames(mergeRawCounts) <- mergeRawCounts[,"gene_id"]
#删除前两列
LUAD_Counts_expMatrix <- mergeRawCounts[,-c(1:2)]
#保存文件
write.csv(LUAD_Counts_expMatrix,file = "../LUAD_Counts_expMatrix.csv")








