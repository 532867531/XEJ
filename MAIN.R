dir="D:/min数据/浙江大学XEJ的5个小鼠肝脏普通转录组建库测序分析结题报告/upload/4_2_singleDE_table"
setwd(dir)
list.files(pattern = "*.csv")
(files_csv=list.files(pattern = "*.csv"))
all_files=list()
for(onefile in files_csv){
  all_files[[onefile]]=read.csv(file = onefile,stringsAsFactors = FALSE)
}
for(onefile in names(all_files)){
  colnames(all_files[[onefile]])[which(grepl(pattern = "ENSMUSG.*",x = all_files[[onefile]][10,],perl = TRUE,ignore.case = T))]='Gene'
}

big_all_file=all_files[[1]]
for(onefile in names(all_files)[c(2:length(names(all_files)))]){
    big_all_file=merge(x = big_all_file,y = all_files[[onefile]],by.x='Gene',by.y='Gene',all = T)
}


big_all_file1=big_all_file[,c(1,12,which(grepl(pattern = ".*count.*|.*normalize.*",x = colnames(big_all_file),ignore.case = T,perl = T)))]
big_all_file1=big_all_file1[,order(colnames(big_all_file1))]
big_all_file1[,"Gene"]=paste(big_all_file1$Gene,big_all_file1$GeneName.x,sep="_")
rownames(big_all_file1)=big_all_file1$Gene
##去除全是NA的列
for(rowindex in c(1:nrow(big_all_file1))){
  big_all_file1[rowindex,which(is.na(big_all_file1[rowindex,]))]=0
}


big_all_file2=big_all_file1[,c("C1_count.x","C1_normalize.x","T1_count.x","T1_normalize.x","T2_count.x","T2_normalize.x","T3_count","T3_normalize","T4_count","T4_normalize")]
big_all_file3=big_all_file2[apply(X = big_all_file2,MARGIN = 1,FUN = function(x){sum(as.numeric(x[which(grepl(pattern = "count",x = colnames(big_all_file2),ignore.case = T,perl = T))])<=10)<4}),c("C1_normalize.x","T1_normalize.x","T2_normalize.x","T3_normalize","T4_normalize")]
(colnames(big_all_file3)=stringi::stri_extract_all(str = colnames(big_all_file3),regex = "\\w+_normalize"))

##
mat=as.matrix(big_all_file3)
library(NGCHM)
mat=t(scale(t(mat),center = T,scale = F))
hm1 <- chmNew ('my-first-ngchm_YX',mat)
setwd("F:/min-labs_paper/work/XEJ")
chmExportToPDF(chm = hm1, filename = "temp.pdf",shaidyMapGen = "F:/min-labs_paper/work/WHX/ShaidyMapGen.jar",overwrite = T)
chmExportToFile(hm1,paste('YANGXIANG','.ngchm',sep = "_"),shaidyMapGen = "F:/min-labs_paper/work/YANGXIANG/ShaidyMapGen.jar",overwrite = T)


