library(RMySQL)
library(dplyr)

# 创建数据库连接
con.lihc = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.lihc', user = 'root', password = 'root')
con.gene = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'annotation', user = 'root', password = 'root')

# 设置编码格式utf8
dbSendQuery(con.gene,'SET NAMES utf8')  

# collagen gene name from wikipedia
v.col.gene = c("COL1A1","COL1A2","COL2A1","COL3A1","COL4A1","COL4A2","COL4A3BP","COL4A3","COL4A4","COL4A5","COL4A6",
               "COL5A1","COL5A2","COL5A3","COL6A1","COL6A2","COL6A3","COL6A4P2","COL6A6","COL7A1","COL8A1","COL8A2",
               "COL9A1","COL9A2","COL9A3","COL10A1","COL11A1","COL11A2","COL12A1","COL13A1","COL14A1","COL15A1",
               "COL16A1","COL17A1","COL18A1","COL19A1","COL20A1","COL21A1","COL22A1","COL23A1","COL24A1","COL25A1",
               "EMID2","COL27A1","COL28A1","COL29A1");

d.gene.hsa = dbReadTable(con.gene, name = 'id_symbol')
d.col.gene = d.gene.hsa[ d.gene.hsa$SYMBOL %in% v.col.gene, ]

dbWriteTable(con.gene, name = 'collagen_sys', value = d.col.gene, row.names = F, overwrite = F)


# Collagen from NCBI
d.col.ncbi = read.delim(file = './data/Collagen.gene.list.NCBI.txt')
dbWriteTable(con.gene, name = 'collagen_ncbi', value = d.col.ncbi, row.names = F, overwrite = F)


# Collagen From genecard
d.col.genecard = read.csv(file = "./data/Collagen.gene.GeneCards.csv", stringsAsFactors = F)[,-1]
dbWriteTable(con.gene, name = 'collagen_genecard', value = d.col.genecard, row.names = F, overwrite = F)

# test Collagen Related Gene
length(intersect(d.col.gene$SYMBOL, d.col.ncbi$Symbol))
length(intersect(d.col.gene$ID, d.col.ncbi$GeneID))

length(intersect(d.col.genecard$ENTREZID, d.col.ncbi$GeneID))

which( d.col.ncbi$GeneID  %in% d.col.genecard$ENTREZID )
which( d.col.genecard$ENTREZID  %in% d.col.ncbi$GeneID)

# get Collagen Related Gene
# A = GeneCard gene relevance.score > 10 
# insect between A and B(NCBI gene)  
d.col.genecard = d.col.genecard[which(d.col.genecard$Relevance.score >= 10), ]
v.col.union = intersect(d.col.genecard$ENTREZID, d.col.ncbi$GeneID)
v.col.union = union(v.col.union, d.col.gene$ID)

d.col.gene = d.gene.hsa[ d.gene.hsa$ID %in% v.col.union, ]
dbWriteTable(con.gene, name = 'collagen_related', value = d.col.gene, row.names = F, overwrite = F)














str.a = 'NM_0001,NM_0002,NM_0003,NM_0004,NM_0005,NM_0006,NM_0007,NM_0008'
result = as.character(strsplit(str.a, split = ',')[[1]])


setwd("E:/CODE/Tumor Grade/CollagenAnalysis/data")
in.data = read.delim(file = './GPL123.txt', stringsAsFactors = F)

# 先去掉GB_List 为空的项
in.data = in.data[nchar(in.data$GB_LIST) > 0,]
v.ID = c()
v.GB = c()
for (i.data in in.data)
{
   GB.list = as.character(strsplit(i.data[2], split = ',')[[1]])
   ID.list = rep(i.data[1], length(GB.list))
   
   if(length(GB.list) == 0) next
   
   v.ID = c(v.ID, ID.list)
   v.GB = c(v.GB, GB.list)
}
out.data = data.frame(ID = v.ID, GB = v.GB)

a = ""
nchar(a)












