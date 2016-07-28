# ***************************************
# * Name: liangsen_2014@163.com
# * Date: Mon Jul 11 19:56:35 2016
# * Description: null 
# *              
# ***************************************

# get ensembl.id vector from orignal TCGA matrix data
v.ensembl.id.dot = rownames(KIRC.FPKM)
v.ensembl.id = sapply(v.ensembl.id.dot, substr,start = 1,stop = 15, USE.NAMES = F)

# map ensembl id to gene symbol and gene entrz_id

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("mygene")

library(mygene)

mapgens = getGenes(v.ensembl.id, scope = 'reporter',fields = c("entrezgene","ensemblgene","symbol","name"), species = 'human', as.return = 'Data.Frame')
save(mapgens, file = 'E:/DATA/annotation/ensembl.map.gene.60000.RDATA')
save(v.ensembl.id,v.ensembl.id.dot, file = 'E:/DATA/annotation/ensembl.gene.60000.RDATA')

entrz.gene = as.data.frame(mapgens[which(mapgens$entrezgene > 0),])


# get annotation from GENECODE website
annotation = read.csv(file = "E:/DATA/annotation/gencode.v24.basic.annotation.gtf", sep = "\t", stringsAsFactors = F, header = F)
colnames(annotation) = c("chromosome_name", "annotation_source","feature_type","genomic_start_location",
                         "genomic_end_location","score","genomic_strand","genomic_phase","additional_info")


v.anno.ensembl = rep(NaN, nrow(annotation))
v.anno.gene.type = rep(NaN, nrow(annotation))

annotation.gene = annotation[which(annotation$V3 == "gene"),]

hm.annotation.grc38 = data.frame(matrix(NA, nrow = nrow(annotation.gene), ncol = 9),stringsAsFactors = F)
colnames(hm.annotation.grc38) = c("chr","Source","feature","gene_start","gene_end","gene_id","gene_type","gene_name","gene_status")

nrow(annotation)


attach(annotation.gene)
for(i in 1:nrow(annotation.gene)){
  info = strsplit(annotation.gene[i,9],split = "; ")[[1]]
  info.list = list();
  for(j in 1:length(info)){
    a.info = strsplit(info[j], split = ' ')[[1]]
    info.list[[a.info[1]]] = a.info[2]
  }

  hm.annotation.grc38[i,] = c(annotation.gene[i,1:5],info.list$gene_id,info.list$gene_type,info.list$gene_name,info.list$gene_status);
}
detach(annotation.gene)



annotation.gene.protein.coding = hm.annotation.grc38[which(hm.annotation.grc38$gene_type == 'protein_coding'),]

v.gene.type = unique(hm.annotation.grc38$gene_type)

# clear away dot
hm.annotation.grc38$gene_id = sapply(hm.annotation.grc38$gene_id, substr,start = 1,stop = 15, USE.NAMES = F)
rownames(hm.annotation.grc38) = hm.annotation.grc38$gene_id


mapgens.uni = unique(x = mapgens)

# entrze.gene dataframe
v.ensembl.gene.type = hm.annotation.grc38[v.ensembl.id, 'gene_type']

write.csv(v.ensembl.id, file = "E:/DATA/annotation/v.ensembl.id.csv")
write.csv(v.ensembl.id, file = "E:/DATA/annotation/map")
                             
                             
d.ensembl.gene = data.frame(ensembl_gene=v.ensembl.id, gene_type = v.ensembl.gene.type, stringsAsFactors = F)
table(d.ensembl.gene$gene_type)

dv.1 = read.csv(file = "E:/DATA/annotation/enseml-2-entrz-id.txt", header = T, sep = "\t")
dv.2 = read.csv(file = "E:/DATA/annotation/enseml-2-offical-symbol.txt", header = T, sep = "\t")
d.ensembl.dv = data.frame(ensembl_id = dv.1$From, entrez_id = dv.1$To, symbol = dv.2$To, species = dv.1$Species, gene_name = dv.1$Gene.Name, stringsAsFactors = F)
save(d.ensembl.dv, file = 'E:/DATA/annotation/ensembl.id.david.RDATA')


# save data to mysql
con.gene = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'annotation', user = 'root', password = 'root')
dbWriteTable(con.gene, name = "ensembl_id_david",value = d.ensembl.dv, row.names = F, overwrite = F)
dbWriteTable(con.gene, name = "ensembl_id_grc38",value = annotation, row.names = F, overwrite = F)
dbWriteTable(con.gene, name = "ensembl_id_tcga", value = data.frame(ensembl_id = v.ensembl.id), row.names = F, overwrite = F)

############### 2016-07-13
v.grc38.gene.protein.coding = sapply(annotation.gene.protein.coding$gene_id, substr,start = 1,stop = 15, USE.NAMES = F)


hm.annotation.grc38.protein.coding = hm.annotation.grc38[v.grc38.gene.protein.coding,]
write.csv2(x = v.grc38.gene.protein.coding, file = "E:/DATA/annotation/v.grc38.gene.protein.coding.csv", row.names = F)



read.csv(file = "E:/DATA/annotation/mart_export.txt", header = T, sep = "\t" )

############### 2016-07-14
v.tcga.grc38.protein.coding = v.ensembl.id[v.ensembl.id %in% v.grc38.gene.protein.coding]
write.csv2(x = v.tcga.grc38.protein.coding, file = "E:/DATA/annotation/v.tcga.grc38.protein.coding.csv", row.names = F, col.names = F, quote = F)

### from uniport
# input v.tcga.grc38.protein.coding
d.uniport.gene = read.csv(file = "E:/DATA/annotation/uniprot-yourlist.tab", header = T, sep = "\t",)
v.uniport.gene = d.uniport.gene$yourlist

v.uniport.gene.u = unique(v.uniport.gene)


# from db2db
d.db2db.gene = read.csv(file = "E:/DATA/annotation/bioDBnet_db2db_160713222539_374717289.txt", header = T, sep = "\t",stringsAsFactors = F)[,-3]

library(dplyr)
# 三类�?(1)一对多的multi�?(2)映射不了的null，（3）多对一�?
d.db2db.gene = data.frame(d.db2db.gene)
#d.db2db.gene.multi = filter(d.db2db.gene, grep(";", d.db2db.gene$Gene.ID))

d.db2db.gene.multi = d.db2db.gene[grep(";", d.db2db.gene$Gene.ID),]
d.db2db.gene.null = d.db2db.gene[grep("-", d.db2db.gene$Gene.ID),]

v.db2db.gene.multi.null = c(d.db2db.gene.null$Ensembl.Gene.ID, d.db2db.gene.multi$Ensembl.Gene.ID)

v.db2db.gene.avai = setdiff(d.db2db.gene$Ensembl.Gene.ID,v.db2db.gene.multi.null)
d.db2db.gene.avai = d.db2db.gene[d.db2db.gene$Ensembl.Gene.ID %in% v.db2db.gene.avai,]

# save to mysql
con.gene = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'annotation', user = 'root', password = 'root')
dbWriteTable(con.gene, name = "tcga_db2db_gene",value = d.db2db.gene, row.names = F, overwrite = F)
dbWriteTable(con.gene, name = "tcga_db2db_gene_multi",value = d.db2db.gene.multi, row.names = F, overwrite = F)
dbWriteTable(con.gene, name = "tcga_db2db_gene_null",value = d.db2db.gene.null, row.names = F, overwrite = F)
dbWriteTable(con.gene, name = "tcga_db2db_gene_avai",value = d.db2db.gene.avai, row.names = F, overwrite = F)

write.csv(x = d.db2db.gene.multi, file = "E:/DATA/annotation/d.db2db.gene.multi.csv", row.names = F, quote = F)
write.csv2(x = d.db2db.gene.null, file = "E:/DATA/annotation/d.db2db.gene.null.csv", row.names = F, quote = F)

############## 2016-07-15

# gene map to entrz ID
d.gene.map = d.db2db.gene
d.gene.map.multi = read.csv(file = "E:/DATA/annotation/d.db2db.gene.multi.csv", header = T, sep = ",",stringsAsFactors = F)
d.gene.map.null = read.csv(file = "E:/DATA/annotation/d.db2db.gene.null.csv", header = T, sep = ";",stringsAsFactors = F)

rownames(d.gene.map.multi) = d.gene.map.multi$Ensembl.Gene.ID

for(i in 1:nrow(d.gene.map)) {
  if(d.gene.map$Ensembl.Gene.ID[i] %in% d.gene.map.null$Ensembl.Gene.ID) {
    d.gene.map$Gene.ID[i] = d.gene.map$Ensembl.Gene.ID[i]
  } 
  
  else if(d.gene.map$Ensembl.Gene.ID[i] %in% d.gene.map.multi$Ensembl.Gene.ID) {
    d.gene.map$Gene.ID[i] = d.gene.map.multi[d.gene.map$Ensembl.Gene.ID[i],"Gene.ID"]
  }
}

# save to mysql
con.gene = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'annotation', user = 'root', password = 'root')
dbWriteTable(con.gene, name = "tcga_gene_map",value = d.gene.map, row.names = F, overwrite = F)
dbWriteTable(con.gene, name = "tcga_gene_map_multi",value = d.gene.map.multi, row.names = F, overwrite = F)
dbWriteTable(con.gene, name = "tcga_gene_map_null",value = d.gene.map.null, row.names = F, overwrite = F)

dbWriteTable(con.gene, name = "ensembl_id_grc38_gene", value = annotation.gene, row.names = F, overwrit = F)
dbWriteTable(con.gene, name = "ensembl_id_grc38_gene_protein_coding", value = annotation.gene.protein.coding, row.names = F, overwrit = F)


save(d.gene.map,d.gene.map.multi,d.gene.map.null, file = "E:/DATA/annotation/tcga.gene.map.RData")

########## make sure 
gene.id.multi = d.gene.map[duplicated(d.gene.map$Gene.ID),]
gene.id.multi = gene.id.multi[order(gene.id.multi$Gene.ID),]

v.gene.id.multi = unique(gene.id.multi$Gene.ID)

v.gene.id.2 = set(unique(d.gene.map$Gene.ID),d.gene.map$Gene.ID)
v.gene.id.multi
a = duplicated(d.gene.map$Gene.ID)

d.gene.map.null = read.csv(file = "E:/DATA/annotation/d.db2db.gene.null.csv", header = T, sep = ";",stringsAsFactors = F)


d.gene.map.id.not.unique = d.gene.map[d.gene.map$Gene.ID %in% v.gene.id.multi,]
d.gene.map.id.unique = d.gene.map[!d.gene.map$Gene.ID %in% v.gene.id.multi,]

dbWriteTable(con.gene, name = "tcga_gene_map_id_not_unique2",value = d.gene.map.id.not.unique2, row.names = F, overwrite = F)
dbWriteTable(con.gene, name = "tcga_gene_map_id_not_unique",value = d.gene.map.id.not.unique, row.names = F, overwrite = F)
dbWriteTable(con.gene, name = "tcga_gene_map_id_unique",value = d.gene.map.id.unique, row.names = F, overwrite = F)

write.csv(v.gene.id.multi, file = "E:/DATA/annotation/v.gene.id.multi.csv", row.names = F, quote = F)
write.csv(d.gene.map.id.unique, file = "E:/DATA/annotation/d.gene.map.id.unique.csv", row.names = F, quote = F)
write.csv(d.gene.map.id.not.unique, file = "E:/DATA/annotation/d.gene.map.id.not.unique.csv", row.names = F, quote = F)

d.gene.map.id.not.unique2 = read.csv(file = "E:/DATA/annotation/d.gene.map.id.not.unique.csv", header = T, sep = ",",stringsAsFactors = F)
d.gene.map.all = rbind(d.gene.map.id.unique, d.gene.map.id.not.unique2)

d.gene.map.all = d.gene.map.all[order(d.gene.map.all$Gene.ID),]
rownames(d.gene.map.all) = d.gene.map.all$Gene.ID
dbWriteTable(con.gene, name = "tcga_gene_map_id_ok",value = d.gene.map.all, row.names = F, overwrite = F)

tcga.gene.map = tcga.gene.map[order(tcga.gene.map$Gene.ID),]
rownames(d.gene.map.all) = d.gene.map.all$Gene.ID
dbWriteTable(con.gene, name = "tcga_gene_map",value = tcga.gene.map, row.names = F, overwrite = T)

######################### 2016-07-28
con.gene = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'annotation', user = 'root', password = 'root')
save(d.gene.map,d.gene.map.multi,d.gene.map.null, file = "E:/DATA/annotation/tcga.gene.map.RData")

