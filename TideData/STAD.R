# ***************************************
# * Name: liangsen_2014@163.com
# * Date: Sun Jul 10 19:51:20 2016
# * Description: STAD
# *              
# ***************************************

# define two file paths
datapath = 'E:/DATA/'

# read manifest file
d.manifest = read.csv(file = 'E:/DATA/TCGA-manifest-file-14/STAD-gdc_manifest.2016-07-10T13-24-40.103392.tsv', sep = '\t',header = T)

# get a test files of XXX.FPKM.txt.gz
FPKM.gz.path = paste(datapath, 'TCGA-Data-14-20160704/',d.manifest$id[1],'/',d.manifest$filename[1], sep = '')
testdata = read.csv(file = FPKM.gz.path, sep = '\t', header = F)

# define a data.frame for all data 
STAD.FPKM = data.frame(matrix(NA, nrow = nrow(testdata), ncol = nrow(d.manifest)), stringsAsFactors = F)
colnames(STAD.FPKM) = apply(X = d.manifest, MARGIN = 1, FUN = function(aRow){substr(aRow[2], 1, 36)})
rownames(STAD.FPKM) = testdata[,1]

# get all data form files
for(i.id in 1:nrow(d.manifest))
{
  
  FPKM.gz.path = paste(datapath, 'TCGA-Data-14-20160704/',d.manifest$id[i.id],'/',d.manifest$filename[i.id], sep = '')
  STAD.FPKM[,i.id] = read.csv(file = FPKM.gz.path, sep = '\t', header = F)[,2]
  
  # remove file to another directions
  file.copy(from = FPKM.gz.path, to = 'E:/DATA/TCGA-Data/STAD/',overwrite = F)
  print(i.id)
}

save(STAD.FPKM, file = 'E:/DATA/TCGA-RDATA/STAD.FPKM.RData')

########### 2016-07-15
load('E:/DATA/TCGA-RDATA/STAD.FPKM.RData')

#save to mysql
con.stad = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.stad', user = 'root', password = 'root')

dbWriteTable(con.stad, name = "tcga_stad_fpkm", value = STAD.FPKM, row.names = T, overwrite = F)

# get tcga gene list
#con.anno = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'annotation', user = 'root', password = 'root')
#tcga.gene.list = dbReadTable(conn = con.anno, name= "ensembl_id_tcga")[,1]

# get protein_coding list
#tcga.gene.map = dbReadTable(conn = con.anno, name= "tcga_gene_map")

# set name
rownames(STAD.FPKM) = tcga.gene.list

# filter
STAD.FPKM.GENE = STAD.FPKM[tcga.gene.map$Ensembl.Gene.ID,] 
rownames(STAD.FPKM.GENE) = tcga.gene.map$Gene.ID

# save to mysql and RDATA
#con.stad = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.stad', user = 'root', password = 'root')
dbWriteTable(con.stad, name = "tcga_stad_fpkm_gene", value = STAD.FPKM.GENE, row.names = T, overwrite = T)
save(STAD.FPKM.GENE, file = "E:/DATA/TCGA-RDATA/tcga_stad_fpkm_gene.RData")

##################### 2016-07-28: tidy clinical files

# read clinical manifest
d.manifest = read.csv(file = "E:/DATA/TCGA-Data/STAD/clinical-manifest.tsv", sep = '\t',header = T)

# copy file to clinical direction
apply(d.manifest, 1, copyfiles <- function(aRow) {
  filepath = paste0("E:/DATA/TCGA-Data/STAD/Clinical1/",aRow[1], '/',aRow[2])
  file.copy(from = filepath, to = "E:/DATA/TCGA-Data/STAD/Clinical/")
  #file.remove(filepath)
} )

############ 2016-07-28: read and make clinical data.frame
library(xml2)
library(RMySQL)

# read clinical files
v.file.names = dir(path = "E:/DATA/TCGA-Data/STAD/Clinical", pattern = "*.xml", recursive = T, full.names = T)

# get a sample file to get clinical data.frame's colnames
root = as_list(read_xml(v.file.names[1]))
v.attr.name = names(root[[2]]) 
v.attr.value = vector()

d.stad.clinical = data.frame(matrix(NA, nrow = length(v.file.names), ncol = length(v.attr.name)), stringsAsFactors = F)
colnames(d.stad.clinical) = v.attr.name

# first loop: retrive clinical files
for (i.file in 1:length(v.file.names)) {
  
  # # first step: comment the second step, execute this step to find which file is invalid
  # tryCatch(expr = {root = as_list(read_xml(v.file.names[i.file]))},
  #          warning = function(war){
  #            print(war)
  #          },error = function(err){
  #            print(c(i.file, v.file.names[i.file]))
  #            next()
  #          },finally = {next()})
  
  # second step: retrive every attribute and add to data.frame
  root = as_list(read_xml(v.file.names[i.file]))
  v.attr.name = names(root[[2]])
  for( name in v.attr.name) {
    if(length(root[[2]][[name]]) == 1 && length(root[[2]][[name]][[1]]) == 1) {
      d.stad.clinical[i.file,name] = root[[2]][[name]][[1]]
    }else{
      d.stad.clinical[i.file,name] = "NULL"
    }
  }
}

#d.stad.clinical[i.file,]
#v.file.names[i.file]

############ 2016-07-22: read metadata files, map patients id to mRNA expression file uuid
library(jsonlite)
library(dplyr)

metadata = fromJSON("E:/DATA/TCGA-Data/STAD/mRNA-metadata.json")

d.stad.uuidmap = as.data.frame(t(data.frame(lapply(metadata$associated_entities, unlist))), stringsAsFactors = F)
d.stad.uuidmap = mutate(d.stad.uuidmap, mRNA_file_uuid = sapply(metadata$file_name,substr,1, 36,USE.NAMES = F))
rownames(d.stad.uuidmap) = NULL
## ***check****

load("E:/DATA/TCGA-RDATA/tcga_stad_fpkm_gene.RData")

## make the case same
d.stad.clinical$bcr_patient_uuid = tolower(d.stad.clinical$bcr_patient_uuid)
d.stad.uuidmap$case_id = tolower(d.stad.uuidmap$case_id)
colnames(STAD.FPKM.GENE) = tolower(colnames(STAD.FPKM.GENE))

## add barcode information to d.stad.uuidmap
d.stad.uuidmap = mutate(d.stad.uuidmap, 
                        patient_barcode = substr(d.stad.uuidmap$entity_submitter_id,1,12), 
                        tumor_sample = substr(d.stad.uuidmap$entity_submitter_id,14,15))

## split to 2 data.frame (Normal and Tumor)
d.stad.uuidmap.normal = filter(d.stad.uuidmap, tumor_sample == "11")
d.stad.uuidmap.tumor = filter(d.stad.uuidmap, tumor_sample == "01")

d.stad.fpkm.normal = STAD.FPKM.GENE[,d.stad.uuidmap.normal$mRNA_file_uuid]
d.stad.fpkm.tumor = STAD.FPKM.GENE[,d.stad.uuidmap.tumor$mRNA_file_uuid]

# change rowname (file uuid) by patient id
colnames(d.stad.fpkm.normal) = sapply(colnames(d.stad.fpkm.normal), function(a){filter(d.stad.uuidmap.normal, mRNA_file_uuid == a)$patient_barcode}, USE.NAMES = F)
colnames(d.stad.fpkm.tumor) = sapply(colnames(d.stad.fpkm.tumor), function(a){filter(d.stad.uuidmap.tumor, mRNA_file_uuid == a)$patient_barcode}, USE.NAMES = F)

## test if colnames have copies
nameids = colnames(d.stad.fpkm.tumor)
length(nameids) != length(unique(nameids))
if(length(nameids) != length(unique(nameids)))
{
  print("have copies!!!!")
  
  # get copies id
  copyids = unique(nameids[duplicated(nameids)])
  
  id.col = which(nameids %in% copyids)
  d.id = d.stad.fpkm.tumor[,id.col]
  d.stad.fpkm.tumor = d.stad.fpkm.tumor[,-id.col]
  
  #orde by colnames
  d.id = d.id[,order(colnames(d.id))]
  copyids = sort(copyids)
  # deal it one by one 
  print(copyids)
  print(colnames(d.id))
  
  # # get v.id by hand
  # v.id = list(c(1,2),c(3,4,5),c(6,7),c(8,9,10),c(11,12,13),c(14,15,16),c(17,18,19),c(20,21,22),c(23,24),c(25,26),c(27,28,29),c(30,31,32),c(33,34,35))
  # 
  # # combind each copies to data.frame
  # for(i in 1:13)
  # {
  #   an = paste0("a", i)
  #   d.stad.fpkm.tumor = cbind.data.frame(d.stad.fpkm.tumor, "an" = apply(d.id[,v.id[[i]]], 1, mean))
  #   colnames(d.stad.fpkm.tumor)[ which(colnames(d.stad.fpkm.tumor) == "an")] = copyids[i]
  # }
  
}

nameids = colnames(d.stad.fpkm.normal)
if(length(nameids) != length(unique(nameids)))
{
  print("have copies!!!!")
}

# reduce clinical data
# if the column:  NULL data > 90% data, then remove the 

# save data 
save(d.stad.clinical, file = 'E:/DATA/TCGA-RDATA/STAD/stad.clinical.RData')
save(d.stad.fpkm.normal, file = 'E:/DATA/TCGA-RDATA/STAD/stad.fpkm.normal.RData')
save(d.stad.fpkm.tumor, file = 'E:/DATA/TCGA-RDATA/STAD/stad.fpkm.tumor.RData')
save(d.stad.uuidmap, file = 'E:/DATA/TCGA-RDATA/STAD/stad.uuidmap.RData')

con.stad = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.stad', user = 'root', password = 'root')
dbWriteTable(con.stad, name = "tcga_stad_clinical", value = d.stad.clinical, row.names = F, overwrite = T)
dbWriteTable(con.stad, name = "tcga_stad_fpkm_tumor", value = d.stad.fpkm.tumor, row.names = T, overwrite = T)
dbWriteTable(con.stad, name = "tcga_stad_fpkm_normal", value = d.stad.fpkm.normal, row.names = T, overwrite = T)
dbWriteTable(con.stad, name = "tcga_stad_uuidmap", value = d.stad.uuidmap, row.names = F, overwrite = T)
