# ***************************************
# * Name: liangsen_2014@163.com
# * Date: Sun Jul 10 19:51:20 2016
# * Description: LIHC
# *              
# ***************************************

# define two file paths
datapath = 'E:/DATA/'

# read manifest file
d.manifest = read.csv(file = 'E:/DATA/TCGA-manifest-file-14/LIHC-gdc_manifest.2016-07-10T13-26-18.842797.tsv', sep = '\t',header = T)

# get a test files of XXX.FPKM.txt.gz
FPKM.gz.path = paste(datapath, 'TCGA-Data-14-20160704/',d.manifest$id[1],'/',d.manifest$filename[1], sep = '')
testdata = read.csv(file = FPKM.gz.path, sep = '\t', header = F)

# define a data.frame for all data 
LIHC.FPKM = data.frame(matrix(NA, nrow = nrow(testdata), ncol = nrow(d.manifest)), stringsAsFactors = F)
colnames(LIHC.FPKM) = apply(X = d.manifest, MARGIN = 1, FUN = function(aRow){substr(aRow[2], 1, 36)})
rownames(LIHC.FPKM) = testdata[,1]

# get all data form files
for(i.id in 1:nrow(d.manifest))
{
  
  FPKM.gz.path = paste(datapath, 'TCGA-Data-14-20160704/',d.manifest$id[i.id],'/',d.manifest$filename[i.id], sep = '')
  LIHC.FPKM[,i.id] = read.csv(file = FPKM.gz.path, sep = '\t', header = F)[,2]
  
  # remove file to another directions
  file.copy(from = FPKM.gz.path, to = 'E:/DATA/TCGA-Data/LIHC/',overwrite = F)
  print(i.id)
}

save(LIHC.FPKM, file = 'E:/DATA/TCGA-RDATA/LIHC.FPKM.RData')

########### 2016-07-15
load('E:/DATA/TCGA-RDATA/LIHC.FPKM.RData')

#save to mysql
con.lihc = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.lihc', user = 'root', password = 'root')

dbWriteTable(con.lihc, name = "tcga_lihc_fpkm", value = LIHC.FPKM, row.names = T, overwrite = F)

# get tcga gene list
#con.anno = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'annotation', user = 'root', password = 'root')
#tcga.gene.list = dbReadTable(conn = con.anno, name= "ensembl_id_tcga")[,1]

# get protein_coding list
#tcga.gene.map = dbReadTable(conn = con.anno, name= "tcga_gene_map")

# set name
rownames(LIHC.FPKM) = tcga.gene.list

# filter
LIHC.FPKM.GENE = LIHC.FPKM[tcga.gene.map$Ensembl.Gene.ID,] 
rownames(LIHC.FPKM.GENE) = tcga.gene.map$Gene.ID

# save to mysql and RDATA
#con.lihc = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.lihc', user = 'root', password = 'root')
dbWriteTable(con.lihc, name = "tcga_lihc_fpkm_gene", value = LIHC.FPKM.GENE, row.names = T, overwrite = T)
save(LIHC.FPKM.GENE, file = "E:/DATA/TCGA-RDATA/tcga_lihc_fpkm_gene.RData")

##################### 2016-07-28: tidy clinical files

# read clinical manifest
d.manifest = read.csv(file = "E:/DATA/TCGA-Data/LIHC/clinical-manifest.txt", sep = '\t',header = T)

# copy file to clinical direction
apply(d.manifest, 1, copyfiles <- function(aRow) {
  filepath = paste0("E:/DATA/TCGA-20160723/",aRow[1], '/',aRow[2])
  file.copy(from = filepath, to = "E:/DATA/TCGA-Data/LIHC/Clinical/")
  file.remove(filepath)
} )

############ 2016-07-28: read and make clinical data.frame
library(xml2)
library(RMySQL)

# read clinical files
v.file.names = dir(path = "E:/DATA/TCGA-Data/LIHC/Clinical", pattern = "*.xml", recursive = T, full.names = T)

# get a sample file to get clinical data.frame's colnames
root = as_list(read_xml(v.file.names[1]))
v.attr.name = names(root[[2]]) 
v.attr.value = vector()

d.lihc.clinical = data.frame(matrix(NA, nrow = length(v.file.names), ncol = length(v.attr.name)), stringsAsFactors = F)
colnames(d.lihc.clinical) = v.attr.name

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
      d.lihc.clinical[i.file,name] = root[[2]][[name]][[1]]
    }else{
      d.lihc.clinical[i.file,name] = "NULL"
    }
  }
}

#d.lihc.clinical[i.file,]
#v.file.names[i.file]

############ 2016-07-22: read metadata files, map patients id to mRNA expression file uuid
library(jsonlite)
library(dplyr)

metadata = fromJSON("E:/DATA/TCGA-Data/LIHC/mRNA-metadata.json")

d.lihc.uuidmap = as.data.frame(t(data.frame(lapply(metadata$associated_entities, unlist))), stringsAsFactors = F)
d.lihc.uuidmap = mutate(d.lihc.uuidmap, mRNA_file_uuid = sapply(metadata$file_name,substr,1, 36,USE.NAMES = F))
rownames(d.lihc.uuidmap) = NULL
## ***check****

load("E:/DATA/TCGA-RDATA/tcga_lihc_fpkm_gene.RData")

## make the case same
d.lihc.clinical$bcr_patient_uuid = tolower(d.lihc.clinical$bcr_patient_uuid)
d.lihc.uuidmap$case_id = tolower(d.lihc.uuidmap$case_id)
colnames(LIHC.FPKM.GENE) = tolower(colnames(LIHC.FPKM.GENE))

## add barcode information to d.lihc.uuidmap
d.lihc.uuidmap = mutate(d.lihc.uuidmap, 
                        patient_barcode = substr(d.lihc.uuidmap$entity_submitter_id,1,12), 
                        tumor_sample = substr(d.lihc.uuidmap$entity_submitter_id,14,15))

## split to 2 data.frame (Normal and Tumor)
d.lihc.uuidmap.normal = filter(d.lihc.uuidmap, tumor_sample == "11")
d.lihc.uuidmap.tumor = filter(d.lihc.uuidmap, tumor_sample == "01")

d.lihc.fpkm.normal = LIHC.FPKM.GENE[,d.lihc.uuidmap.normal$mRNA_file_uuid]
d.lihc.fpkm.tumor = LIHC.FPKM.GENE[,d.lihc.uuidmap.tumor$mRNA_file_uuid]

# change rowname (file uuid) by patient id
colnames(d.lihc.fpkm.normal) = sapply(colnames(d.lihc.fpkm.normal), function(a){filter(d.lihc.uuidmap.normal, mRNA_file_uuid == a)$patient_barcode}, USE.NAMES = F)
colnames(d.lihc.fpkm.tumor) = sapply(colnames(d.lihc.fpkm.tumor), function(a){filter(d.lihc.uuidmap.tumor, mRNA_file_uuid == a)$patient_barcode}, USE.NAMES = F)

## test if colnames have copies
nameids = colnames(d.lihc.fpkm.tumor)
length(nameids) != length(unique(nameids))
if(length(nameids) != length(unique(nameids)))
{
  print("have copies!!!!")
  
  # get copies id
  copyids = unique(nameids[duplicated(nameids)])
  
  id.col = which(nameids %in% copyids)
  d.id = d.lihc.fpkm.tumor[,id.col]
  d.lihc.fpkm.tumor = d.lihc.fpkm.tumor[,-id.col]
  
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
  #   d.lihc.fpkm.tumor = cbind.data.frame(d.lihc.fpkm.tumor, "an" = apply(d.id[,v.id[[i]]], 1, mean))
  #   colnames(d.lihc.fpkm.tumor)[ which(colnames(d.lihc.fpkm.tumor) == "an")] = copyids[i]
  # }
  
}

nameids = colnames(d.lihc.fpkm.normal)
if(length(nameids) != length(unique(nameids)))
{
  print("have copies!!!!")
}

# reduce clinical data
# if the column:  NULL data > 90% data, then remove the 

# save data 
save(d.lihc.clinical, file = 'E:/DATA/TCGA-RDATA/LIHC/lihc.clinical.RData')
save(d.lihc.fpkm.normal, file = 'E:/DATA/TCGA-RDATA/LIHC/lihc.fpkm.normal.RData')
save(d.lihc.fpkm.tumor, file = 'E:/DATA/TCGA-RDATA/LIHC/lihc.fpkm.tumor.RData')
save(d.lihc.uuidmap, file = 'E:/DATA/TCGA-RDATA/LIHC/lihc.uuidmap.RData')

con.lihc = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.lihc', user = 'root', password = 'root')
dbWriteTable(con.lihc, name = "tcga_lihc_clinical", value = d.lihc.clinical, row.names = F, overwrite = T)
dbWriteTable(con.lihc, name = "tcga_lihc_fpkm_tumor", value = d.lihc.fpkm.tumor, row.names = T, overwrite = T)
dbWriteTable(con.lihc, name = "tcga_lihc_fpkm_normal", value = d.lihc.fpkm.normal, row.names = T, overwrite = T)
dbWriteTable(con.lihc, name = "tcga_lihc_uuidmap", value = d.lihc.uuidmap, row.names = F, overwrite = T)
