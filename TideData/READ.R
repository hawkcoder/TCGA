# ***************************************
# * Name: liangsen_2014@163.com
# * Date: Sun Jul 10 19:51:20 2016
# * Description: READ
# *              
# ***************************************

# define two file paths
datapath = 'E:/DATA/'

# read manifest file
d.manifest = read.csv(file = 'E:/DATA/TCGA-Data/READ/mRNA-manifest.tsv', sep = '\t',header = T)

# copy file to mRNA-FPKM direction
apply(d.manifest, 1, copyfiles <- function(aRow) {
  filepath = paste0("E:/DATA/TCGA-Data/READ/mRNA-FPKM1/",aRow[1], '/',aRow[2])
  file.copy(from = filepath, to = "E:/DATA/TCGA-Data/READ/mRNA-FPKM/")
  #file.remove(filepath)
} )


# get a test files of XXX.FPKM.txt.gz
FPKM.gz.path = paste0('E:/DATA/TCGA-Data/READ/mRNA-FPKM/',d.manifest$filename[1])
testdata = read.csv(file = FPKM.gz.path, sep = '\t', header = F)

# define a data.frame for all data 
d.read.fpkm.all.gene = data.frame(matrix(NA, nrow = nrow(testdata), ncol = nrow(d.manifest)), stringsAsFactors = F)
colnames(d.read.fpkm.all.gene) = apply(X = d.manifest, MARGIN = 1, FUN = function(aRow){substr(aRow[2], 1, 36)})
rownames(d.read.fpkm.all.gene) = testdata[,1]

# get all data form files
for(i.id in 1:nrow(d.manifest))
{
  
  FPKM.gz.path = paste0('E:/DATA/TCGA-Data/READ/mRNA-FPKM/',d.manifest$filename[i.id])
  d.read.fpkm.all.gene[,i.id] = read.csv(file = FPKM.gz.path, sep = '\t', header = F)[,2]
  print(i.id)
}

save(d.read.fpkm.all.gene, file = 'E:/DATA/TCGA-RDATA/READ/read.fpkm.all.gene.RData')

########### 2016-07-15
#load('E:/DATA/TCGA-RDATA/READ/read.fpkm.all.gene.RData')

#save to mysql
con.read = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.read', user = 'root', password = 'root')
dbWriteTable(con.read, name = "tcga_read_fpkm", value = d.read.fpkm.all.gene, row.names = T, overwrite = F)

# get tcga gene list
con.anno = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'annotation', user = 'root', password = 'root')
tcga.gene.list = dbReadTable(conn = con.anno, name= "ensembl_id_tcga")[,1]

# get protein_coding list
tcga.gene.map = dbReadTable(conn = con.anno, name= "tcga_gene_map")

# set name
rownames(d.read.fpkm.all.gene) = tcga.gene.list

# filter
d.read.fpkm = d.read.fpkm.all.gene[tcga.gene.map$Ensembl.Gene.ID,] 
rownames(d.read.fpkm) = tcga.gene.map$Gene.ID

# save to mysql and RDATA
#con.read = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.read', user = 'root', password = 'root')
dbWriteTable(con.read, name = "tcga_read_fpkm_gene", value = d.read.fpkm, row.names = T, overwrite = T)
save(d.read.fpkm, file = "E:/DATA/TCGA-RDATA/READ/read.fpkm.RData")

##################### 2016-07-28: tidy clinical files

# read clinical manifest
d.manifest = read.csv(file = "E:/DATA/TCGA-Data/READ/clinical-manifest.tsv", sep = '\t',header = T)

# copy file to clinical direction
apply(d.manifest, 1, copyfiles <- function(aRow) {
  filepath = paste0("E:/DATA/TCGA-Data/READ/Clinical1/",aRow[1], '/',aRow[2])
  file.copy(from = filepath, to = "E:/DATA/TCGA-Data/READ/Clinical/")
  #file.remove(filepath)
} )

# ############ 2016-07-28: read and make clinical data.frame
# library(xml2)
# library(RMySQL)
# 
# # read clinical files
# v.file.names = dir(path = "E:/DATA/TCGA-Data/READ/Clinical", pattern = "*.xml", recursive = T, full.names = T)
# 
# # get a sample file to get clinical data.frame's colnames
# root = as_list(read_xml(v.file.names[1]))
# v.attr.name = names(root[[2]]) 
# v.attr.value = vector()
# 
# d.read.clinical = data.frame(matrix(NA, nrow = length(v.file.names), ncol = length(v.attr.name)), stringsAsFactors = F)
# colnames(d.read.clinical) = v.attr.name
# 
# # first loop: retrive clinical files
# for (i.file in 1:length(v.file.names)) {
#   
#   # # first step: comment the second step, execute this step to find which file is invalid
#   # tryCatch(expr = {root = as_list(read_xml(v.file.names[i.file]))},
#   #          warning = function(war){
#   #            print(war)
#   #          },error = function(err){
#   #            print(c(i.file, v.file.names[i.file]))
#   #            next()
#   #          },finally = {next()})
#   
#   # second step: retrive every attribute and add to data.frame
#   root = as_list(read_xml(v.file.names[i.file]))
#   v.attr.name = names(root[[2]])
#   for( name in v.attr.name) {
#     if(length(root[[2]][[name]]) == 1 && length(root[[2]][[name]][[1]]) == 1) {
#       d.read.clinical[i.file,name] = root[[2]][[name]][[1]]
#     }else{
#       d.read.clinical[i.file,name] = "NULL"
#     }
#   }
# }
# 
# #d.read.clinical[i.file,]
# #v.file.names[i.file]
############################ 2016-07-29: read and make clinical data.frame
library(xml2)
library(RMySQL)

# read clinical files
v.file.names = dir(path = "E:/DATA/TCGA-Data/READ/Clinical", pattern = "*.xml", recursive = T, full.names = T)

## get a sample file to get clinical data.frame's colnames
xml.nodes = xml_find_all(read_xml(v.file.names[1]), "//*")
v.name = vector()
for (lst in xml.nodes) {
  # if it is a leaf element, then leaves
  if(length(xml_find_all(lst, ".//*")) == 0)
  {
    v.name = c(v.name, xml_name(lst))
  }
}
# delete name which have copies
copies = unique(v.name[duplicated(v.name)])
v.name = v.name[which(!v.name %in% copies)]


v.clinical.name = v.name
# get common attribute names
for (file in v.file.names) {
  xml.nodes = xml_find_all(read_xml(file), "//*")
  v.xml.name = xml_name(xml.nodes)
  copies = unique(v.xml.name[duplicated(v.xml.name)])
  
  v.xml.name = v.xml.name[which(!v.xml.name %in% copies)]
  v.clinical.name = intersect(v.clinical.name, v.xml.name)
}

## start bulid a data.frame
d.read.clinical = data.frame(matrix(NA, nrow = length(v.file.names), ncol = length(v.clinical.name)), stringsAsFactors = F)
colnames(d.read.clinical) = v.clinical.name


# first loop: retrive clinical files
for (i.file in 1:length(v.file.names)) {
  
  # first step: comment the second step, execute this step to find which file is invalid
  # 
  # second step: retrive every attribute and add to data.frame
  xml.nodes = xml_find_all(read_xml(v.file.names[i.file]), "//*")
  v.xml.name = xml_name(xml.nodes)
  v.xml.value = xml_text(xml.nodes)
  
  order1 =  which(v.xml.name %in% v.clinical.name)
  v.xml.name = v.xml.name[order1]
  v.xml.value = v.xml.value[order1]
  names(v.xml.value) = v.xml.name
  
  # add to data.frame
  d.read.clinical[i.file,] = v.xml.value[v.clinical.name]
}

# save file
save(d.read.clinical, file = 'E:/DATA/TCGA-RDATA/READ/read.clinical.RData')

con.read = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.read', user = 'root', password = 'root')
dbWriteTable(con.read, name = "tcga_read_clinical", value = d.read.clinical, row.names = F, overwrite = T)

############ 2016-07-22: read metadata files, map patients id to mRNA expression file uuid
library(jsonlite)
library(dplyr)

metadata = fromJSON("E:/DATA/TCGA-Data/READ/mRNA-metadata.json")

d.read.uuidmap = as.data.frame(t(data.frame(lapply(metadata$associated_entities, unlist))), stringsAsFactors = F)
d.read.uuidmap = mutate(d.read.uuidmap, mRNA_file_uuid = sapply(metadata$file_name,substr,1, 36,USE.NAMES = F))
rownames(d.read.uuidmap) = NULL
## ***check****

#load("E:/DATA/TCGA-RDATA/READ/tcga_read_fpkm_gene.RData")

## make the case same
d.read.clinical$bcr_patient_uuid = tolower(d.read.clinical$bcr_patient_uuid)
d.read.uuidmap$case_id = tolower(d.read.uuidmap$case_id)
colnames(d.read.fpkm) = tolower(colnames(d.read.fpkm))

## add barcode information to d.read.uuidmap
d.read.uuidmap = mutate(d.read.uuidmap, 
                        patient_barcode = substr(d.read.uuidmap$entity_submitter_id,1,12), 
                        tumor_sample = substr(d.read.uuidmap$entity_submitter_id,14,15))

## split to 2 data.frame (Normal and Tumor)
d.read.uuidmap.normal = filter(d.read.uuidmap, tumor_sample == "11")
d.read.uuidmap.tumor = filter(d.read.uuidmap, tumor_sample == "01")

d.read.fpkm.normal = d.read.fpkm[,d.read.uuidmap.normal$mRNA_file_uuid]
d.read.fpkm.tumor = d.read.fpkm[,d.read.uuidmap.tumor$mRNA_file_uuid]

# change rowname (file uuid) by patient id
colnames(d.read.fpkm.normal) = sapply(colnames(d.read.fpkm.normal), function(a){filter(d.read.uuidmap.normal, mRNA_file_uuid == a)$patient_barcode}, USE.NAMES = F)
colnames(d.read.fpkm.tumor) = sapply(colnames(d.read.fpkm.tumor), function(a){filter(d.read.uuidmap.tumor, mRNA_file_uuid == a)$patient_barcode}, USE.NAMES = F)

## test if colnames have copies
nameids = colnames(d.read.fpkm.tumor)
length(nameids) != length(unique(nameids))
if(length(nameids) != length(unique(nameids)))
{
  print("have copies!!!!")
  
  # get copies id
  copyids = unique(nameids[duplicated(nameids)])
  
  id.col = which(nameids %in% copyids)
  d.id = d.read.fpkm.tumor[,id.col]
  d.read.fpkm.tumor = d.read.fpkm.tumor[,-id.col]
  
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
  #   d.read.fpkm.tumor = cbind.data.frame(d.read.fpkm.tumor, "an" = apply(d.id[,v.id[[i]]], 1, mean))
  #   colnames(d.read.fpkm.tumor)[ which(colnames(d.read.fpkm.tumor) == "an")] = copyids[i]
  # }
  
}

nameids = colnames(d.read.fpkm.normal)
if(length(nameids) != length(unique(nameids)))
{
  print("have copies!!!!")
}

# reduce clinical data
# if the column:  NULL data > 90% data, then remove the 

# save data 
save(d.read.clinical, file = 'E:/DATA/TCGA-RDATA/READ/read.clinical.RData')
save(d.read.fpkm.normal, file = 'E:/DATA/TCGA-RDATA/READ/read.fpkm.normal.RData')
save(d.read.fpkm.tumor, file = 'E:/DATA/TCGA-RDATA/READ/read.fpkm.tumor.RData')
save(d.read.uuidmap, file = 'E:/DATA/TCGA-RDATA/READ/read.uuidmap.RData')

con.read = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.read', user = 'root', password = 'root')
dbWriteTable(con.read, name = "tcga_read_clinical", value = d.read.clinical, row.names = F, overwrite = T)
dbWriteTable(con.read, name = "tcga_read_fpkm_tumor", value = d.read.fpkm.tumor, row.names = T, overwrite = T)
dbWriteTable(con.read, name = "tcga_read_fpkm_normal", value = d.read.fpkm.normal, row.names = T, overwrite = T)
dbWriteTable(con.read, name = "tcga_read_uuidmap", value = d.read.uuidmap, row.names = F, overwrite = T)

