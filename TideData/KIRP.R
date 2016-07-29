# ***************************************
# * Name: liangsen_2014@163.com
# * Date: Sun Jul 10 19:51:20 2016
# * Description: KIRP
# *              
# ***************************************

# define two file paths
datapath = 'E:/DATA/'

# read manifest file
d.manifest = read.csv(file = 'E:/DATA/TCGA-manifest-file-14/KIRP-gdc_manifest.2016-07-10T13-26-49.194029.tsv', sep = '\t',header = T)

# get a test files of XXX.FPKM.txt.gz
FPKM.gz.path = paste(datapath, 'TCGA-Data-14-20160704/',d.manifest$id[1],'/',d.manifest$filename[1], sep = '')
testdata = read.csv(file = FPKM.gz.path, sep = '\t', header = F)

# define a data.frame for all data 
KIRP.FPKM = data.frame(matrix(NA, nrow = nrow(testdata), ncol = nrow(d.manifest)), stringsAsFactors = F)
colnames(KIRP.FPKM) = apply(X = d.manifest, MARGIN = 1, FUN = function(aRow){substr(aRow[2], 1, 36)})
rownames(KIRP.FPKM) = testdata[,1]

# get all data form files
for(i.id in 1:nrow(d.manifest))
{
  
  FPKM.gz.path = paste(datapath, 'TCGA-Data-14-20160704/',d.manifest$id[i.id],'/',d.manifest$filename[i.id], sep = '')
  KIRP.FPKM[,i.id] = read.csv(file = FPKM.gz.path, sep = '\t', header = F)[,2]
  
  # remove file to another directions
  file.copy(from = FPKM.gz.path, to = 'E:/DATA/TCGA-Data/KIRP/',overwrite = F)
  print(i.id)
}

save(KIRP.FPKM, file = 'E:/DATA/TCGA-RDATA/KIRP.FPKM.RData')

########### 2016-07-15
load('E:/DATA/TCGA-RDATA/KIRP.FPKM.RData')

#save to mysql
con.kirp = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.kirp', user = 'root', password = 'root')

dbWriteTable(con.kirp, name = "tcga_kirp_fpkm", value = KIRP.FPKM, row.names = T, overwrite = F)

# get tcga gene list
#con.anno = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'annotation', user = 'root', password = 'root')
#tcga.gene.list = dbReadTable(conn = con.anno, name= "ensembl_id_tcga")[,1]

# get protein_coding list
#tcga.gene.map = dbReadTable(conn = con.anno, name= "tcga_gene_map")

# set name
rownames(KIRP.FPKM) = tcga.gene.list

# filter
KIRP.FPKM.GENE = KIRP.FPKM[tcga.gene.map$Ensembl.Gene.ID,] 
rownames(KIRP.FPKM.GENE) = tcga.gene.map$Gene.ID

# save to mysql and RDATA
#con.kirp = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.kirp', user = 'root', password = 'root')
dbWriteTable(con.kirp, name = "tcga_kirp_fpkm_gene", value = KIRP.FPKM.GENE, row.names = T, overwrite = T)
save(KIRP.FPKM.GENE, file = "E:/DATA/TCGA-RDATA/tcga_kirp_fpkm_gene.RData")

##################### 2016-07-28: tidy clinical files

# read clinical manifest
d.manifest = read.csv(file = "E:/DATA/TCGA-Data/KIRP/clinical-manifest.txt", sep = '\t',header = T)

# copy file to clinical direction
apply(d.manifest, 1, copyfiles <- function(aRow) {
  filepath = paste0("E:/DATA/TCGA-20160723/",aRow[1], '/',aRow[2])
  file.copy(from = filepath, to = "E:/DATA/TCGA-Data/KIRP/Clinical/")
  file.remove(filepath)
} )

# ############ 2016-07-28: read and make clinical data.frame
# library(xml2)
# library(RMySQL)
# 
# # read clinical files
# v.file.names = dir(path = "E:/DATA/TCGA-Data/KIRP/Clinical", pattern = "*.xml", recursive = T, full.names = T)
# 
# # get a sample file to get clinical data.frame's colnames
# root = as_list(read_xml(v.file.names[1]))
# v.attr.name = names(root[[2]]) 
# v.attr.value = vector()
# 
# d.kirp.clinical = data.frame(matrix(NA, nrow = length(v.file.names), ncol = length(v.attr.name)), stringsAsFactors = F)
# colnames(d.kirp.clinical) = v.attr.name
# 
# # first loop: retrive clinical files
# for (i.file in 1:length(v.file.names)) {
#   root = as_list(read_xml(v.file.names[i.file]))
#   # tryCatch(expr = {root = as_list(read_xml(v.file.names[i.file]))},
#   #          warning = function(war){
#   #            print(war)
#   #          },error = function(err){
#   #            print(v.file.names[i.file])
#   #            next()
#   #          },finally = {next()})
#   
#   v.attr.name = names(root[[2]])
#   
#   # second loop: retrive every attribute and add to data.frame
#   for( name in v.attr.name) {
#     if(length(root[[2]][[name]]) == 1 && length(root[[2]][[name]][[1]]) == 1) {
#       d.kirp.clinical[i.file,name] = root[[2]][[name]][[1]]
#     }else{
#       d.kirp.clinical[i.file,name] = "NULL"
#     }
#   }
# }
# 
# #d.kirp.clinical[i.file,]
# #v.file.names[i.file]

############################ 2016-07-29: read and make clinical data.frame
library(xml2)
library(RMySQL)

# read clinical files
v.file.names = dir(path = "E:/DATA/TCGA-Data/KIRP/Clinical", pattern = "*.xml", recursive = T, full.names = T)

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
d.kirp.clinical = data.frame(matrix(NA, nrow = length(v.file.names), ncol = length(v.clinical.name)), stringsAsFactors = F)
colnames(d.kirp.clinical) = v.clinical.name


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
  d.kirp.clinical[i.file,] = v.xml.value[v.clinical.name]
}

# save file
save(d.kirp.clinical, file = 'E:/DATA/TCGA-RDATA/KIRP/kirp.clinical.RData')

con.kirp = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.kirp', user = 'root', password = 'root')
dbWriteTable(con.kirp, name = "tcga_kirp_clinical", value = d.kirp.clinical, row.names = F, overwrite = T)


############ 2016-07-22: read metadata files, map patients id to mRNA expression file uuid
library(jsonlite)
library(dplyr)

metadata = fromJSON("E:/DATA/TCGA-Data/KIRP/mRNA-metadata.json")

d.kirp.uuidmap = as.data.frame(t(data.frame(lapply(metadata$associated_entities, unlist))), stringsAsFactors = F)
d.kirp.uuidmap = mutate(d.kirp.uuidmap, mRNA_file_uuid = sapply(metadata$file_name,substr,1, 36,USE.NAMES = F))
rownames(d.kirp.uuidmap) = NULL
## ***check****

load("E:/DATA/TCGA-RDATA/tcga_kirp_fpkm_gene.RData")

## make the case same
d.kirp.clinical$bcr_patient_uuid = tolower(d.kirp.clinical$bcr_patient_uuid)
d.kirp.uuidmap$case_id = tolower(d.kirp.uuidmap$case_id)
colnames(KIRP.FPKM.GENE) = tolower(colnames(KIRP.FPKM.GENE))

## add barcode information to d.kirp.uuidmap
d.kirp.uuidmap = mutate(d.kirp.uuidmap, 
                        patient_barcode = substr(d.kirp.uuidmap$entity_submitter_id,1,12), 
                        tumor_sample = substr(d.kirp.uuidmap$entity_submitter_id,14,15))

## split to 2 data.frame (Normal and Tumor)
d.kirp.uuidmap.normal = filter(d.kirp.uuidmap, tumor_sample == "11")
d.kirp.uuidmap.tumor = filter(d.kirp.uuidmap, tumor_sample == "01")

d.kirp.fpkm.normal = KIRP.FPKM.GENE[,d.kirp.uuidmap.normal$mRNA_file_uuid]
d.kirp.fpkm.tumor = KIRP.FPKM.GENE[,d.kirp.uuidmap.tumor$mRNA_file_uuid]

# change rowname (file uuid) by patient id
colnames(d.kirp.fpkm.normal) = sapply(colnames(d.kirp.fpkm.normal), function(a){filter(d.kirp.uuidmap.normal, mRNA_file_uuid == a)$patient_barcode}, USE.NAMES = F)
colnames(d.kirp.fpkm.tumor) = sapply(colnames(d.kirp.fpkm.tumor), function(a){filter(d.kirp.uuidmap.tumor, mRNA_file_uuid == a)$patient_barcode}, USE.NAMES = F)

## test if colnames have copies
nameids = colnames(d.kirp.fpkm.tumor)
length(nameids) != length(unique(nameids))
if(length(nameids) != length(unique(nameids)))
{
  print("have copies!!!!")
  
  # get copies id
  copyids = unique(nameids[duplicated(nameids)])
  
  id.col = which(nameids %in% copyids)
  d.id = d.kirp.fpkm.tumor[,id.col]
  d.kirp.fpkm.tumor = d.kirp.fpkm.tumor[,-id.col]
  
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
  #   d.kirp.fpkm.tumor = cbind.data.frame(d.kirp.fpkm.tumor, "an" = apply(d.id[,v.id[[i]]], 1, mean))
  #   colnames(d.kirp.fpkm.tumor)[ which(colnames(d.kirp.fpkm.tumor) == "an")] = copyids[i]
  # }
  
}

nameids = colnames(d.kirp.fpkm.normal)
if(length(nameids) != length(unique(nameids)))
{
  print("have copies!!!!")
}

# reduce clinical data
# if the column:  NULL data > 90% data, then remove the 

# save data 
save(d.kirp.clinical, file = 'E:/DATA/TCGA-RDATA/KIRP/kirp.clinical.RData')
save(d.kirp.fpkm.normal, file = 'E:/DATA/TCGA-RDATA/KIRP/kirp.fpkm.normal.RData')
save(d.kirp.fpkm.tumor, file = 'E:/DATA/TCGA-RDATA/KIRP/kirp.fpkm.tumor.RData')
save(d.kirp.uuidmap, file = 'E:/DATA/TCGA-RDATA/KIRP/kirp.uuidmap.RData')

con.kirp = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.kirp', user = 'root', password = 'root')
dbWriteTable(con.kirp, name = "tcga_kirp_clinical", value = d.kirp.clinical, row.names = F, overwrite = T)
dbWriteTable(con.kirp, name = "tcga_kirp_fpkm_tumor", value = d.kirp.fpkm.tumor, row.names = T, overwrite = T)
dbWriteTable(con.kirp, name = "tcga_kirp_fpkm_normal", value = d.kirp.fpkm.normal, row.names = T, overwrite = T)
dbWriteTable(con.kirp, name = "tcga_kirp_uuidmap", value = d.kirp.uuidmap, row.names = F, overwrite = T)
