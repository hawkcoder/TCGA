# ***************************************
# * Name: liangsen_2014@163.com
# * Date: Sun Jul 10 19:51:20 2016
# * Description: BLCA
# *              
# ***************************************

# read manifest file
d.manifest = read.csv(file = 'E:/DATA/TCGA-manifest-file-14/BLCA-gdc_manifest.2016-07-10T13-25-42.227541.tsv', sep = '\t',header = T)

# get a test files of XXX.FPKM.txt.gz
FPKM.gz.path = paste(datapath, 'TCGA-Data-14-20160704/',d.manifest$id[1],'/',d.manifest$filename[1], sep = '')
testdata = read.csv(file = FPKM.gz.path, sep = '\t', header = F)

# define a data.frame for all data 
BLCA.FPKM = data.frame(matrix(NA, nrow = nrow(testdata), ncol = nrow(d.manifest)), stringsAsFactors = F)
colnames(BLCA.FPKM) = apply(X = d.manifest, MARGIN = 1, FUN = function(aRow){substr(aRow[2], 1, 36)})
rownames(BLCA.FPKM) = testdata[,1]

# get all data form files
for(i.id in 1:nrow(d.manifest))
{
  
  FPKM.gz.path = paste(datapath, 'TCGA-Data-14-20160704/',d.manifest$id[i.id],'/',d.manifest$filename[i.id], sep = '')
  BLCA.FPKM[,i.id] = read.csv(file = FPKM.gz.path, sep = '\t', header = F)[,2]
  
  # remove file to another directions
  file.copy(from = FPKM.gz.path, to = 'E:/DATA/TCGA-Data/BLCA/',overwrite = F)
  print(i.id)
}

save(BLCA.FPKM, file = 'E:/DATA/TCGA-RDATA/BLCA.FPKM.RData')




########### 2016-07-15
load('E:/DATA/TCGA-RDATA/BLCA.FPKM.RData')

#save to mysql
con.blca = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.blca', user = 'root', password = 'root')
dbWriteTable(con.blca, name = "tcga_blca_fpkm", value = BLCA.FPKM, row.names = T, overwrite = T)

# get tcga gene list
con.anno = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'annotation', user = 'root', password = 'root')
tcga.gene.list = dbReadTable(conn = con.anno, name= "ensembl_id_tcga")[,1]

# get protein_coding list
tcga.gene.map = dbReadTable(conn = con.anno, name= "tcga_gene_map")

# set name
rownames(BLCA.FPKM) = tcga.gene.list

# filter
BLCA.FPKM.GENE = BLCA.FPKM[tcga.gene.map$Ensembl.Gene.ID,] 
rownames(BLCA.FPKM.GENE) = tcga.gene.map$Gene.ID

# save to mysql and RDATA
con.blca = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.blca', user = 'root', password = 'root')
dbWriteTable(con.blca, name = "tcga_blca_fpkm_gene", value = BLCA.FPKM.GENE, row.names = T, overwrite = T)
save(BLCA.FPKM.GENE, file = "E:/DATA/TCGA-RDATA/tcga_blca_fpkm_gene.RData")

##################### 2016-07-23: tidy clinical files

# read clinical manifest
d.manifest = read.csv(file = "E:/DATA/TCGA-Data/BLCA/clinical-manifest.txt", sep = '\t',header = T)

# copy file to clinical direction
apply(d.manifest, 1, copyfiles <- function(aRow) {
  filepath = paste0("E:/DATA/TCGA-20160723/",aRow[1], '/',aRow[2])
  file.copy(from = filepath, to = "E:/DATA/TCGA-Data/BLCA/Clinical/")
  file.remove(filepath)
} )


############ 2016-07-22: read and make clinical data.frame
# library(xml2)
# library(RMySQL)
# 
# # read clinical files
# v.file.names = dir(path = "E:/DATA/TCGA-Data/BLCA/Clinical", pattern = "*.xml", recursive = T, full.names = T)
# 
# # get a sample file to get clinical data.frame's colnames
# root = as_list(read_xml(v.file.names[1]))
# v.attr.name = names(root[[2]]) 
# v.attr.value = vector()
# 
# d.blca.clinical = data.frame(matrix(NA, nrow = length(v.file.names), ncol = length(v.attr.name)), stringsAsFactors = F)
# colnames(d.blca.clinical) = v.attr.name
# 
# # first loop: retrive clinical files
# for (i.file in 1:length(v.file.names)) {
#   root = as_list(read_xml(v.file.names[i.file]))
#   v.attr.name = names(root[[2]])
#   
#   # second loop: retrive every attribute and add to data.frame
#   for( name in v.attr.name) {
#     if(length(root[[2]][[name]]) == 1 && length(root[[2]][[name]][[1]]) == 1) {
#       d.blca.clinical[i.file,name] = root[[2]][[name]][[1]]
#     }else{
#       d.blca.clinical[i.file,name] = "NULL"
#     }
#   }
# }

############################ 2016-07-29: read and make clinical data.frame
library(xml2)
library(RMySQL)

# read clinical files
v.file.names = dir(path = "E:/DATA/TCGA-Data/BLCA/Clinical", pattern = "*.xml", recursive = T, full.names = T)

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
d.blca.clinical = data.frame(matrix(NA, nrow = length(v.file.names), ncol = length(v.clinical.name)), stringsAsFactors = F)
colnames(d.blca.clinical) = v.clinical.name


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
  d.blca.clinical[i.file,] = v.xml.value[v.clinical.name]
}

# save file
save(d.blca.clinical, file = 'E:/DATA/TCGA-RDATA/BLCA/blca.clinical.RData')

con.blca = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.blca', user = 'root', password = 'root')
dbWriteTable(con.blca, name = "tcga_blca_clinical", value = d.blca.clinical, row.names = F, overwrite = T)


############ 2016-07-22: read metadata files, map patients id to mRNA expression file uuid
library(jsonlite)
library(dplyr)

metadata = fromJSON("E:/DATA/TCGA-Data/BLCA/mRNA-metadata.json")

d.blca.uuidmap = as.data.frame(t(data.frame(lapply(metadata$associated_entities, unlist))), stringsAsFactors = F)
d.blca.uuidmap = mutate(d.blca.uuidmap, mRNA_file_uuid = sapply(metadata$file_name,substr,1, 36,USE.NAMES = F))
rownames(d.blca.uuidmap) = NULL
## ***check****

load("E:/DATA/TCGA-RDATA/tcga_blca_fpkm_gene.RData")

## make the case same
d.blca.clinical$bcr_patient_uuid = tolower(d.blca.clinical$bcr_patient_uuid)
d.blca.uuidmap$case_id = tolower(d.blca.uuidmap$case_id)
colnames(BLCA.FPKM.GENE) = tolower(colnames(BLCA.FPKM.GENE))

## add barcode information to d.blca.uuidmap
d.blca.uuidmap = mutate(d.blca.uuidmap, 
                   patient_barcode = substr(d.blca.uuidmap$entity_submitter_id,1,12), 
                   tumor_sample = substr(d.blca.uuidmap$entity_submitter_id,14,15))

## split to 2 data.frame (Normal and Tumor)
d.blca.uuidmap.normal = filter(d.blca.uuidmap, tumor_sample == "11")
d.blca.uuidmap.tumor = filter(d.blca.uuidmap, tumor_sample == "01")

d.blca.fpkm.normal = BLCA.FPKM.GENE[,d.blca.uuidmap.normal$mRNA_file_uuid]
d.blca.fpkm.tumor = BLCA.FPKM.GENE[,d.blca.uuidmap.tumor$mRNA_file_uuid]

# change rowname (file uuid) by patient id
colnames(d.blca.fpkm.normal) = sapply(colnames(d.blca.fpkm.normal), function(a){filter(d.blca.uuidmap.normal, mRNA_file_uuid == a)$patient_barcode}, USE.NAMES = F)
colnames(d.blca.fpkm.tumor) = sapply(colnames(d.blca.fpkm.tumor), function(a){filter(d.blca.uuidmap.tumor, mRNA_file_uuid == a)$patient_barcode}, USE.NAMES = F)

## test if colnames have copies
nameids = colnames(d.blca.fpkm.tumor)

if(length(nameids) != length(unique(nameids)))
{
  print("have copies!!!!")
  
  # get copies id
  copyids = unique(nameids[duplicated(nameids)])
  
  id.col = which(nameids %in% copyids)
  d.id = d.blca.fpkm.tumor[,id.col]
  d.blca.fpkm.tumor = d.blca.fpkm.tumor[,-id.col]
  
  # deal it one by one 
  id = copyids[1]
  aa = d.id[,c(2,3,6)] # "TCGA-BL-A0C8"
  bb = d.id[,c(1,5,9)] # "TCGA-BL-A13J"
  cc = d.id[,c(4,7,8)] # "TCGA-BL-A13I"
  
  d.blca.fpkm.tumor = cbind.data.frame(d.blca.fpkm.tumor, "TCGA-BL-A0C8" = apply(aa, 1, mean))
  d.blca.fpkm.tumor = cbind.data.frame(d.blca.fpkm.tumor, "TCGA-BL-A13J" = apply(bb, 1, mean))
  d.blca.fpkm.tumor = cbind.data.frame(d.blca.fpkm.tumor, "TCGA-BL-A13I" = apply(cc, 1, mean))
}

nameids = colnames(d.blca.fpkm.normal)
if(length(nameids) != length(unique(nameids)))
{
  print("have copies!!!!")
}


# save data 
save(d.blca.clinical, file = 'E:/DATA/TCGA-RDATA/BLCA/blca.clinical.RData')
save(d.blca.fpkm.normal, file = 'E:/DATA/TCGA-RDATA/BLCA/blca.fpkm.normal.RData')
save(d.blca.fpkm.tumor, file = 'E:/DATA/TCGA-RDATA/BLCA/blca.fpkm.tumor.RData')
save(d.blca.uuidmap, file = 'E:/DATA/TCGA-RDATA/BLCA/blca.uuidmap.RData')

save(d.blca.clinical, file = 'E:/DATA/TCGA-RDATA/BLCA/blca.clinical.RData')
con.blca = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.blca', user = 'root', password = 'root')
dbWriteTable(con.blca, name = "tcga_blca_clinical", value = d.blca.clinical, row.names = F, overwrite = T)
dbWriteTable(con.blca, name = "tcga_blca_fpkm_normal", value = d.blca.fpkm.normal, row.names = T, overwrite = T)
dbWriteTable(con.blca, name = "tcga_blca_fpkm_tumor", value = d.blca.fpkm.tumor, row.names = T, overwrite = T)
dbWriteTable(con.blca, name = "tcga_blca_uuidmap", value = d.blca.uuidmap, row.names = F, overwrite = T)




























