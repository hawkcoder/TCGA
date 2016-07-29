# ***************************************
# * Name: liangsen_2014@163.com
# * Date: Sun Jul 10 19:51:20 2016
# * Description: BRCA
# *              
# ***************************************

# read manifest file
d.manifest = read.csv(file = 'E:/DATA/TCGA-manifest-file-14/BRCA-gdc_manifest.2016-07-10T13-18-55.630543.tsv', sep = '\t',header = T)

# get a test files of XXX.FPKM.txt.gz
FPKM.gz.path = paste(datapath, 'TCGA-Data-14-20160704/',d.manifest$id[i.id],'/',d.manifest$filename[i.id], sep = '')
testdata = read.csv(file = FPKM.gz.path, sep = '\t', header = F)

# define a data.frame for all data 
BRCA.FPKM = data.frame(matrix(NA, nrow = nrow(testdata), ncol = nrow(d.manifest)), stringsAsFactors = F)
colnames(BRCA.FPKM) = apply(X = d.manifest, MARGIN = 1, FUN = function(aRow){substr(aRow[2], 1, 36)})
rownames(BRCA.FPKM) = testdata[,1]

# get all data form files
for(i.id in 1:nrow(d.manifest))
{
  
  FPKM.gz.path = paste(datapath, 'TCGA-Data-14-20160704/',d.manifest$id[i.id],'/',d.manifest$filename[i.id], sep = '')
  BRCA.FPKM[,i.id] = read.csv(file = FPKM.gz.path, sep = '\t', header = F)[,2]
  
  # remove file to another directions
  file.copy(from = FPKM.gz.path, to = 'E:/DATA/TCGA-Data/BRCA/',overwrite = F)
  print(i.id)
}

save(BRCA.FPKM, file = 'E:/DATA/TCGA-RDATA/BRCA.FPKM.RData')


########### 2016-07-15
load('E:/DATA/TCGA-RDATA/BRCA.FPKM.RData')

#save to mysql
con.brca = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.brca', user = 'root', password = 'root')

dbWriteTable(con.brca, name = "tcga_brca_fpkm", value = BRCA.FPKM, row.names = T, overwrite = F)

# get tcga gene list
#con.anno = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'annotation', user = 'root', password = 'root')
#tcga.gene.list = dbReadTable(conn = con.anno, name= "ensembl_id_tcga")[,1]

# get protein_coding list
#tcga.gene.map = dbReadTable(conn = con.anno, name= "tcga_gene_map")

# set name
rownames(BRCA.FPKM) = tcga.gene.list

# filter
BRCA.FPKM.GENE = BRCA.FPKM[tcga.gene.map$Ensembl.Gene.ID,] 
rownames(BRCA.FPKM.GENE) = tcga.gene.map$Gene.ID

# save to mysql and RDATA
#con.brca = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.brca', user = 'root', password = 'root')
dbWriteTable(con.brca, name = "tcga_brca_fpkm_gene", value = BRCA.FPKM.GENE, row.names = T, overwrite = T)
save(BRCA.FPKM.GENE, file = "E:/DATA/TCGA-RDATA/tcga_brca_fpkm_gene.RData")


##################### 2016-07-23: tidy clinical files

# read clinical manifest
d.manifest = read.csv(file = "E:/DATA/TCGA-Data/BRCA/clinical-manifest.txt", sep = '\t',header = T)

# copy file to clinical direction
apply(d.manifest, 1, copyfiles <- function(aRow) {
  filepath = paste0("E:/DATA/TCGA-20160723/",aRow[1], '/',aRow[2])
  file.copy(from = filepath, to = "E:/DATA/TCGA-Data/BRCA/Clinical/")
  file.remove(filepath)
} )


# ############ 2016-07-22: read and make clinical data.frame
# library(xml2)
# library(RMySQL)
# 
# # read clinical files
# v.file.names = dir(path = "E:/DATA/TCGA-Data/BRCA/Clinical", pattern = "*.xml", recursive = T, full.names = T)
# 
# # get a sample file to get clinical data.frame's colnames
# root = as_list(read_xml(v.file.names[1]))
# v.attr.name = names(root[[2]]) 
# v.attr.value = vector()
# 
# d.brca.clinical = data.frame(matrix(NA, nrow = length(v.file.names), ncol = length(v.attr.name)), stringsAsFactors = F)
# colnames(d.brca.clinical) = v.attr.name
# 
# # first loop: retrive clinical files
# for (i.file in 1:length(v.file.names)) {
#   #root = as_list(read_xml(v.file.names[i.file]))
#   tryCatch(expr = {root = as_list(read_xml(v.file.names[i.file]))},
#            warning = function(war){
#              print(war)
#             },error = function(err){
#               print(v.file.names[i.file])
#               next()
#            },finally = {})
#   
#   v.attr.name = names(root[[2]])
#   
#   # second loop: retrive every attribute and add to data.frame
#   for( name in v.attr.name) {
#     if(length(root[[2]][[name]]) == 1 && length(root[[2]][[name]][[1]]) == 1) {
#       d.brca.clinical[i.file,name] = root[[2]][[name]][[1]]
#     }else{
#       d.brca.clinical[i.file,name] = "NULL"
#     }
#   }
# }
# 
# #d.brca.clinical[i.file,]
# #v.file.names[i.file]

############################ 2016-07-29: read and make clinical data.frame
library(xml2)
library(RMySQL)

# read clinical files
v.file.names = dir(path = "E:/DATA/TCGA-Data/BRCA/Clinical", pattern = "*.xml", recursive = T, full.names = T)

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
d.brca.clinical = data.frame(matrix(NA, nrow = length(v.file.names), ncol = length(v.clinical.name)), stringsAsFactors = F)
colnames(d.brca.clinical) = v.clinical.name


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
  d.brca.clinical[i.file,] = v.xml.value[v.clinical.name]
}

# save file
save(d.brca.clinical, file = 'E:/DATA/TCGA-RDATA/BRCA/brca.clinical.RData')

con.brca = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.brca', user = 'root', password = 'root')
dbWriteTable(con.brca, name = "tcga_brca_clinical", value = d.brca.clinical, row.names = F, overwrite = T)

############ 2016-07-22: read metadata files, map patients id to mRNA expression file uuid
library(jsonlite)
library(dplyr)

metadata = fromJSON("E:/DATA/TCGA-Data/BRCA/mRNA-metadata.json")

d.brca.uuidmap = as.data.frame(t(data.frame(lapply(metadata$associated_entities, unlist))), stringsAsFactors = F)
d.brca.uuidmap = mutate(d.brca.uuidmap, mRNA_file_uuid = sapply(metadata$file_name,substr,1, 36,USE.NAMES = F))
rownames(d.brca.uuidmap) = NULL
## ***check****

load("E:/DATA/TCGA-RDATA/tcga_brca_fpkm_gene.RData")

## make the case same
d.brca.clinical$bcr_patient_uuid = tolower(d.brca.clinical$bcr_patient_uuid)
d.brca.uuidmap$case_id = tolower(d.brca.uuidmap$case_id)
colnames(BRCA.FPKM.GENE) = tolower(colnames(BRCA.FPKM.GENE))

## add barcode information to d.brca.uuidmap
d.brca.uuidmap = mutate(d.brca.uuidmap, 
                        patient_barcode = substr(d.brca.uuidmap$entity_submitter_id,1,12), 
                        tumor_sample = substr(d.brca.uuidmap$entity_submitter_id,14,15))

## split to 2 data.frame (Normal and Tumor)
d.brca.uuidmap.normal = filter(d.brca.uuidmap, tumor_sample == "11")
d.brca.uuidmap.tumor = filter(d.brca.uuidmap, tumor_sample == "01")

d.brca.fpkm.normal = BRCA.FPKM.GENE[,d.brca.uuidmap.normal$mRNA_file_uuid]
d.brca.fpkm.tumor = BRCA.FPKM.GENE[,d.brca.uuidmap.tumor$mRNA_file_uuid]

# change rowname (file uuid) by patient id
colnames(d.brca.fpkm.normal) = sapply(colnames(d.brca.fpkm.normal), function(a){filter(d.brca.uuidmap.normal, mRNA_file_uuid == a)$patient_barcode}, USE.NAMES = F)
colnames(d.brca.fpkm.tumor) = sapply(colnames(d.brca.fpkm.tumor), function(a){filter(d.brca.uuidmap.tumor, mRNA_file_uuid == a)$patient_barcode}, USE.NAMES = F)

## test if colnames have copies
nameids = colnames(d.brca.fpkm.tumor)

if(length(nameids) != length(unique(nameids)))
{
  print("have copies!!!!")
  
  # get copies id
  copyids = unique(nameids[duplicated(nameids)])
  
  id.col = which(nameids %in% copyids)
  d.id = d.brca.fpkm.tumor[,id.col]
  d.brca.fpkm.tumor = d.brca.fpkm.tumor[,-id.col]
  
  #orde by colnames
  d.id = d.id[,order(colnames(d.id))]
  copyids = sort(copyids)
  # deal it one by one 
  print(copyids)
  print(colnames(d.id))
  a1 = d.id[,c(1,2,3)] # "TCGA-A7-A0DB"
  a2 = d.id[,c(4,5)] # "TCGA-A7-A0DC"
  a3 = d.id[,c(6,7,8)] # "TCGA-A7-A13D"
  a4 = d.id[,c(9,10,11)] # "TCGA-A7-A13E"
  a5 = d.id[,c(12,13,14)] # "TCGA-A7-A26E"
  a6 = d.id[,c(15,16,17)] # "TCGA-A7-A26J"
  
  
  d.brca.fpkm.tumor = cbind.data.frame(d.brca.fpkm.tumor, "TCGA-A7-A0DB" = apply(a1, 1, mean))
  d.brca.fpkm.tumor = cbind.data.frame(d.brca.fpkm.tumor, "TCGA-A7-A0DC" = apply(a2, 1, mean))
  d.brca.fpkm.tumor = cbind.data.frame(d.brca.fpkm.tumor, "TCGA-A7-A13D" = apply(a3, 1, mean))
  d.brca.fpkm.tumor = cbind.data.frame(d.brca.fpkm.tumor, "TCGA-A7-A13E" = apply(a4, 1, mean))
  d.brca.fpkm.tumor = cbind.data.frame(d.brca.fpkm.tumor, "TCGA-A7-A26E" = apply(a5, 1, mean))
  d.brca.fpkm.tumor = cbind.data.frame(d.brca.fpkm.tumor, "TCGA-A7-A26J" = apply(a6, 1, mean))
  
}

nameids = colnames(d.brca.fpkm.normal)
if(length(nameids) != length(unique(nameids)))
{
  print("have copies!!!!")
}

# reduce clinical data
# if the column:  NULL data > 90% data, then remove the 

# save data 
save(d.brca.clinical, file = 'E:/DATA/TCGA-RDATA/BRCA/brca.clinical.RData')
save(d.brca.fpkm.normal, file = 'E:/DATA/TCGA-RDATA/BRCA/brca.fpkm.normal.RData')
save(d.brca.fpkm.tumor, file = 'E:/DATA/TCGA-RDATA/BRCA/brca.fpkm.tumor.RData')
save(d.brca.uuidmap, file = 'E:/DATA/TCGA-RDATA/BRCA/brca.uuidmap.RData')

con.brca = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.brca', user = 'root', password = 'root')
dbWriteTable(con.brca, name = "tcga_brca_clinical", value = d.brca.clinical, row.names = F, overwrite = T)
dbWriteTable(con.brca, name = "tcga_brca_fpkm_tumor", value = d.brca.fpkm.tumor, row.names = T, overwrite = T)


dbWriteTable(con.brca, name = "tcga_brca_fpkm_normal", value = d.brca.fpkm.normal, row.names = T, overwrite = T)
dbWriteTable(con.brca, name = "tcga_brca_uuidmap", value = d.brca.uuidmap, row.names = F, overwrite = T)


























