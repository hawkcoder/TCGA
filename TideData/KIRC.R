# ***************************************
# * Name: liangsen_2014@163.com
# * Date: Sun Jul 10 19:51:20 2016
# * Description: KIRC
# *              
# ***************************************

# define two file paths
datapath = 'E:/DATA/'
codepath = 'E:/CODE/TCGA/TideData/'

# read manifest file
d.manifest = read.csv(file = 'E:/DATA/TCGA-manifest-file-14/TCGA-KIRC-gdc_manifest.2016-07-10T07-46-16.482076.tsv', sep = '\t',header = T)

# get a test files of XXX.FPKM.txt.gz
FPKM.gz.path = paste(datapath, 'TCGA-Data-14-20160704/',d.manifest$id[i.id],'/',d.manifest$filename[i.id], sep = '')
testdata = read.csv(file = FPKM.gz.path, sep = '\t', header = F)

# define a data.frame for all data 
KIRC.FPKM = data.frame(matrix(NA, nrow = nrow(testdata), ncol = nrow(d.manifest)), stringsAsFactors = F)
colnames(KIRC.FPKM) = apply(X = d.manifest, MARGIN = 1, FUN = function(aRow){substr(aRow[2], 1, 36)})
rownames(KIRC.FPKM) = testdata[,1]

# get all data form files
for(i.id in 2:nrow(d.manifest))
{
  
  FPKM.gz.path = paste(datapath, 'TCGA-Data-14-20160704/',d.manifest$id[i.id],'/',d.manifest$filename[i.id], sep = '')
  KIRC.FPKM[,i.id] = read.csv(file = FPKM.gz.path, sep = '\t', header = F)[,2]
  
  # remove file to another directions
  file.copy(from = FPKM.gz.path, to = 'E:/DATA/TCGA-Data/KIRC/',overwrite = F)
  print(i.id)
}

save(KIRC.FPKM, file = 'E:/DATA/TCGA-RDATA/KIRC.FPKM.RData')

########### 2016-07-15
load('E:/DATA/TCGA-RDATA/KIRC.FPKM.RData')

#save to mysql
con.kirc = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.kirc', user = 'root', password = 'root')

dbWriteTable(con.kirc, name = "tcga_kirc_fpkm", value = KIRC.FPKM, row.names = T, overwrite = F)

# get tcga gene list
#con.anno = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'annotation', user = 'root', password = 'root')
#tcga.gene.list = dbReadTable(conn = con.anno, name= "ensembl_id_tcga")[,1]

# get protein_coding list
#tcga.gene.map = dbReadTable(conn = con.anno, name= "tcga_gene_map")

# set name
rownames(KIRC.FPKM) = tcga.gene.list

# filter
KIRC.FPKM.GENE = KIRC.FPKM[tcga.gene.map$Ensembl.Gene.ID,] 
rownames(KIRC.FPKM.GENE) = tcga.gene.map$Gene.ID

# save to mysql and RDATA
#con.kirc = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.kirc', user = 'root', password = 'root')
dbWriteTable(con.kirc, name = "tcga_kirc_fpkm_gene", value = KIRC.FPKM.GENE, row.names = T, overwrite = T)
save(KIRC.FPKM.GENE, file = "E:/DATA/TCGA-RDATA/tcga_kirc_fpkm_gene.RData")

# ############ 2016-07-22: read and make clinical data.frame
# library(xml2)
# library(RMySQL)
# 
# datapath = "E:/DATA/"
# 
# # read clinical files
# v.file.names = dir(path = paste0(datapath, "TCGA-Data/KIRC/Clinical"), pattern = "*.xml", recursive = T, full.names = T)
# 
# # get a sample file to get clinical data.frame's colnames
# root = as_list(read_xml(v.file.names[1]))
# v.attr.name = names(root[[2]]) 
# v.attr.value = vector()
# 
# d.kirc.clinical = data.frame(matrix(NA, nrow = length(v.file.names), ncol = length(v.attr.name)), stringsAsFactors = F)
# colnames(d.kirc.clinical) = v.attr.name
# 
# # first loop: retrive clinical files
# for (i.file in 1:length(v.file.names)) {
#   root = as_list(read_xml(v.file.names[i.file]))
#   v.attr.name = names(root[[2]])
#   
#   # second loop: retrive every attribute and add to data.frame
#   for( name in v.attr.name) {
#     if(length(root[[2]][[name]]) == 1 && length(root[[2]][[name]][[1]]) == 1) {
#       d.kirc.clinical[i.file,name] = root[[2]][[name]][[1]]
#     }else{
#       d.kirc.clinical[i.file,name] = "NULL"
#     }
#   }
# }
# 
# # save data 
# save(d.kirc.clinical, file = 'E:/DATA/TCGA-RDATA/KIRC.CLINICAL.RData')
# 
# con.kirc = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.kirc', user = 'root', password = 'root')
# dbWriteTable(con.kirc, name = "tcga_kirc_clinical", value = d.kirc.clinical, row.names = F, overwrite = T)

############################ 2016-07-29: read and make clinical data.frame
library(xml2)
library(RMySQL)

# read clinical files
v.file.names = dir(path = "E:/DATA/TCGA-Data/KIRC/Clinical", pattern = "*.xml", recursive = T, full.names = T)

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
d.kirc.clinical = data.frame(matrix(NA, nrow = length(v.file.names), ncol = length(v.clinical.name)), stringsAsFactors = F)
colnames(d.kirc.clinical) = v.clinical.name


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
  d.kirc.clinical[i.file,] = v.xml.value[v.clinical.name]
}

# save file
save(d.kirc.clinical, file = 'E:/DATA/TCGA-RDATA/KIRC/kirc.clinical.RData')

con.kirc = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.kirc', user = 'root', password = 'root')
dbWriteTable(con.kirc, name = "tcga_kirc_clinical", value = d.kirc.clinical, row.names = F, overwrite = T)

############ 2016-07-22: read metadata files, map patients id to mRNA expression file uuid
install.packages("jsonlite")
library(jsonlite)

metadata = fromJSON("E:/DATA/TCGA-Data-20160722/metadata.cart.2016-07-22T08-57-29.078449.json")

d.uuidmap = as.data.frame(t(data.frame(lapply(metadata$associated_entities, unlist))), stringsAsFactors = F)
d.uuidmap = mutate(d.uuidmap, mRNA_file_uuid = sapply(metadata$file_name,substr,1, 36,USE.NAMES = F))
rownames(d.uuidmap) = NULL


load("E:/DATA/TCGA-RDATA/tcga_kirc_fpkm_gene.RData")

## make the case same
d.kirc.clinical$bcr_patient_uuid = tolower(d.kirc.clinical$bcr_patient_uuid)
d.uuidmap$case_id = tolower(d.uuidmap$case_id)
colnames(KIRC.FPKM.GENE) = tolower(colnames(KIRC.FPKM.GENE))

## add barcode information to d.uuidmap
library(dplyr)
d.uuidmap = tbl_df(d.uuidmap)

d.uuidmap = mutate(d.uuidmap, 
                   patient_barcode = substr(d.uuidmap$entity_submitter_id,1,12), 
                   tumor_sample = substr(d.uuidmap$entity_submitter_id,14,15))

## split to 2 data.frame (Normal and Tumor)
d.uuidmap.normal = filter(d.uuidmap, tumor_sample >= 10)
d.uuidmap.tumor = filter(d.uuidmap, tumor_sample < 10)

kirc.fpkm.normal = KIRC.FPKM.GENE[,d.uuidmap.normal$mRNA_file_uuid]
kirc.fpkm.tumor = KIRC.FPKM.GENE[,d.uuidmap.tumor$mRNA_file_uuid]

# change rowname (file uuid) by patient id
colnames(kirc.fpkm.normal) = sapply(colnames(kirc.fpkm.normal), function(a){filter(d.uuidmap.normal, mRNA_file_uuid == a)$patient_barcode}, USE.NAMES = F)
colnames(kirc.fpkm.tumor) = sapply(colnames(kirc.fpkm.tumor), function(a){filter(d.uuidmap.tumor, mRNA_file_uuid == a)$patient_barcode}, USE.NAMES = F)

# "TCGA-B2-5635" have 3 FPKM files, Now combine it
TCGA.B2.5635 = d.kirc.fpkm.tumor[,c(45,231,352)]
d.kirc.fpkm.tumor = d.kirc.fpkm.tumor[,-c(231,352)]
d.kirc.fpkm.tumor[,45] = apply(TCGA.B2.5635, MARGIN = 1,FUN = mean)

# save data 
d.kirc.uuidmap = d.uuidmap
d.kirc.fpkm.normal = kirc.fpkm.normal
d.kirc.fpkm.tumor = kirc.fpkm.tumor

save(d.kirc.clinical, file = 'E:/DATA/TCGA-RDATA/KIRC/kirc.clinical.RData')
save(d.kirc.fpkm.normal, file = 'E:/DATA/TCGA-RDATA/KIRC/kirc.fpkm.normal.RData')
save(d.kirc.fpkm.tumor, file = 'E:/DATA/TCGA-RDATA/KIRC/kirc.fpkm.tumor.RData')
save(d.kirc.uuidmap, file = 'E:/DATA/TCGA-RDATA/KIRC/kirc.uuidmap.RData')

con.kirc = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.kirc', user = 'root', password = 'root')
dbWriteTable(con.kirc, name = "tcga_kirc_clinical", value = d.kirc.clinical, row.names = F, overwrite = T)
dbWriteTable(con.kirc, name = "tcga_kirc_fpkm_normal", value = d.kirc.fpkm.normal, row.names = T, overwrite = T)
dbWriteTable(con.kirc, name = "tcga_kirc_fpkm_tumor", value = d.kirc.fpkm.tumor, row.names = T, overwrite = T)
dbWriteTable(con.kirc, name = "tcga_kirc_uuidmap", value = d.kirc.uuidmap, row.names = F, overwrite = T)
