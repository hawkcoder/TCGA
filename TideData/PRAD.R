# ***************************************
# * Name: liangsen_2014@163.com
# * Date: Sun Jul 10 19:51:20 2016
# * Description: PRAD
# *              
# ***************************************

# define two file paths
datapath = 'E:/DATA/'

# read manifest file
d.manifest = read.csv(file = 'E:/DATA/TCGA-manifest-file-14/PRAD-gdc_manifest.2016-07-10T13-23-59.221838.tsv', sep = '\t',header = T)

# get a test files of XXX.FPKM.txt.gz
FPKM.gz.path = paste(datapath, 'TCGA-Data-14-20160704/',d.manifest$id[1],'/',d.manifest$filename[1], sep = '')
testdata = read.csv(file = FPKM.gz.path, sep = '\t', header = F)

# define a data.frame for all data 
PRAD.FPKM = data.frame(matrix(NA, nrow = nrow(testdata), ncol = nrow(d.manifest)), stringsAsFactors = F)
colnames(PRAD.FPKM) = apply(X = d.manifest, MARGIN = 1, FUN = function(aRow){substr(aRow[2], 1, 36)})
rownames(PRAD.FPKM) = testdata[,1]

# get all data form files
for(i.id in 1:nrow(d.manifest))
{
  
  FPKM.gz.path = paste(datapath, 'TCGA-Data-14-20160704/',d.manifest$id[i.id],'/',d.manifest$filename[i.id], sep = '')
  PRAD.FPKM[,i.id] = read.csv(file = FPKM.gz.path, sep = '\t', header = F)[,2]
  
  # remove file to another directions
  file.copy(from = FPKM.gz.path, to = 'E:/DATA/TCGA-Data/PRAD/',overwrite = F)
  print(i.id)
}

save(PRAD.FPKM, file = 'E:/DATA/TCGA-RDATA/PRAD.FPKM.RData')

########### 2016-07-15
load('E:/DATA/TCGA-RDATA/PRAD.FPKM.RData')

#save to mysql
con.prad = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.prad', user = 'root', password = 'root')

dbWriteTable(con.prad, name = "tcga_prad_fpkm", value = PRAD.FPKM, row.names = T, overwrite = F)

# get tcga gene list
#con.anno = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'annotation', user = 'root', password = 'root')
#tcga.gene.list = dbReadTable(conn = con.anno, name= "ensembl_id_tcga")[,1]

# get protein_coding list
#tcga.gene.map = dbReadTable(conn = con.anno, name= "tcga_gene_map")

# set name
rownames(PRAD.FPKM) = tcga.gene.list

# filter
PRAD.FPKM.GENE = PRAD.FPKM[tcga.gene.map$Ensembl.Gene.ID,] 
rownames(PRAD.FPKM.GENE) = tcga.gene.map$Gene.ID

# save to mysql and RDATA
#con.prad = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.prad', user = 'root', password = 'root')
dbWriteTable(con.prad, name = "tcga_prad_fpkm_gene", value = PRAD.FPKM.GENE, row.names = T, overwrite = T)
save(PRAD.FPKM.GENE, file = "E:/DATA/TCGA-RDATA/tcga_prad_fpkm_gene.RData")

##################### 2016-07-28: tidy clinical files

# read clinical manifest
d.manifest = read.csv(file = "E:/DATA/TCGA-Data/PRAD/clinical-manifest.tsv", sep = '\t',header = T)

# copy file to clinical direction
apply(d.manifest, 1, copyfiles <- function(aRow) {
  filepath = paste0("E:/DATA/TCGA-20160723/",aRow[1], '/',aRow[2])
  file.copy(from = filepath, to = "E:/DATA/TCGA-Data/PRAD/Clinical/")
  file.remove(filepath)
} )

############ 2016-07-28: read and make clinical data.frame
library(xml2)
library(RMySQL)

# read clinical files
v.file.names = dir(path = "E:/DATA/TCGA-Data/PRAD/Clinical", pattern = "*.xml", recursive = T, full.names = T)

# get a sample file to get clinical data.frame's colnames
root = as_list(read_xml(v.file.names[1]))
v.attr.name = names(root[[2]]) 
v.attr.value = vector()

d.prad.clinical = data.frame(matrix(NA, nrow = length(v.file.names), ncol = length(v.attr.name)), stringsAsFactors = F)
colnames(d.prad.clinical) = v.attr.name

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
      d.prad.clinical[i.file,name] = root[[2]][[name]][[1]]
    }else{
      d.prad.clinical[i.file,name] = "NULL"
    }
  }
}
# d.manifest[192,]
# v.file.names[192]

############ 2016-07-22: read metadata files, map patients id to mRNA expression file uuid
library(jsonlite)
library(dplyr)

metadata = fromJSON("E:/DATA/TCGA-Data/PRAD/mRNA-metadata.json")

d.prad.uuidmap = as.data.frame(t(data.frame(lapply(metadata$associated_entities, unlist))), stringsAsFactors = F)
d.prad.uuidmap = mutate(d.prad.uuidmap, mRNA_file_uuid = sapply(metadata$file_name,substr,1, 36,USE.NAMES = F))
rownames(d.prad.uuidmap) = NULL
## ***check****

load("E:/DATA/TCGA-RDATA/tcga_prad_fpkm_gene.RData")

## make the case same
d.prad.clinical$bcr_patient_uuid = tolower(d.prad.clinical$bcr_patient_uuid)
d.prad.uuidmap$case_id = tolower(d.prad.uuidmap$case_id)
colnames(PRAD.FPKM.GENE) = tolower(colnames(PRAD.FPKM.GENE))

## add barcode information to d.prad.uuidmap
d.prad.uuidmap = mutate(d.prad.uuidmap, 
                        patient_barcode = substr(d.prad.uuidmap$entity_submitter_id,1,12), 
                        tumor_sample = substr(d.prad.uuidmap$entity_submitter_id,14,15))

## split to 2 data.frame (Normal and Tumor)
d.prad.uuidmap.normal = filter(d.prad.uuidmap, tumor_sample == "11")
d.prad.uuidmap.tumor = filter(d.prad.uuidmap, tumor_sample == "01")

d.prad.fpkm.normal = PRAD.FPKM.GENE[,d.prad.uuidmap.normal$mRNA_file_uuid]
d.prad.fpkm.tumor = PRAD.FPKM.GENE[,d.prad.uuidmap.tumor$mRNA_file_uuid]

# change rowname (file uuid) by patient id
colnames(d.prad.fpkm.normal) = sapply(colnames(d.prad.fpkm.normal), function(a){filter(d.prad.uuidmap.normal, mRNA_file_uuid == a)$patient_barcode}, USE.NAMES = F)
colnames(d.prad.fpkm.tumor) = sapply(colnames(d.prad.fpkm.tumor), function(a){filter(d.prad.uuidmap.tumor, mRNA_file_uuid == a)$patient_barcode}, USE.NAMES = F)

## test if colnames have copies
nameids = colnames(d.prad.fpkm.tumor)
length(nameids) != length(unique(nameids))
if(length(nameids) != length(unique(nameids)))
{
  print("have copies!!!!")
  
  # get copies id
  copyids = unique(nameids[duplicated(nameids)])
  
  id.col = which(nameids %in% copyids)
  d.id = d.prad.fpkm.tumor[,id.col]
  d.prad.fpkm.tumor = d.prad.fpkm.tumor[,-id.col]
  
  #orde by colnames
  d.id = d.id[,order(colnames(d.id))]
  copyids = sort(copyids)
  # deal it one by one 
  print(copyids)
  print(colnames(d.id))
  
  # get v.id by hand
  v.id = list(c(1,2),c(3,4),c(6,5))
  
  # combind each copies to data.frame
  for(i in 1:3)
  {
    an = paste0("a", i)
    d.prad.fpkm.tumor = cbind.data.frame(d.prad.fpkm.tumor, "an" = apply(d.id[,v.id[[i]]], 1, mean))
    colnames(d.prad.fpkm.tumor)[ which(colnames(d.prad.fpkm.tumor) == "an")] = copyids[i]
  }
  
}

nameids = colnames(d.prad.fpkm.normal)
if(length(nameids) != length(unique(nameids)))
{
  print("have copies!!!!")
}

# reduce clinical data
# if the column:  NULL data > 90% data, then remove the 

# save data 
save(d.prad.clinical, file = 'E:/DATA/TCGA-RDATA/PRAD/prad.clinical.RData')
save(d.prad.fpkm.normal, file = 'E:/DATA/TCGA-RDATA/PRAD/prad.fpkm.normal.RData')
save(d.prad.fpkm.tumor, file = 'E:/DATA/TCGA-RDATA/PRAD/prad.fpkm.tumor.RData')
save(d.prad.uuidmap, file = 'E:/DATA/TCGA-RDATA/PRAD/prad.uuidmap.RData')

con.prad = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.prad', user = 'root', password = 'root')
dbWriteTable(con.prad, name = "tcga_prad_clinical", value = d.prad.clinical, row.names = F, overwrite = T)
dbWriteTable(con.prad, name = "tcga_prad_fpkm_tumor", value = d.prad.fpkm.tumor, row.names = T, overwrite = T)
dbWriteTable(con.prad, name = "tcga_prad_fpkm_normal", value = d.prad.fpkm.normal, row.names = T, overwrite = T)
dbWriteTable(con.prad, name = "tcga_prad_uuidmap", value = d.prad.uuidmap, row.names = F, overwrite = T)

