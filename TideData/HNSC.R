# ***************************************
# * Name: liangsen_2014@163.com
# * Date: Sun Jul 10 19:51:20 2016
# * Description: HNSC
# *              
# ***************************************

# define two file paths
datapath = 'E:/DATA/'

# read manifest file
d.manifest = read.csv(file = 'E:/DATA/TCGA-manifest-file-14/HNSC-gdc_manifest.2016-07-10T13-21-38.558570.tsv', sep = '\t',header = T)

# get a test files of XXX.FPKM.txt.gz
FPKM.gz.path = paste(datapath, 'TCGA-Data-14-20160704/',d.manifest$id[1],'/',d.manifest$filename[1], sep = '')
testdata = read.csv(file = FPKM.gz.path, sep = '\t', header = F)

# define a data.frame for all data 
HNSC.FPKM = data.frame(matrix(NA, nrow = nrow(testdata), ncol = nrow(d.manifest)), stringsAsFactors = F)
colnames(HNSC.FPKM) = apply(X = d.manifest, MARGIN = 1, FUN = function(aRow){substr(aRow[2], 1, 36)})
rownames(HNSC.FPKM) = testdata[,1]

# get all data form files
for(i.id in 1:nrow(d.manifest))
{
  
  FPKM.gz.path = paste(datapath, 'TCGA-Data-14-20160704/',d.manifest$id[i.id],'/',d.manifest$filename[i.id], sep = '')
  HNSC.FPKM[,i.id] = read.csv(file = FPKM.gz.path, sep = '\t', header = F)[,2]
  
  # remove file to another directions
  file.copy(from = FPKM.gz.path, to = 'E:/DATA/TCGA-Data/HNSC/',overwrite = F)
  print(i.id)
}

save(HNSC.FPKM, file = 'E:/DATA/TCGA-RDATA/HNSC.FPKM.RData')

########### 2016-07-15
load('E:/DATA/TCGA-RDATA/HNSC.FPKM.RData')

#save to mysql
con.hnsc = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.hnsc', user = 'root', password = 'root')

dbWriteTable(con.hnsc, name = "tcga_hnsc_fpkm", value = HNSC.FPKM, row.names = T, overwrite = F)

# get tcga gene list
#con.anno = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'annotation', user = 'root', password = 'root')
#tcga.gene.list = dbReadTable(conn = con.anno, name= "ensembl_id_tcga")[,1]

# get protein_coding list
#tcga.gene.map = dbReadTable(conn = con.anno, name= "tcga_gene_map")

# set name
rownames(HNSC.FPKM) = tcga.gene.list

# filter
HNSC.FPKM.GENE = HNSC.FPKM[tcga.gene.map$Ensembl.Gene.ID,] 
rownames(HNSC.FPKM.GENE) = tcga.gene.map$Gene.ID

# save to mysql and RDATA
#con.hnsc = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.hnsc', user = 'root', password = 'root')
dbWriteTable(con.hnsc, name = "tcga_hnsc_fpkm_gene", value = HNSC.FPKM.GENE, row.names = T, overwrite = T)
save(HNSC.FPKM.GENE, file = "E:/DATA/TCGA-RDATA/tcga_hnsc_fpkm_gene.RData")

##################### 2016-07-28: tidy clinical files

# read clinical manifest
d.manifest = read.csv(file = "E:/DATA/TCGA-Data/HNSC/clinical-manifest.txt", sep = '\t',header = T)

# copy file to clinical direction
apply(d.manifest, 1, copyfiles <- function(aRow) {
  filepath = paste0("E:/DATA/TCGA-20160723/",aRow[1], '/',aRow[2])
  file.copy(from = filepath, to = "E:/DATA/TCGA-Data/HNSC/Clinical/")
  file.remove(filepath)
} )

############ 2016-07-28: read and make clinical data.frame
library(xml2)
library(RMySQL)

# read clinical files
v.file.names = dir(path = "E:/DATA/TCGA-Data/HNSC/Clinical", pattern = "*.xml", recursive = T, full.names = T)

# get a sample file to get clinical data.frame's colnames
root = as_list(read_xml(v.file.names[1]))
v.attr.name = names(root[[2]]) 
v.attr.value = vector()

d.hnsc.clinical = data.frame(matrix(NA, nrow = length(v.file.names), ncol = length(v.attr.name)), stringsAsFactors = F)
colnames(d.hnsc.clinical) = v.attr.name

# first loop: retrive clinical files
for (i.file in 1:length(v.file.names)) {
  root = as_list(read_xml(v.file.names[i.file]))
  # tryCatch(expr = {root = as_list(read_xml(v.file.names[i.file]))},
  #          warning = function(war){
  #            print(war)
  #          },error = function(err){
  #            print(v.file.names[i.file])
  #            next()
  #          },finally = {next()})
  
  v.attr.name = names(root[[2]])
  
  # second loop: retrive every attribute and add to data.frame
  for( name in v.attr.name) {
    if(length(root[[2]][[name]]) == 1 && length(root[[2]][[name]][[1]]) == 1) {
      d.hnsc.clinical[i.file,name] = root[[2]][[name]][[1]]
    }else{
      d.hnsc.clinical[i.file,name] = "NULL"
    }
  }
}

#d.hnsc.clinical[i.file,]
#v.file.names[i.file]

############ 2016-07-22: read metadata files, map patients id to mRNA expression file uuid
library(jsonlite)
library(dplyr)

metadata = fromJSON("E:/DATA/TCGA-Data/HNSC/mRNA-metadata.json")

d.hnsc.uuidmap = as.data.frame(t(data.frame(lapply(metadata$associated_entities, unlist))), stringsAsFactors = F)
d.hnsc.uuidmap = mutate(d.hnsc.uuidmap, mRNA_file_uuid = sapply(metadata$file_name,substr,1, 36,USE.NAMES = F))
rownames(d.hnsc.uuidmap) = NULL
## ***check****

load("E:/DATA/TCGA-RDATA/tcga_hnsc_fpkm_gene.RData")

## make the case same
d.hnsc.clinical$bcr_patient_uuid = tolower(d.hnsc.clinical$bcr_patient_uuid)
d.hnsc.uuidmap$case_id = tolower(d.hnsc.uuidmap$case_id)
colnames(HNSC.FPKM.GENE) = tolower(colnames(HNSC.FPKM.GENE))

## add barcode information to d.hnsc.uuidmap
d.hnsc.uuidmap = mutate(d.hnsc.uuidmap, 
                        patient_barcode = substr(d.hnsc.uuidmap$entity_submitter_id,1,12), 
                        tumor_sample = substr(d.hnsc.uuidmap$entity_submitter_id,14,15))

## split to 2 data.frame (Normal and Tumor)
d.hnsc.uuidmap.normal = filter(d.hnsc.uuidmap, tumor_sample == "11")
d.hnsc.uuidmap.tumor = filter(d.hnsc.uuidmap, tumor_sample == "01")

d.hnsc.fpkm.normal = HNSC.FPKM.GENE[,d.hnsc.uuidmap.normal$mRNA_file_uuid]
d.hnsc.fpkm.tumor = HNSC.FPKM.GENE[,d.hnsc.uuidmap.tumor$mRNA_file_uuid]

# change rowname (file uuid) by patient id
colnames(d.hnsc.fpkm.normal) = sapply(colnames(d.hnsc.fpkm.normal), function(a){filter(d.hnsc.uuidmap.normal, mRNA_file_uuid == a)$patient_barcode}, USE.NAMES = F)
colnames(d.hnsc.fpkm.tumor) = sapply(colnames(d.hnsc.fpkm.tumor), function(a){filter(d.hnsc.uuidmap.tumor, mRNA_file_uuid == a)$patient_barcode}, USE.NAMES = F)

## test if colnames have copies
nameids = colnames(d.hnsc.fpkm.tumor)

if(length(nameids) != length(unique(nameids)))
{
  print("have copies!!!!")
  
  # get copies id
  copyids = unique(nameids[duplicated(nameids)])
  
  id.col = which(nameids %in% copyids)
  d.id = d.hnsc.fpkm.tumor[,id.col]
  d.hnsc.fpkm.tumor = d.hnsc.fpkm.tumor[,-id.col]
  
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
  #   d.hnsc.fpkm.tumor = cbind.data.frame(d.hnsc.fpkm.tumor, "an" = apply(d.id[,v.id[[i]]], 1, mean))
  #   colnames(d.hnsc.fpkm.tumor)[ which(colnames(d.hnsc.fpkm.tumor) == "an")] = copyids[i]
  # }
  
}

nameids = colnames(d.hnsc.fpkm.normal)
if(length(nameids) != length(unique(nameids)))
{
  print("have copies!!!!")
}

# reduce clinical data
# if the column:  NULL data > 90% data, then remove the 

# save data 
save(d.hnsc.clinical, file = 'E:/DATA/TCGA-RDATA/HNSC/hnsc.clinical.RData')
save(d.hnsc.fpkm.normal, file = 'E:/DATA/TCGA-RDATA/HNSC/hnsc.fpkm.normal.RData')
save(d.hnsc.fpkm.tumor, file = 'E:/DATA/TCGA-RDATA/HNSC/hnsc.fpkm.tumor.RData')
save(d.hnsc.uuidmap, file = 'E:/DATA/TCGA-RDATA/HNSC/hnsc.uuidmap.RData')

con.hnsc = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.hnsc', user = 'root', password = 'root')
dbWriteTable(con.hnsc, name = "tcga_hnsc_clinical", value = d.hnsc.clinical, row.names = F, overwrite = T)
dbWriteTable(con.hnsc, name = "tcga_hnsc_fpkm_tumor", value = d.hnsc.fpkm.tumor, row.names = T, overwrite = T)
dbWriteTable(con.hnsc, name = "tcga_hnsc_fpkm_normal", value = d.hnsc.fpkm.normal, row.names = T, overwrite = T)
dbWriteTable(con.hnsc, name = "tcga_hnsc_uuidmap", value = d.hnsc.uuidmap, row.names = F, overwrite = T)
