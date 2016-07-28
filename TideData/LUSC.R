# ***************************************
# * Name: liangsen_2014@163.com
# * Date: Sun Jul 10 19:51:20 2016
# * Description: LUSC
# *              
# ***************************************

# define two file paths
datapath = 'E:/DATA/'

# read manifest file
d.manifest = read.csv(file = 'E:/DATA/TCGA-manifest-file-14/LUSC-gdc_manifest.2016-07-10T13-23-17.816264.tsv', sep = '\t',header = T)

# get a test files of XXX.FPKM.txt.gz
FPKM.gz.path = paste(datapath, 'TCGA-Data-14-20160704/',d.manifest$id[1],'/',d.manifest$filename[1], sep = '')
testdata = read.csv(file = FPKM.gz.path, sep = '\t', header = F)

# define a data.frame for all data 
LUSC.FPKM = data.frame(matrix(NA, nrow = nrow(testdata), ncol = nrow(d.manifest)), stringsAsFactors = F)
colnames(LUSC.FPKM) = apply(X = d.manifest, MARGIN = 1, FUN = function(aRow){substr(aRow[2], 1, 36)})
rownames(LUSC.FPKM) = testdata[,1]

# get all data form files
for(i.id in 1:nrow(d.manifest))
{
  
  FPKM.gz.path = paste(datapath, 'TCGA-Data-14-20160704/',d.manifest$id[i.id],'/',d.manifest$filename[i.id], sep = '')
  LUSC.FPKM[,i.id] = read.csv(file = FPKM.gz.path, sep = '\t', header = F)[,2]
  
  # remove file to another directions
  file.copy(from = FPKM.gz.path, to = 'E:/DATA/TCGA-Data/LUSC/',overwrite = F)
  #print(i.id)
}

save(LUSC.FPKM, file = 'E:/DATA/TCGA-RDATA/LUSC.FPKM.RData')

########### 2016-07-15
#load('E:/DATA/TCGA-RDATA/LUSC.FPKM.RData')

#save to mysql
con.lusc = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.lusc', user = 'root', password = 'root')

dbWriteTable(con.lusc, name = "tcga_lusc_fpkm", value = LUSC.FPKM, row.names = T, overwrite = F)

# get tcga gene list
#con.anno = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'annotation', user = 'root', password = 'root')
#tcga.gene.list = dbReadTable(conn = con.anno, name= "ensembl_id_tcga")[,1]

# get protein_coding list
#tcga.gene.map = dbReadTable(conn = con.anno, name= "tcga_gene_map")

# set name
rownames(LUSC.FPKM) = tcga.gene.list

# filter
LUSC.FPKM.GENE = LUSC.FPKM[tcga.gene.map$Ensembl.Gene.ID,] 
rownames(LUSC.FPKM.GENE) = tcga.gene.map$Gene.ID

# save to mysql and RDATA
#con.lusc = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.lusc', user = 'root', password = 'root')
dbWriteTable(con.lusc, name = "tcga_lusc_fpkm_gene", value = LUSC.FPKM.GENE, row.names = T, overwrite = T)
save(LUSC.FPKM.GENE, file = "E:/DATA/TCGA-RDATA/tcga_lusc_fpkm_gene.RData")

##################### 2016-07-28: tidy clinical files

# read clinical manifest
d.manifest = read.csv(file = "E:/DATA/TCGA-Data/LUSC/clinical-manifest.tsv", sep = '\t',header = T)

# copy file to clinical direction
apply(d.manifest, 1, copyfiles <- function(aRow) {
  filepath = paste0("E:/DATA/TCGA-20160723/",aRow[1], '/',aRow[2])
  file.copy(from = filepath, to = "E:/DATA/TCGA-Data/LUSC/Clinical/")
  file.remove(filepath)
} )

############ 2016-07-28: read and make clinical data.frame
library(xml2)
library(RMySQL)

# read clinical files
v.file.names = dir(path = "E:/DATA/TCGA-Data/LUSC/Clinical", pattern = "*.xml", recursive = T, full.names = T)

# get a sample file to get clinical data.frame's colnames
root = as_list(read_xml(v.file.names[1]))
v.attr.name = names(root[[2]]) 
v.attr.value = vector()

d.lusc.clinical = data.frame(matrix(NA, nrow = length(v.file.names), ncol = length(v.attr.name)), stringsAsFactors = F)
colnames(d.lusc.clinical) = v.attr.name

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
      d.lusc.clinical[i.file,name] = root[[2]][[name]][[1]]
    }else{
      d.lusc.clinical[i.file,name] = "NULL"
    }
  }
}

#d.lusc.clinical[i.file,]
#v.file.names[i.file]

############ 2016-07-22: read metadata files, map patients id to mRNA expression file uuid
library(jsonlite)
library(dplyr)

metadata = fromJSON("E:/DATA/TCGA-Data/LUSC/mRNA-metadata.json")

d.lusc.uuidmap = as.data.frame(t(data.frame(lapply(metadata$associated_entities, unlist))), stringsAsFactors = F)
d.lusc.uuidmap = mutate(d.lusc.uuidmap, mRNA_file_uuid = sapply(metadata$file_name,substr,1, 36,USE.NAMES = F))
rownames(d.lusc.uuidmap) = NULL
## ***check****

load("E:/DATA/TCGA-RDATA/tcga_lusc_fpkm_gene.RData")

## make the case same
d.lusc.clinical$bcr_patient_uuid = tolower(d.lusc.clinical$bcr_patient_uuid)
d.lusc.uuidmap$case_id = tolower(d.lusc.uuidmap$case_id)
colnames(LUSC.FPKM.GENE) = tolower(colnames(LUSC.FPKM.GENE))

## add barcode information to d.lusc.uuidmap
d.lusc.uuidmap = mutate(d.lusc.uuidmap, 
                        patient_barcode = substr(d.lusc.uuidmap$entity_submitter_id,1,12), 
                        tumor_sample = substr(d.lusc.uuidmap$entity_submitter_id,14,15))

## split to 2 data.frame (Normal and Tumor)
d.lusc.uuidmap.normal = filter(d.lusc.uuidmap, tumor_sample == "11")
d.lusc.uuidmap.tumor = filter(d.lusc.uuidmap, tumor_sample == "01")

d.lusc.fpkm.normal = LUSC.FPKM.GENE[,d.lusc.uuidmap.normal$mRNA_file_uuid]
d.lusc.fpkm.tumor = LUSC.FPKM.GENE[,d.lusc.uuidmap.tumor$mRNA_file_uuid]

# change rowname (file uuid) by patient id
colnames(d.lusc.fpkm.normal) = sapply(colnames(d.lusc.fpkm.normal), function(a){filter(d.lusc.uuidmap.normal, mRNA_file_uuid == a)$patient_barcode}, USE.NAMES = F)
colnames(d.lusc.fpkm.tumor) = sapply(colnames(d.lusc.fpkm.tumor), function(a){filter(d.lusc.uuidmap.tumor, mRNA_file_uuid == a)$patient_barcode}, USE.NAMES = F)

## test if colnames have copies
nameids = colnames(d.lusc.fpkm.tumor)
length(nameids) != length(unique(nameids))
if(length(nameids) != length(unique(nameids)))
{
  print("have copies!!!!")
  
  # get copies id
  copyids = unique(nameids[duplicated(nameids)])
  
  id.col = which(nameids %in% copyids)
  d.id = d.lusc.fpkm.tumor[,id.col]
  d.lusc.fpkm.tumor = d.lusc.fpkm.tumor[,-id.col]
  
  #orde by colnames
  d.id = d.id[,order(colnames(d.id))]
  copyids = sort(copyids)
  # deal it one by one 
  print(copyids)
  print(colnames(d.id))
  
  # get v.id by hand
  v.id = list(c(1,2))
  # combind each copies to data.frame
  for(i in 1:1)
  {
  #  i=1
    an = paste0("a", i)
    d.lusc.fpkm.tumor = cbind.data.frame(d.lusc.fpkm.tumor, "an" = apply(d.id[,v.id[[i]]], 1, mean))
    colnames(d.lusc.fpkm.tumor)[ which(colnames(d.lusc.fpkm.tumor) == "an")] = copyids[i]
  }
  
}

nameids = colnames(d.lusc.fpkm.normal)
if(length(nameids) != length(unique(nameids)))
{
  print("have copies!!!!")
}

# reduce clinical data
# if the column:  NULL data > 90% data, then remove the 

# save data 
save(d.lusc.clinical, file = 'E:/DATA/TCGA-RDATA/LUSC/lusc.clinical.RData')
save(d.lusc.fpkm.normal, file = 'E:/DATA/TCGA-RDATA/LUSC/lusc.fpkm.normal.RData')
save(d.lusc.fpkm.tumor, file = 'E:/DATA/TCGA-RDATA/LUSC/lusc.fpkm.tumor.RData')
save(d.lusc.uuidmap, file = 'E:/DATA/TCGA-RDATA/LUSC/lusc.uuidmap.RData')

con.lusc = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.lusc', user = 'root', password = 'root')
dbWriteTable(con.lusc, name = "tcga_lusc_clinical", value = d.lusc.clinical, row.names = F, overwrite = T)
dbWriteTable(con.lusc, name = "tcga_lusc_fpkm_tumor", value = d.lusc.fpkm.tumor, row.names = T, overwrite = T)
dbWriteTable(con.lusc, name = "tcga_lusc_fpkm_normal", value = d.lusc.fpkm.normal, row.names = T, overwrite = T)
dbWriteTable(con.lusc, name = "tcga_lusc_uuidmap", value = d.lusc.uuidmap, row.names = F, overwrite = T)


