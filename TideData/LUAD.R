# ***************************************
# * Name: liangsen_2014@163.com
# * Date: Sun Jul 10 19:51:20 2016
# * Description: LUAD
# *              
# ***************************************

# define two file paths
datapath = 'E:/DATA/'

# read manifest file
d.manifest = read.csv(file = 'E:/DATA/TCGA-manifest-file-14/LUAD-gdc_manifest.2016-07-10T13-20-28.096456.tsv', sep = '\t',header = T)

# get a test files of XXX.FPKM.txt.gz
FPKM.gz.path = paste(datapath, 'TCGA-Data-14-20160704/',d.manifest$id[1],'/',d.manifest$filename[1], sep = '')
testdata = read.csv(file = FPKM.gz.path, sep = '\t', header = F)

# define a data.frame for all data 
LUAD.FPKM = data.frame(matrix(NA, nrow = nrow(testdata), ncol = nrow(d.manifest)), stringsAsFactors = F)
colnames(LUAD.FPKM) = apply(X = d.manifest, MARGIN = 1, FUN = function(aRow){substr(aRow[2], 1, 36)})
rownames(LUAD.FPKM) = testdata[,1]

# get all data form files
for(i.id in 1:nrow(d.manifest))
{
  
  FPKM.gz.path = paste(datapath, 'TCGA-Data-14-20160704/',d.manifest$id[i.id],'/',d.manifest$filename[i.id], sep = '')
  LUAD.FPKM[,i.id] = read.csv(file = FPKM.gz.path, sep = '\t', header = F)[,2]
  
  # remove file to another directions
  file.copy(from = FPKM.gz.path, to = 'E:/DATA/TCGA-Data/LUAD/',overwrite = F)
  print(i.id)
}

save(LUAD.FPKM, file = 'E:/DATA/TCGA-RDATA/LUAD.FPKM.RData')

########### 2016-07-15
load('E:/DATA/TCGA-RDATA/LUAD.FPKM.RData')

#save to mysql
con.luad = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.luad', user = 'root', password = 'root')

dbWriteTable(con.luad, name = "tcga_luad_fpkm", value = LUAD.FPKM, row.names = T, overwrite = F)

# get tcga gene list
#con.anno = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'annotation', user = 'root', password = 'root')
#tcga.gene.list = dbReadTable(conn = con.anno, name= "ensembl_id_tcga")[,1]

# get protein_coding list
#tcga.gene.map = dbReadTable(conn = con.anno, name= "tcga_gene_map")

# set name
rownames(LUAD.FPKM) = tcga.gene.list

# filter
LUAD.FPKM.GENE = LUAD.FPKM[tcga.gene.map$Ensembl.Gene.ID,] 
rownames(LUAD.FPKM.GENE) = tcga.gene.map$Gene.ID

# save to mysql and RDATA
#con.luad = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.luad', user = 'root', password = 'root')
dbWriteTable(con.luad, name = "tcga_luad_fpkm_gene", value = LUAD.FPKM.GENE, row.names = T, overwrite = T)
save(LUAD.FPKM.GENE, file = "E:/DATA/TCGA-RDATA/tcga_luad_fpkm_gene.RData")

##################### 2016-07-28: tidy clinical files

# read clinical manifest
d.manifest = read.csv(file = "E:/DATA/TCGA-Data/LUAD/clinical-manifest.tsv", sep = '\t',header = T)

# copy file to clinical direction
apply(d.manifest, 1, copyfiles <- function(aRow) {
  filepath = paste0("E:/DATA/TCGA-20160723/",aRow[1], '/',aRow[2])
  file.copy(from = filepath, to = "E:/DATA/TCGA-Data/LUAD/Clinical/")
  file.remove(filepath)
} )

############ 2016-07-28: read and make clinical data.frame
library(xml2)
library(RMySQL)

# read clinical files
v.file.names = dir(path = "E:/DATA/TCGA-Data/LUAD/Clinical", pattern = "*.xml", recursive = T, full.names = T)

# get a sample file to get clinical data.frame's colnames
root = as_list(read_xml(v.file.names[1]))
v.attr.name = names(root[[2]]) 
v.attr.value = vector()

d.luad.clinical = data.frame(matrix(NA, nrow = length(v.file.names), ncol = length(v.attr.name)), stringsAsFactors = F)
colnames(d.luad.clinical) = v.attr.name

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
      d.luad.clinical[i.file,name] = root[[2]][[name]][[1]]
    }else{
      d.luad.clinical[i.file,name] = "NULL"
    }
  }
}

#d.luad.clinical[i.file,]
#v.file.names[i.file]

############ 2016-07-22: read metadata files, map patients id to mRNA expression file uuid
library(jsonlite)
library(dplyr)

metadata = fromJSON("E:/DATA/TCGA-Data/LUAD/mRNA-metadata.json")

d.luad.uuidmap = as.data.frame(t(data.frame(lapply(metadata$associated_entities, unlist))), stringsAsFactors = F)
d.luad.uuidmap = mutate(d.luad.uuidmap, mRNA_file_uuid = sapply(metadata$file_name,substr,1, 36,USE.NAMES = F))
rownames(d.luad.uuidmap) = NULL
## ***check****

load("E:/DATA/TCGA-RDATA/tcga_luad_fpkm_gene.RData")

## make the case same
d.luad.clinical$bcr_patient_uuid = tolower(d.luad.clinical$bcr_patient_uuid)
d.luad.uuidmap$case_id = tolower(d.luad.uuidmap$case_id)
colnames(LUAD.FPKM.GENE) = tolower(colnames(LUAD.FPKM.GENE))

## add barcode information to d.luad.uuidmap
d.luad.uuidmap = mutate(d.luad.uuidmap, 
                        patient_barcode = substr(d.luad.uuidmap$entity_submitter_id,1,12), 
                        tumor_sample = substr(d.luad.uuidmap$entity_submitter_id,14,15))

## split to 2 data.frame (Normal and Tumor)
d.luad.uuidmap.normal = filter(d.luad.uuidmap, tumor_sample == "11")
d.luad.uuidmap.tumor = filter(d.luad.uuidmap, tumor_sample == "01")

d.luad.fpkm.normal = LUAD.FPKM.GENE[,d.luad.uuidmap.normal$mRNA_file_uuid]
d.luad.fpkm.tumor = LUAD.FPKM.GENE[,d.luad.uuidmap.tumor$mRNA_file_uuid]

# change rowname (file uuid) by patient id
colnames(d.luad.fpkm.normal) = sapply(colnames(d.luad.fpkm.normal), function(a){filter(d.luad.uuidmap.normal, mRNA_file_uuid == a)$patient_barcode}, USE.NAMES = F)
colnames(d.luad.fpkm.tumor) = sapply(colnames(d.luad.fpkm.tumor), function(a){filter(d.luad.uuidmap.tumor, mRNA_file_uuid == a)$patient_barcode}, USE.NAMES = F)

## test if colnames have copies
nameids = colnames(d.luad.fpkm.tumor)
length(nameids) != length(unique(nameids))
if(length(nameids) != length(unique(nameids)))
{
  print("have copies!!!!")
  
  # get copies id
  copyids = unique(nameids[duplicated(nameids)])
  
  id.col = which(nameids %in% copyids)
  d.id = d.luad.fpkm.tumor[,id.col]
  d.luad.fpkm.tumor = d.luad.fpkm.tumor[,-id.col]
  
  #orde by colnames
  d.id = d.id[,order(colnames(d.id))]
  copyids = sort(copyids)
  # deal it one by one 
  print(copyids)
  print(colnames(d.id))
  
  # get v.id by hand
  v.id = list(c(1,2,3),c(4,5,6),c(7,8,9),c(10,11,12),c(14,13),c(17,15,16),c(18,19),c(20,21,22),c(23,24,25),c(27,28,26),c(30,31,29))

  # combind each copies to data.frame
  for(i in 1:11)
  {
    an = paste0("a", i)
    d.luad.fpkm.tumor = cbind.data.frame(d.luad.fpkm.tumor, "an" = apply(d.id[,v.id[[i]]], 1, mean))
    colnames(d.luad.fpkm.tumor)[ which(colnames(d.luad.fpkm.tumor) == "an")] = copyids[i]
  }
  
}

nameids = colnames(d.luad.fpkm.normal)
if(length(nameids) != length(unique(nameids)))
{
  print("have copies!!!!")
}

# reduce clinical data
# if the column:  NULL data > 90% data, then remove the 

# save data 
save(d.luad.clinical, file = 'E:/DATA/TCGA-RDATA/LUAD/luad.clinical.RData')
save(d.luad.fpkm.normal, file = 'E:/DATA/TCGA-RDATA/LUAD/luad.fpkm.normal.RData')
save(d.luad.fpkm.tumor, file = 'E:/DATA/TCGA-RDATA/LUAD/luad.fpkm.tumor.RData')
save(d.luad.uuidmap, file = 'E:/DATA/TCGA-RDATA/LUAD/luad.uuidmap.RData')

con.luad = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.luad', user = 'root', password = 'root')
dbWriteTable(con.luad, name = "tcga_luad_clinical", value = d.luad.clinical, row.names = F, overwrite = T)
dbWriteTable(con.luad, name = "tcga_luad_fpkm_tumor", value = d.luad.fpkm.tumor, row.names = T, overwrite = T)
dbWriteTable(con.luad, name = "tcga_luad_fpkm_normal", value = d.luad.fpkm.normal, row.names = T, overwrite = T)
dbWriteTable(con.luad, name = "tcga_luad_uuidmap", value = d.luad.uuidmap, row.names = F, overwrite = T)

