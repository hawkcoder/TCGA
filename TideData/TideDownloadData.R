# ***************************************
# * Name: liangsen.2014@163.com
# * Date: 2016-07-10 14:24:54 
# * Description: null 
# *              
# ***************************************

# define two file paths
datapath = 'E:/DATA/'
codepath = 'E:/CODE/TCGA/TideData/'

# 14 tumors' name
v.tumor.name = c("BRCA","LUAD","KIRC","HNSC","THCA","LUSC","PRAD","STAD","COAD","BLCA","LIHC","KIRP","KICH","ESCA") 

# new directions file
for (name in v.tumor.name) {
  dir.create(paste0('E:/DATA/TCGA-RData/',name), showWarnings = T)
}



# read manifest file
v.manifest.file.names = dir(path = paste(datapath,'TCGA-manifest-file-14/', sep = ''))

# First Loop: each manifest file
for(i.maf in 1:nrow(v.manifest.file.names))
{
  
  tumor.abbrev = unlist(strsplit(v.manifest.file.names[i.maf], '-'))[2]    # tumor name abbreviate like "KIRC","LIHC"
  manifest.file.path = paste(datapath,'TCGA-manifest-file-14/',v.manifest.file.names[i.maf], sep = "")
  d.manifest = read.csv(file = manifest.file.path, sep = '\t',header = T)
  
}

# get a test files of XXX.FPKM.txt.gz
FPKM.gz.path = paste(datapath, 'TCGA-Data-14-20160704/',d.manifest$id[i.id],'/',d.manifest$filename[i.id], sep = '')
testdata = read.csv(file = FPKM.gz.path, sep = '\t', header = F)

# define a data.frame for all data 
KIRC.FPKM = data.frame(matrix(NA, nrow = nrow(testdata), ncol = nrow(d.manifest)), stringsAsFactors = F)
colnames(KIRC.FPKM) = apply(X = d.manifest, MARGIN = 1, FUN = function(aRow){substr(aRow[2], 1, 36)})
rownames(KIRC.FPKM) = testdata[,1]

# get all data form files
for(i.id in 1:nrow(d.manifest))
{
  FPKM.gz.path = paste(datapath, 'TCGA-Data-14-20160704/',d.manifest$id[i.id],'/',d.manifest$filename[i.id], sep = '')
  KIRC.FPKM[,i.id] = read.csv(file = FPKM.gz.path, sep = '\t', header = F)[,2]
}


##################### 2016-07-23: tidy clinical files

# get all clinical files 
v.clinical.files = dir(path = "E:/DATA/TCGA-Data/READ/Clinical/")
d.manifest = read.csv(file = "E:/DATA/TCGA-Data/READ/Clinical/clinical-manifest.tsv", header = T, sep = "\t", stringsAsFactors = F)

d.manifest = d.manifest[which( !(d.manifest$id %in% v.clinical.files)), ]
write.table(d.manifest, file = "E:/DATA/TCGA-Data/READ/Clinical/clinical-manifest-2.tsv", quote = F, row.names = F, sep = "\t")

# get all mRNA files 
v.clinical.files = dir(path = "E:/DATA/TCGA-Data/STAD/Clinical/")
d.manifest = read.csv(file = "E:/DATA/TCGA-Data/STAD/clinical-manifest.tsv", header = T, sep = "\t", stringsAsFactors = F)

d.manifest = d.manifest[which( !(d.manifest$id %in% v.clinical.files)), ]
write.table(d.manifest, file = "E:/DATA/TCGA-Data/STAD/clinical-manifest-1.tsv", quote = F, row.names = F, sep = "\t")



#################### 2016-07-28
v.tumor.name = c("BRCA","LUAD","KIRC","HNSC","THCA","LUSC","PRAD","STAD","COAD","BLCA","LIHC","KIRP","KICH","ESCA") 

v.tumor.name = c("BRCA","COAD","ESCA","HNSC","KICH","KIRP","LIHC","LUAD")
for(name in v.tumor.name)
{
  xxxx.fpkm.all.gene = paste0(tolower(name), ".fpkm.all.gene")
  xxxx.fpkm = paste0(tolower(name), ".fpkm")
  
  load( paste0("E:/DATA/TCGA-RData/tcga_",tolower(name), "_fpkm_gene.RData"))
  load( paste0("E:/DATA/TCGA-RData/",name,".FPKM.RData"))
  
  assign(xxxx.fpkm.all.gene, get(paste0(name, ".FPKM")))
  assign(xxxx.fpkm, get(paste0(name,".FPKM.GENE")))
  
  get(xxxx.fpkm.all.gene) %>% save(file = paste0("E:/DATA/TCGA-RData/",name,"/", xxxx.fpkm.all.gene, ".RData"))
  get(xxxx.fpkm) %>% save(file = paste0("E:/DATA/TCGA-RData/",name,"/", xxxx.fpkm, ".RData"))
  
  file.remove(paste0("E:/DATA/TCGA-RData/tcga_",tolower(name), "_fpkm_gene.RData"))
  file.remove(paste0("E:/DATA/TCGA-RData/",name,".FPKM.RData"))
}

