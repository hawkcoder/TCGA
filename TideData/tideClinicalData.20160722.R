library(xml2)

data <- read_xml("E:/DATA/TCGA-Data-20160722/gdc_download_20160722_051255/00b31387-a3b4-4181-9c23-f85438c580f0/nationwidechildrens.org_clinical.TCGA-B0-4843.xml")

data2 = xmlParse("E:/DATA/TCGA-Data-20160722/gdc_download_20160722_051255/8e14979d-4f1b-4694-8d87-a475bf993e9a/nationwidechildrens.org_clinical.TCGA-CJ-4638.xml")


root = as_list(data)

names(root[[1]])
names(root[[2]][["days_to_last_known_alive"]])
root[[2]][["days_to_last_known_alive"]][[1]]
root[[2]][["patient_id"]][[1]]

mode(root[[2]][["patient_id"]])
mode(root[[2]][["days_to_last_known_alive"]])

v.attr.name = names(root[[2]])
v.attr.value = character(length = length(v.attr.name))
names(v.attr.value) = v.attr.name

for( name in v.attr.name) {
  if(length(root[[2]][[name]]) == 1 && length(root[[2]][[name]][[1]]) == 1) {
    v.attr.value[name] = root[[2]][[name]][[1]]
  }else{
    v.attr.value[name] = "NULL"
  }
}



v.attr.value

length(root[[2]][["new_tumor_events"]])

root[[2]][["new_tumor_events"]][[1]]


xml_data$patient$patient_id
a = "patient_id"
xml_data$patient[[a]][["text"]]
xml_data$patient[["additional_studies"]]

v.attr.value[name] = "dd"
xml_data[[name]][["text"]]
name

############################### 2016-07-29 read a clinical file with all its elements
library(xml2)
data.xml = read_xml("E:/DATA/TCGA-Data/BRCA/Clinical/nationwidechildrens.org_clinical.TCGA-3C-AAAU.xml")

all.node = xml_find_all(data.xml, "//*")

v.name = vector()
v.value = vector()
for (lst in all.node) {
  # if it is a leaf element, then leaves
  if(length(xml_find_all(lst, ".//*")) == 0)
  {
    v.name = c(v.name, xml_name(lst))
    v.value = c(v.value, xml_text(lst))
  }
}
d.clinical = data.frame(name = v.name, value = v.value, stringsAsFactors = F)

# delete name which have copies
copies = unique(v.name[duplicated(v.name)])

v.name.uq = v.name[which(!v.name %in% copies)]
v.value.uq = v.value[which(!v.name %in% copies)]

a1 = xml_path(all.node)
cs = unique(a1[duplicated(a1)])

a1.uq = a1[which(!v.name %in% copies)]

a3 = xml_name(all.node)
a3%in% v.name.uq

d.clinical = d.clinical[which(!d.clinical$name %in% copies),]

a1.da = xml_find_all(data.xml, a1.uq[1])

############################ 2016-07-29: read and make clinical data.frame
library(xml2)
library(RMySQL)

# read clinical files
v.file.names = dir(path = "E:/DATA/TCGA-Data/COAD/Clinical", pattern = "*.xml", recursive = T, full.names = T)

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
d.coad.clinical = data.frame(matrix(NA, nrow = length(v.file.names), ncol = length(v.clinical.name)), stringsAsFactors = F)
colnames(d.coad.clinical) = v.clinical.name


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
  d.coad.clinical[i.file,] = v.xml.value[v.clinical.name]
}

# save file
library(RMySQL)
save(d.coad.clinical, file = 'E:/DATA/TCGA-RDATA/COAD/coad.clinical.RData')

con.coad = dbConnect(MySQL(), host = '127.0.0.1', dbname = 'tcga.coad', user = 'root', password = 'root')
dbWriteTable(con.coad, name = "tcga_coad_clinical", value = d.coad.clinical, row.names = F, overwrite = T)
