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
