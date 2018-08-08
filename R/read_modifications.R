#' @author Yiming Zhao
#'
#' @title read_modifications
#'
#' @description read modifications.xml file

#' @param xmlurl modifications.xml file path


#' @return modifications

#' @import stringr
#' @import xml2
#' @importFrom utils read.table
#' @export


read_modifications <-function(xmlurl){
  xmlread = read_xml(xmlurl)
  f1 = xml_find_all(xmlread,"//modifications/modification")
  f2 = xml_find_all(xmlread,"//modifications/modification/modification_site")
  modfication = xml_attr(f1,"composition")
  mod = xml_attr(f1,"title")
  mod = substr(mod,1,2)
  mod = tolower(mod)
  Modified_amino_acids = xml_attr(f2,"site")

  modifications = cbind(modfication,mod,Modified_amino_acids)
  modifications = data.frame(modifications)
  return(modifications)
}
