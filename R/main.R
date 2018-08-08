#' @author Yiming Zhao
#' @title main
#' @description main function in this package
#'
#' @param input_table_path input_table_path
#'
#' @param xmlurl xmlurl

#'
#' @return NULL
#' @import stringr
#' @import TeachingDemos
#' @import xml2
#' @importFrom MSnbase readMSData pickPeaks
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline axis lines par plot segments text
#' @importFrom utils read.table

#' @export

# 自动化绘图处理
main <-function(input_table_path,xmlurl){
  input_table <- read.table(input_table_path, na.strings = "NA", sep = "\t",
                            check.names = FALSE, fill = TRUE, header = TRUE,
                            stringsAsFactors = FALSE)
  mz2 = MSnbase::readMSData(input_table$`Raw file`[1],mode="onDisk") # open file
  for (i in 1:nrow(input_table)) {
    if(i==1){
      msexp <- pickPeaks(mz2)
    }else{
      if(input_table$`Raw file`[i] != input_table$`Raw file`[i-1]){
        mz2 = MSnbase::readMSData(input_table$`Raw file`[i],mode="onDisk")
        msexp <- pickPeaks(mz2)
      }
    }
    if(input_table$mod[i] != ""){

      mod = unlist(strsplit(input_table$mod[i],":"))[2]
      input_table$`Modifications`[i] = mod
      mod=paste("(",mod,")", sep="")
      mod_index = unlist(strsplit(input_table$mod[i],":"))[1]
      mod_index = as.numeric(substr(mod_index,2,nchar(mod_index)))
      input_table$`Modified sequence`[i] =
        paste(substr(input_table$Sequence[i],1,mod_index),mod,
              substr(input_table$Sequence[i],mod_index,
                     nchar(input_table$Sequence[i])),sep = "")
    }
    if(input_table$mod[i] == ""){
      input_table$`Modifications`[i] = "Unmodified"
    }

    mod_or_unmod = msexp[[input_table$`Scan number`[i]]]

    input_table$Charge[i] = mod_or_unmod@precursorCharge
    input_table$`Retention time`[i] = mod_or_unmod@rt
    input_table$`m/z`[i] = mod_or_unmod@precursorMz


    mod_or_unmod = as.data.frame(mod_or_unmod)
    for(line in 1:nrow(mod_or_unmod)){
      input_table$mz[i] <- paste(input_table$mz[i],mod_or_unmod[line,1],sep=";")
      input_table$intensity[i] <- paste(input_table$intensity[i],
                                        mod_or_unmod[line,2],sep=";")
    }

  }
  # 去掉开始的;号
  input_table$mz = substr(input_table$mz, 2, nchar(input_table$mz))
  input_table$intensity = substr(input_table$intensity, 2,
                                 nchar(input_table$intensity))

  # 数据处理完后开始以label为单位绘图
  input_table_copy = input_table
  labels = unique(input_table$label)
  for (i in 1:length(labels)) {
    input_table_copy = input_table[c(input_table$label==labels[i]),]
    # 若labels[i]在input_table数量是2，调用plot_mirror(),否调用plot_parallel()
    rawfile=unlist(
      strsplit(input_table_copy[1,1], "\\\\|/|;|=", fixed = FALSE))[
        length(unlist(strsplit(input_table_copy[1,1], "\\\\|/|;|=",
                               fixed = FALSE)))]

    if(sum(unlist(input_table_copy$label)==labels[i])==2){
      plot_mirror(paste(rawfile,"_mirror",labels[i],Sys.time()),
                  input_table_copy,xmlurl)
    }
    if(sum(unlist(input_table_copy$label)==labels[i])>2){
      plot_parallel(paste(rawfile,"_parallel",labels[i],Sys.time()),
                    input_table_copy,xmlurl)
    }

  }
}
