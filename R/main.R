#' @author Yiming Zhao
#' @title main
#' @description main function in this package
#'
#' @param input_table_path input_table_path
#'
#' @param xmlurl xmlurl
#' 
#' @param min_intensity min_intensity=100
#'
#' @param cex cex
#'
#' @param srt srt
#'
#' @param ppm ppm
#'
#' @param PPM_denominator PPM_denominator
#'
#' @param pdf_width pdf_width
#'
#' @param pdf_height pdf_height
#' @param y_ion_col y ion's color
#' @param b_ion_col b ion's color
#' @return NULL
#' @import stringr
#' @import TeachingDemos
#' @import xml2
#' @importFrom MSnbase readMSData pickPeaks
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline axis lines par plot segments text
#' @importFrom data.table fread
#' @export

# Automatic drawing processing
main <-function(input_table_path,xmlurl,
                min_intensity=100,cex=1,srt=0,ppm=20,PPM_denominator=1E6,pdf_width=20,pdf_height=10,y_ion_col="red",b_ion_col="blue"){
  input_table <- data.table::fread(input_table_path, na.strings = "NA", sep = "\t",
                            check.names = FALSE, fill = TRUE, header = TRUE,
                            stringsAsFactors = FALSE)
  input_table <- as.data.frame(input_table)
  
  
  
  # check input_table
  if(nrow(input_table) == 0){stop("The variable input_table_path do NOT contain any MS2 information! Please see example files.")} # 
  input_table_split = split(input_table, input_table$label) # split table by label
  
  labelgroup=sapply(input_table_split, "[",1) # extract the first item in the split list 
  
  numbers_in_each_labelgroup = sapply(labelgroup, length) # MS2 Number for each label should be >=2 (mirror or groups)
  if(all(numbers_in_each_labelgroup >=2) == F){stop("Each group should contains at least 2 MS2 IDs!")}
  unique_number_in_each_labelgroup = length(sapply(labelgroup, unique)) #
  if(unique_number_in_each_labelgroup != length(input_table_split)){stop("MS2 IDs in each group should derive from the same raw MS file!")}
  
  # 
  
  mz2 = MSnbase::readMSData(input_table$`Raw file`[1],mode="onDisk") # open file
  
  

 
  
  # if(nrow(input_table)==0){
  #     stop("input_table is empty.")
  # }
  # # if(all(input_table$label == input_table$label[1])){
  # #     stop("Lable should not be unique.")
  # # }
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
          input_table$`Modified sequence`[i] = paste(
              substr(input_table$Sequence[i],1,mod_index),
              mod,substr(input_table$Sequence[i],
                         mod_index,nchar(input_table$Sequence[i])),sep = "")
      }
      if(input_table$mod[i] == ""){
          input_table$`Modifications`[i] = "Unmodified"
          input_table$`Modified sequence`[i] = NA
      }

      mod_or_unmod = msexp[[input_table$`Scan number`[i]]]

      input_table$Charge[i] = mod_or_unmod@precursorCharge
      input_table$`Retention time`[i] = mod_or_unmod@rt
      input_table$`m/z`[i] = mod_or_unmod@precursorMz


      mod_or_unmod = as.data.frame(mod_or_unmod)
      for(line in 1:nrow(mod_or_unmod)){
          input_table$mz[i] <- paste(input_table$mz[i],
                                     mod_or_unmod[line,1],sep=";")
          input_table$intensity[i] <- paste(input_table$intensity[i],
                                            mod_or_unmod[line,2],sep=";")
      }

  }
  # Remove the ; of the beginning
  input_table$mz = substr(input_table$mz, 2, nchar(input_table$mz))
  input_table$intensity = substr(input_table$intensity, 2,
                                 nchar(input_table$intensity))

  # After data processing, it begins to draw in units of label
  input_table_copy = input_table
  labels = unique(input_table$label)
  for (i in 1:length(labels)) {
    input_table_copy = input_table[c(input_table$label==labels[i]),]
    # If labels[i] is 2 in input_table, calling plot_mirror (), else call plot_parallel ()
    # rawfile=unlist(
    #   strsplit(input_table_copy[1,1], "\\\\|/|;|=", fixed = FALSE))[
    #     length(unlist(strsplit(input_table_copy[1,1], "\\\\|/|;|=",
    #                            fixed = FALSE)))]

    if(sum(unlist(input_table_copy$label)==labels[i])==2){
      plot_mirror(paste(input_table_copy$Sequence[1],"_mirror_",labels[i],sep = ""),
                  input_table_copy,xmlurl,min_intensity=min_intensity,cex=cex,srt=srt,ppm=ppm,PPM_denominator=PPM_denominator,pdf_width=pdf_width,pdf_height=pdf_height,y_ion_col=y_ion_col,b_ion_col=b_ion_col)
    }
    if(sum(unlist(input_table_copy$label)==labels[i])>2){
      plot_parallel(paste(input_table_copy$Sequence[1],"_parallel_",labels[i],sep = ""),
                    input_table_copy,xmlurl,min_intensity=min_intensity,cex=cex,srt=srt,ppm=ppm,PPM_denominator=PPM_denominator,pdf_width=pdf_width,pdf_height=pdf_height,y_ion_col=y_ion_col,b_ion_col=b_ion_col)
    }

  }
}
