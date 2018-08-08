#' @author Yiming Zhao
#' @title plot_mirror
#' @description plot_mirror function in this package
#'
#' @param fileName filename
#'
#' @param f.msms f.msms
#'
#' @param xmlurl modifications.xml path
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



# 镜像绘图方法
plot_mirror <- function(fileName,f.msms,xmlurl,min_intensity=100,cex=1,srt=0,
                        ppm=20,PPM_denominator=1E6,pdf_width=20,pdf_height=10){
  pdf(file=paste(fileName,sep="",".pdf"), width=pdf_width, height=pdf_height)

  # 处理 modfication.xml 文件，
  # 并读取 site ，title ，composition 处理后并入f.unimod
  modifications<-read_modifications(xmlurl)

  f.atom_MW = structure(list(Element = structure(c(20L, 4L, 24L, 9L, 11L, 1L,
  28L, 2L, 31L, 3L, 18L, 29L, 25L, 6L, 32L, 36L, 14L, 23L, 12L,
  16L, 26L, 19L, 30L, 15L, 17L, 38L, 7L, 10L, 37L, 27L, 35L, 33L,
  5L, 13L, 22L, 34L, 8L, 21L), .Label = c("13C", "15N", "18O",
  "2H", "Ag", "Al", "As", "Au", "B", "Br", "C", "Ca", "Cd", "Cl",
  "Co", "Cr", "Cu", "F", "Fe", "H", "Hg", "I", "K", "Li", "Mg",
  "Mn", "Mo", "N", "Na", "Ni", "O", "P", "Pd", "Pt", "Ru", "S",
  "Se", "Zn"), class = "factor"), Monoisotopic = c(1.007825035,
  2.014101779, 7.016003, 11.0093055, 12, 13.00335483, 14.003074,
  15.00010897, 15.99491463, 17.9991603, 18.99840322, 22.9897677,
  23.9850423, 26.9815386, 30.973762, 31.9720707, 34.96885272, 38.9637074,
  39.9625906, 51.9405098, 54.9380471, 55.9349393, 57.9353462, 58.9331976,
  62.9295989, 63.9291448, 74.9215942, 78.9183361, 79.9165196, 97.9054073,
  101.9043485, 105.903478, 106.905092, 113.903357, 126.904473,
  194.964766, 196.966543, 201.970617)), class = "data.frame", row.names = c(NA,
   -38L))

  f.unimod = structure(list(AA = structure(1:25,
  .Label = c("A", "B", "C",
  "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "P", "Q",
  "R", "S", "T", "U", "V", "W", "X", "Y", "Z"), class = "factor"),
  weight = c(71.03711, 111.03203, 103.00919, 115.02694, 129.04259,
 147.06841, 57.02146, 137.05891, 113.08406, 173.03242, 128.09496,
 113.08406, 131.04049, 114.04293, 97.05276, 128.05858, 156.10111,
 87.03203, 101.04768, 113.08406, 99.06841, 186.07931, 113.04768,
 163.06333, 103.00919), location = structure(c(1L, 1L, 1L,
 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = "M", class = "factor")),
 class = "data.frame", row.names = c(NA,-25L))

  modfication = as.character(modifications$modfication[1])
  monoisotopic = f.atom_MW$Monoisotopic[which(f.atom_MW$Element==modfication)]
  mod = modifications$mod[1]
  mod=paste("(",mod,")", sep="")
  f.unimod$weight = as.numeric(f.unimod$weight)

  l = vector()
  for (i in 1:nrow(f.unimod)){
    for (j in 1:nrow(modifications)){
      if(f.unimod$AA[i] == modifications$Modified_amino_acids[j]){
        l = c(l,i)
      }
    }
  }
  f = f.unimod[l,]
  f$AA =  as.factor(paste(f$AA,mod,sep = ''))
  f$weight = f$weight + monoisotopic
  f$location = f$location
  f.unimod = rbind(f.unimod,f)

  H=f.atom_MW$Monoisotopic[which(f.atom_MW$Element=="H")]
  O=f.atom_MW$Monoisotopic[which(f.atom_MW$Element=="O")]
  H2O=H*2+O


  f.msms_u = subset(f.msms, Modifications %in% "Unmodified")
  f.msms_m = subset(f.msms, ! Modifications %in% "Unmodified")

  l1.masses = strsplit(f.msms_m$mz,split=";")
  l1.intensity = strsplit(f.msms_m$intensity,split=";")
  mod.ms2 = data.frame(l1.masses, l1.intensity, stringsAsFactors = F)
  colnames(mod.ms2) = c('mass','intensity')

  l1.masses = strsplit(f.msms_u$mz,split=";")
  l1.intensity = strsplit(f.msms_u$intensity,split=";")
  unmod.ms2 = data.frame(l1.masses, l1.intensity, stringsAsFactors = F)
  colnames(unmod.ms2) = c('mass','intensity')
  mod.ms2$mass<-as.numeric(mod.ms2$mass) # 保持精度
  mod.ms2$intensity<-as.numeric(mod.ms2$intensity)
  unmod.ms2$mass<-as.numeric(unmod.ms2$mass) # 保持精度
  unmod.ms2$intensity<-as.numeric(unmod.ms2$intensity)

  min_intensity=max(max(mod.ms2$intensity),max(unmod.ms2$intensity))/
    min_intensity
  # calculate theoretical b and y ions for given sequence pair
  # input amino acid weight, theoretical weight
  # fu: unmod
  # fm: mod
  fu.mz_b_final = data.frame()
  fm.mz_b_final = data.frame()
  fu.mz_y_final = data.frame()
  fm.mz_y_final = data.frame()

  # ### read unimod file for AA mono weights
  #f.unimod = read.table(f_aaweight, sep = "\t", header = TRUE)
  f.unimod$AA_loc = paste(f.unimod$AA,f.unimod$location,sep="_")
  f.unimod=f.unimod[,c("weight", "AA_loc")]

  ### calculate MW for b/y ions for unmod peptide
  fu.seq = data.frame(AA=unlist(strsplit(f.msms_u$Sequence, split="")),
                      index=1:nchar(f.msms_u$Sequence),stringsAsFactors = F)
  fu.seq$loc="M"
  fu.seq$AA_loc = paste(fu.seq$AA,fu.seq$loc,sep="_")
  fu.mz = merge(fu.seq, f.unimod, by="AA_loc",all.x=T)
  # MW for b ion
  fu.mz_b = fu.mz[with(fu.mz, order(index)),]
  fu.mz_b$accu_weight = cumsum(fu.mz_b$weight)
  # MW for y ion
  fu.mz_y = fu.mz[with(fu.mz, order(-index)),]
  fu.mz_y$accu_weight = cumsum(fu.mz_y$weight)

  ### calculate MW for b ions for mod peptide
  fm.seq = data.frame(AA=unlist(strsplit(f.msms_m$Sequence, split="")),
                      index=1:nchar(f.msms_m$Sequence),stringsAsFactors = FALSE)
  fm.mod_length = nchar(unlist(strsplit(f.msms_m$`Modified sequence`,
                                        split=mod, fixed=TRUE))[1])
  fm.seq$AA[fm.mod_length] = paste(fm.seq$AA[fm.mod_length], mod, sep="")
  fm.seq$loc="M"
  fm.seq$AA_loc = paste(fm.seq$AA,fm.seq$loc,sep="_")
  fm.mz = merge(fm.seq, f.unimod, by="AA_loc",all.x=TRUE)
  # MW for b ion
  fm.mz_b = fm.mz[with(fm.mz, order(index)),]
  fm.mz_b$accu_weight = cumsum(fm.mz_b$weight)
  # MW for y ion
  fm.mz_y = fm.mz[with(fm.mz, order(-index)),]
  fm.mz_y$accu_weight = cumsum(fm.mz_y$weight)

  ##############################################
  # 离子 m/z 计算
  # calculate m/z for b/y ions         ***
  # unmodified 和modified charge 必须相同
  ##############################################
  ms2_charge_u = f.msms_u$Charge - 1
  for ( charge in 1:ms2_charge_u ) {
    #ion_b = (residues + H*z)/z
    fu.mz_b$mz_b = (fu.mz_b$accu_weight + 1.007*charge) / charge
    # ion_y = (residues+h2o+H*z)/z
    fu.mz_y$mz_y = (fu.mz_y$accu_weight + H2O + H*charge) / charge

    fu.mz_b$charge = charge
    fu.mz_y$charge = charge

    if( nrow(fu.mz_b_final)==0 ){
      fu.mz_b_final = fu.mz_b
      fu.mz_y_final = fu.mz_y
    }else{
      fu.mz_b_final = rbind(fu.mz_b_final, fu.mz_b)
      fu.mz_y_final = rbind(fu.mz_y_final, fu.mz_y)
    }
  }
  ms2_charge_m = f.msms_m$Charge - 1
  for (charge in 1:ms2_charge_m) {
    fm.mz_b$mz_b = (fm.mz_b$accu_weight + 1.007*charge) / charge
    fm.mz_y$mz_y = (fm.mz_y$accu_weight + H2O + H*charge) / charge
    fm.mz_b$charge = charge
    fm.mz_y$charge = charge
    if( nrow(fm.mz_b_final)==0 ){
      fm.mz_b_final = fm.mz_b
      fm.mz_y_final = fm.mz_y
    }else{
      fm.mz_b_final = rbind(fm.mz_b_final, fm.mz_b)
      fm.mz_y_final = rbind(fm.mz_y_final, fm.mz_y)
    }
  }


  # fu.psm fm.psm 中的x列 中的推理论值与mod unmod 中的实际值比较
  # 差值
  # ppm=(理论准确分子量(单同位素)-实测分子量(单同位素))/
  # 理论准确分子量(单同位素)*1000000
  # m实际值 M理论值
  # m/z = (1 - 20/1000000)*M/z    ||  m/z = (1 + 20/1000000)* M/z
  # 先通过20ppm计算出数值范围，只匹配范围内最高峰

  # 求出实际质量的范围 weight_min weight_max
  fu.mz_b_final$mz_b_min = (1 - ppm/PPM_denominator) * fu.mz_b_final$mz_b
  fu.mz_b_final$mz_b_max = (1 + ppm/PPM_denominator) * fu.mz_b_final$mz_b
  fu.mz_y_final$mz_y_min = (1 - ppm/PPM_denominator) * fu.mz_y_final$mz_y
  fu.mz_y_final$mz_y_max = (1 + ppm/PPM_denominator) * fu.mz_y_final$mz_y

  fm.mz_b_final$mz_b_min = (1 - ppm/PPM_denominator) * fm.mz_b_final$mz_b
  fm.mz_b_final$mz_b_max = (1 + ppm/PPM_denominator) * fm.mz_b_final$mz_b
  fm.mz_y_final$mz_y_min = (1 - ppm/PPM_denominator) * fm.mz_y_final$mz_y
  fm.mz_y_final$mz_y_max = (1 + ppm/PPM_denominator) * fm.mz_y_final$mz_y

  # unmod: fu b/y ion PSM
  fu.mz_b_psm = data.frame(mz = seq(1:nrow(fu.mz_b_final)))
  fu.mz_y_psm = data.frame(mz = seq(1:nrow(fu.mz_y_final)))
  # 2n^2
  for(i in 1:nrow(unmod.ms2)){
    mz = unmod.ms2$mass[i]
    fu.mz_b_psm$mz = ifelse(mz>fu.mz_b_final$mz_b_min & mz <
                              fu.mz_b_final$mz_b_max, unmod.ms2$intensity[i], 0)
    fu.mz_y_psm$mz = ifelse(mz>fu.mz_y_final$mz_y_min & mz <
                              fu.mz_y_final$mz_y_max, unmod.ms2$intensity[i], 0)
    names(fu.mz_b_psm)[names(fu.mz_b_psm) == "mz"] <- mz
    names(fu.mz_y_psm)[names(fu.mz_y_psm) == "mz"] <- mz
  }

  fu.mz_b_psm_colsum = colSums(fu.mz_b_psm)
  fu.mz_y_psm_colsum = colSums(fu.mz_y_psm)
  fu.mz_b_psm_colsum_index = which(fu.mz_b_psm_colsum >0 )
  fu.mz_y_psm_colsum_index = which(fu.mz_y_psm_colsum >0 )
  fu.mz_b_psm = fu.mz_b_psm[fu.mz_b_psm_colsum_index]
  fu.mz_y_psm = fu.mz_y_psm[fu.mz_y_psm_colsum_index]

  # leave the MZ with the largest intensity if multiple MZs match the same ion
  # find the max per row and retain the columns with the max
  rowmax <- apply(fu.mz_b_psm, 1, max)
  index_b = vector()
  for(i in 1:length(rowmax)){
    if(rowmax[i]>0){
      index = which(fu.mz_b_psm[i,] == rowmax[i] )
      index_b = c(index_b, index)
    }
  }
  fu.mz_b_psm = fu.mz_b_psm[index_b]

  rowmax <- apply(fu.mz_y_psm, 1, max)
  index_y = vector()

  for(i in 1:length(rowmax)){
    if(rowmax[i]>0){
      index = which(fu.mz_y_psm[i,] == rowmax[i] )
      index_y = c(index_y, index)
    }
  }
  fu.mz_y_psm = fu.mz_y_psm[index_y]

  # 去除多余列 不影响后续
  l = vector()
  for (i in 1:length(fu.mz_b_psm)) {
    if(lengths(strsplit(names(fu.mz_b_psm[i]),'\\.'))>2){
      #print(names(fu.mz_b_psm[i]))
      #fu.mz_b_psm = fu.mz_b_psm[,-i]
      l <- c(l,i)
    }
  }
  l = sort(l,decreasing = TRUE)
  for (i in l) {
    fu.mz_b_psm = fu.mz_b_psm[,-i]
  }

  l = vector()
  for (i in 1:length(fu.mz_y_psm)) {
    if(lengths(strsplit(names(fu.mz_y_psm[i]),'\\.'))>2){
      #print(names(fu.mz_y_psm[i]))
      l <- c(l,i)
    }
  }
  l = sort(l,decreasing = TRUE)
  for (i in l) {
    fu.mz_y_psm = fu.mz_y_psm[,-i]
  }


  ### combine
  fu.mz_b_final = cbind(fu.mz_b_final, fu.mz_b_psm)
  fu.mz_y_final = cbind(fu.mz_y_final, fu.mz_y_psm)

  #####################################
  #maxquant  unmod from f.msms_old
  #
  #y1;y3;y4;y5;y6;y7;y8;y9;
  #fu.mz_y_final 匹配9个；正常
  #b2;b3;b4
  #fu.mz_b_final 匹配3个；正常

  ##################################
  # mod: fm b/y ion PSM
  fm.mz_b_psm = data.frame(mz = seq(1:nrow(fm.mz_b_final)))
  fm.mz_y_psm = data.frame(mz = seq(1:nrow(fm.mz_y_final)))
  for(i in 1:nrow(mod.ms2)){
    mz = mod.ms2$mass[i]
    fm.mz_b_psm$mz = ifelse(mz>fm.mz_b_final$mz_b_min & mz <
                              fm.mz_b_final$mz_b_max, mod.ms2$intensity[i], 0)
    fm.mz_y_psm$mz = ifelse(mz>fm.mz_y_final$mz_y_min & mz <
                              fm.mz_y_final$mz_y_max, mod.ms2$intensity[i], 0)
    names(fm.mz_b_psm)[names(fm.mz_b_psm) == "mz"] <- mz
    names(fm.mz_y_psm)[names(fm.mz_y_psm) == "mz"] <- mz
  }

  fm.mz_b_psm_colsum = colSums(fm.mz_b_psm)
  fm.mz_y_psm_colsum = colSums(fm.mz_y_psm)
  fm.mz_b_psm_colsum_index = which(fm.mz_b_psm_colsum >0 )
  fm.mz_y_psm_colsum_index = which(fm.mz_y_psm_colsum >0 )
  fm.mz_b_psm = fm.mz_b_psm[fm.mz_b_psm_colsum_index]
  fm.mz_y_psm = fm.mz_y_psm[fm.mz_y_psm_colsum_index]

  # leave the MZ with the largest intensity if multiple MZs match the same ion
  # find the max per row and retain the columns with the max
  rowmax <- apply(fm.mz_b_psm, 1, max)
  index_b = vector()
  for(i in 1:length(rowmax)){
    if(rowmax[i]>0){
      index = which(fm.mz_b_psm[i,] == rowmax[i] )
      index_b = c(index_b, index)
    }
  }
  fm.mz_b_psm = fm.mz_b_psm[index_b]

  rowmax <- apply(fm.mz_y_psm, 1, max)
  index_y = vector()
  for(i in 1:length(rowmax)){
    if(rowmax[i]>0){
      index = which(fm.mz_y_psm[i,] == rowmax[i] )
      index_y = c(index_y, index)
    }
  }
  fm.mz_y_psm = fm.mz_y_psm[index_y]

  # 去除多余列
  l = vector()
  for (i in 1:length(fm.mz_b_psm)) {
    if(lengths(strsplit(names(fm.mz_b_psm[i]),'\\.'))>2){
      #print(names(fm.mz_b_psm[i]))
      l <- c(l,i)
    }
  }
  l = sort(l,decreasing = TRUE)
  for (i in l) {
    fm.mz_b_psm = fm.mz_b_psm[,-i]
  }

  l = vector()
  for (i in 1:length(fm.mz_y_psm)) {
    if(lengths(strsplit(names(fm.mz_y_psm[i]),'\\.'))>2){
      #print(names(fm.mz_y_psm[i]))
      l <- c(l,i)
    }
  }
  l = sort(l,decreasing = TRUE)
  for (i in l) {
    fm.mz_y_psm = fm.mz_y_psm[,-i]
  }




  ### combine
  fm.mz_b_final = cbind(fm.mz_b_final, fm.mz_b_psm)
  fm.mz_y_final = cbind(fm.mz_y_final, fm.mz_y_psm)


  # intensity太低的m/z就不显示，也不和b/y ion匹配
  # intensity太低的缺省定义是最高峰的1%。
  ## 去除掉峰度太低的匹配的b y离子
  # fu.mz_b_final
  l = vector()
  for (j in 11:length(fu.mz_b_final)) {
    for (i in 1:nrow(fu.mz_b_final)) {
      if(fu.mz_b_final[i,j] < min_intensity && fu.mz_b_final[i,j] >0){
        l <- c(l,j)
      }
    }
  }
  l = sort(l,decreasing = TRUE)
  for (i in l) {
    fu.mz_b_final = fu.mz_b_final[,-i]
  }
  # fu.mz_y_final
  l = vector()
  for (j in 11:length(fu.mz_y_final)) {
    for (i in 1:nrow(fu.mz_y_final)) {
      if(fu.mz_y_final[i,j] < min_intensity && fu.mz_y_final[i,j] >0){
        l <- c(l,j)
      }
    }
  }
  l = sort(l,decreasing = TRUE)
  for (i in l) {
    fu.mz_y_final = fu.mz_y_final[,-i]
  }

  # fm.mz_b_final
  l = vector()
  for (j in 11:length(fm.mz_b_final)) {
    for (i in 1:nrow(fm.mz_b_final)) {
      if(fm.mz_b_final[i,j] < min_intensity && fm.mz_b_final[i,j] >0){
        l <- c(l,j)
      }
    }
  }
  l = sort(l,decreasing = TRUE)
  for (i in l) {
    fm.mz_b_final = fm.mz_b_final[,-i]
  }

  # fm.mz_y_final
  l = vector()
  for (j in 11:length(fm.mz_y_final)) {
    for (i in 1:nrow(fm.mz_y_final)) {
      if(fm.mz_y_final[i,j] < min_intensity && fm.mz_y_final[i,j] >0){
        l <- c(l,j)
      }
    }
  }
  l = sort(l,decreasing = TRUE)
  for (i in l) {
    fm.mz_y_final = fm.mz_y_final[,-i]
  }
  ## 去除掉峰度太低的ms2
  l = vector()
  for (i in 1:nrow(mod.ms2)) {
    if(mod.ms2$intensity[i]<min_intensity){
      l <- c(l,i)
    }
  }
  l = sort(l,decreasing = TRUE)
  for (i in l) {
    mod.ms2 = mod.ms2[-i,]
  }

  l = vector()
  for (i in 1:nrow(unmod.ms2)) {
    if(unmod.ms2$intensity[i]<min_intensity){
      l <- c(l,i)
    }
  }
  l = sort(l,decreasing = TRUE)
  for (i in l) {
    unmod.ms2 = unmod.ms2[-i,]
  }

  # msms中的数值不需要画出
  # 先进行mod.ms2和unmod.ms2等比缩放以保持上下画幅一致，
  # 并将unmod.ms2的intensity变为负数
  if(max(unmod.ms2$intensity)>max(mod.ms2$intensity)){
    ratio = max(unmod.ms2$intensity)/max(mod.ms2$intensity)
    unmod.ms2$intensity_adj = unmod.ms2$intensity/ratio * -1  #调整为负值
    mod.ms2$intensity_adj = mod.ms2$intensity
    fu.mz_b_final[,11:length(fu.mz_b_final)] =
      fu.mz_b_final[,11:length(fu.mz_b_final)]/ratio *-1
    fu.mz_y_final[,11:length(fu.mz_y_final)] =
      fu.mz_y_final[,11:length(fu.mz_y_final)]/ratio *-1
  }else{
    ratio = max(mod.ms2$intensity)/max(unmod.ms2$intensity)
    mod.ms2$intensity_adj = mod.ms2$intensity/ratio
    unmod.ms2$intensity_adj = unmod.ms2$intensity * -1
    fm.mz_b_final[,11:length(fm.mz_b_final)] =
      fm.mz_b_final[,11:length(fm.mz_b_final)]/ratio
    fm.mz_y_final[,11:length(fm.mz_y_final)] =
      fm.mz_y_final[,11:length(fm.mz_y_final)]/ratio
    fu.mz_b_final[,11:length(fu.mz_b_final)] =
      fu.mz_b_final[,11:length(fu.mz_b_final)] * -1
    fu.mz_y_final[,11:length(fu.mz_y_final)] =
      fu.mz_y_final[,11:length(fu.mz_y_final)] * -1
  }
  f.psm = rbind(mod.ms2, unmod.ms2)

  # max of MS1 masses
  n.ms1.mass = max(max(f.psm$mass),max(fm.mz_b_final$mz_b),
                   max(fm.mz_y_final$mz_y),max(fu.mz_b_final$mz_b),
                   max(fu.mz_y_final$mz_y))
  # max of MS2 mod
  n.ms2.max.intensity = max(f.psm$intensity_adj)
  # min of MS2 unmod
  n.ms2.min.intensity = min(f.psm$intensity_adj)

  ############################################
  # 画图
  options(scipen=22)
  #outer margin
  par(oma=c(0,4,0,4))
  # margin
  #par(mar=c(4,2,2,2))
  par(xaxs = "i", yaxs = "i")
  graphics::plot(f.psm$mass, f.psm$intensity_adj, type = "h",las=1,xlab="m/z",
                 ylab="",yaxt="n",xaxt="n", xlim=c(0, n.ms1.mass),
                 ylim=c(n.ms2.min.intensity*1.7, n.ms2.max.intensity*1.7))

  options(scipen=200)

  abline(h = 0)
  # label x axis
  axis(side = 1,at=seq(0, n.ms1.mass, by=200))

  # label y axis (positive)
  axis_y_positive = seq(0,signif(n.ms2.max.intensity,digits=2),
                        signif(n.ms2.max.intensity,digits=2)/4)
  axis_y_positive_real = seq(0,n.ms2.max.intensity,n.ms2.max.intensity/4)
  axis(side = 2,las = 1,at = axis_y_positive,labels =
         format(axis_y_positive,scientific=F))
  axis(side = 4,las = 1,at = axis_y_positive_real,
       labels = c("0%","25%","50%","75%","100%"))

  # label y axis (negative)
  axis_y_negative = seq(0,signif(n.ms2.min.intensity,digits=2),
                        signif(n.ms2.min.intensity,digits=2)/4)
  axis_y_negative2 = seq(0,n.ms2.min.intensity,n.ms2.min.intensity/4)
  axis_y_negative_real = seq(0,signif(max(unmod.ms2$intensity),digits=2),
                             signif(max(unmod.ms2$intensity),digits=2)/4)
  axis_y_negative = axis_y_negative[2:length(axis_y_negative)]
  axis_y_negative_real = axis_y_negative_real[2:length(axis_y_negative_real)]
  axis_y_negative2 = axis_y_negative2[2:length(axis_y_negative2)]
  axis(side = 2,las = 1,at = axis_y_negative,labels =
         format(axis_y_negative_real,scientific=F))
  axis(side = 4,las = 1,at = axis_y_negative2,labels =
         c("25%","50%","75%","100%"))

  # 画出匹配的b y 离子，并在柱顶标注
  # y离子 index 顺序调整
  for (i in 1:nrow(fu.mz_y_final)) {
    fu.mz_y_final$index[i] = nrow(fu.mz_y_final)/
      max(fu.mz_y_final$charge) -fu.mz_y_final$index[i] +1
  }
  for (i in 1:nrow(fm.mz_y_final)) {
    fm.mz_y_final$index[i] = nrow(fm.mz_y_final)/
      max(fm.mz_y_final$charge) -fm.mz_y_final$index[i] +1
  }
  # 先遍历列再遍历行
  for (i in 11:length(fm.mz_b_final)) {
    for (j in 1:nrow(fm.mz_b_final)) {
      if(fm.mz_b_final[j,i] != 0){
        if(fm.mz_b_final$charge[j]==1){
          lines(as.numeric(names(fm.mz_b_final[i])),
                fm.mz_b_final[j,i], type = "h",las=1, col="green")
          text(as.numeric(names(fm.mz_b_final[i])),fm.mz_b_final[j,i],
               paste('b',fm.mz_b_final$index[j],'+',sep = ''),
               cex = cex,col="green",pos = 3,srt=srt)
        }else if(fm.mz_b_final$charge[j]==2){
          lines(as.numeric(names(fm.mz_b_final[i])), fm.mz_b_final[j,i],
                type = "h",las=1, col="green")
          text(as.numeric(names(fm.mz_b_final[i])),fm.mz_b_final[j,i],
               paste('b',fm.mz_b_final$index[j],'++',sep = ''),cex = cex,
               col="green",pos = 3,srt=srt)
        }else if(fm.mz_b_final$charge[j]==3){
          lines(as.numeric(names(fm.mz_b_final[i])), fm.mz_b_final[j,i],
                type = "h",las=1, col="green")
          text(as.numeric(names(fm.mz_b_final[i])),fm.mz_b_final[j,i],
               paste('b',fm.mz_b_final$index[j],'+++',sep = ''),
               cex = cex,col="green",pos = 3,srt=srt)
        }
      }
    }
  }
  for (i in 11:length(fm.mz_y_final)) {
    for (j in 1:nrow(fm.mz_y_final)) {
      if(fm.mz_y_final[j,i] != 0){
        if(fm.mz_y_final$charge[j]==1){
          lines(as.numeric(names(fm.mz_y_final[i])), fm.mz_y_final[j,i],
                type = "h",las=1, col="orange")
          text(as.numeric(names(fm.mz_y_final[i])),fm.mz_y_final[j,i],
               paste('y',fm.mz_y_final$index[j],'+ ',sep = ''),
               cex = cex, col="orange",pos = 3,srt=srt)
        }else if(fm.mz_y_final$charge[j]==2){
          lines(as.numeric(names(fm.mz_y_final[i])), fm.mz_y_final[j,i],
                type = "h",las=1, col="orange")
          text(as.numeric(names(fm.mz_y_final[i])),fm.mz_y_final[j,i],
               paste('y',fm.mz_y_final$index[j],'++',sep = ''),
               cex = cex, col="orange",pos = 3,srt=srt)
        }else if(fm.mz_y_final$charge[j]==3){
          lines(as.numeric(names(fm.mz_y_final[i])), fm.mz_y_final[j,i],
                type = "h",las=1, col="orange")
          text(as.numeric(names(fm.mz_y_final[i])),fm.mz_y_final[j,i],
               paste('y',fm.mz_y_final$index[j],'+++',sep = ''),
               cex = cex, col="orange",pos = 3,srt=srt)
        }
      }
    }
  }
  for (i in 11:length(fu.mz_b_final)) {
    for (j in 1:nrow(fu.mz_b_final)) {
      if(fu.mz_b_final[j,i] != 0){
        if(fu.mz_b_final$charge[j]==1){
          lines(as.numeric(names(fu.mz_b_final[i])), fu.mz_b_final[j,i],
                type = "h",las=1, col="green")
          text(as.numeric(names(fu.mz_b_final[i])),fu.mz_b_final[j,i],
               paste('b',fu.mz_b_final$index[j],'+',sep = ''),
               cex = cex, col="green",pos = 1,srt=srt)
        }else if(fu.mz_b_final$charge[j]==2){
          lines(as.numeric(names(fu.mz_b_final[i])), fu.mz_b_final[j,i],
                type = "h",las=1, col="green")
          text(as.numeric(names(fu.mz_b_final[i])),fu.mz_b_final[j,i],
               paste('b',fu.mz_b_final$index[j],'++',sep = ''),
               cex = cex, col="green",pos = 1,srt=srt)
        }else if(fu.mz_b_final$charge[j]==3){
          lines(as.numeric(names(fu.mz_b_final[i])), fu.mz_b_final[j,i],
                type = "h",las=1, col="green")
          text(as.numeric(names(fu.mz_b_final[i])),fu.mz_b_final[j,i],
               paste('b',fu.mz_b_final$index[j],'+++',sep = ''),
               cex = cex, col="green",pos = 1,srt=srt)
        }
      }
    }
  }
  for (i in 11:length(fu.mz_y_final)) {
    for (j in 1:nrow(fu.mz_y_final)) {
      if(fu.mz_y_final[j,i] != 0){
        if(fu.mz_y_final$charge[j]==1){
          lines(as.numeric(names(fu.mz_y_final[i])), fu.mz_y_final[j,i],
                type = "h",las=1, col="orange")
          text(as.numeric(names(fu.mz_y_final[i])),fu.mz_y_final[j,i],
               paste('y',fu.mz_y_final$index[j],'+',sep = ''),
               cex = cex, col="orange",pos = 1,srt=srt)
        }else if(fu.mz_y_final$charge[j]==2){
          lines(as.numeric(names(fu.mz_y_final[i])), fu.mz_y_final[j,i],
                type = "h",las=1, col="orange")
          text(as.numeric(names(fu.mz_y_final[i])),fu.mz_y_final[j,i],
               paste('y',fu.mz_y_final$index[j],'++',sep = ''),
               cex = cex, col="orange",pos = 1,srt=srt)
        }else if(fu.mz_y_final$charge[j]==3){
          lines(as.numeric(names(fu.mz_y_final[i])), fu.mz_y_final[j,i],
                type = "h",las=1, col="orange")
          text(as.numeric(names(fu.mz_y_final[i])),fu.mz_y_final[j,i],
               paste('y',fu.mz_y_final$index[j],'+++',sep = ''),
               cex = cex, col="orange",pos = 1,srt=srt)
        }
      }
    }
  }




  # 对于charge为2 的子离子，先删除处理
  # proteomicsdb只对charge为1的进行氨基酸的匹配
  fm.mz_b_final = fm.mz_b_final[which(fm.mz_b_final$charge==1),]
  fm.mz_y_final = fm.mz_y_final[which(fm.mz_y_final$charge==1),]
  fu.mz_b_final = fu.mz_b_final[which(fu.mz_b_final$charge==1),]
  fu.mz_y_final = fu.mz_y_final[which(fu.mz_y_final$charge==1),]


  tag1 <- array(0)
  for(k in 1:nrow(fm.mz_b_final)){
    tag1[k] <- 0
  }
  tag2 <- array(0)
  for(k in 1:nrow(fu.mz_b_final)){
    tag2[k] <- 0
  }




  for(i in 1:nrow(fm.mz_b_final)){  # 140%
    par(new = TRUE)
    if(i<nrow(fm.mz_b_final)){

      segments(x0=fm.mz_b_final$mz_b[i], y0=n.ms2.max.intensity*1.35,
               x1=fm.mz_b_final$mz_b[i], y1=n.ms2.max.intensity*1.45)
      segments(x0=fu.mz_b_final$mz_b[i], y0=n.ms2.min.intensity*1.35,
               x1=fu.mz_b_final$mz_b[i], y1=n.ms2.min.intensity*1.45)


      # 最后一个手动赋值
      segments(fm.mz_b_final$mz_b[nrow(fm.mz_b_final)-1],
               n.ms2.max.intensity*1.4, fm.mz_b_final$mz_b[
                 nrow(fm.mz_b_final)-1]+0.25*(n.ms1.mass-fm.mz_b_final$mz_b[
                   nrow(fm.mz_b_final)-1]), n.ms2.max.intensity*1.4)
      segments(fm.mz_b_final$mz_b[nrow(fm.mz_b_final)-1]+0.75*(
        n.ms1.mass-fm.mz_b_final$mz_b[nrow(fm.mz_b_final)-1]),
        n.ms2.max.intensity*1.4,n.ms1.mass,n.ms2.max.intensity*1.4)

      segments(fu.mz_b_final$mz_b[nrow(fu.mz_b_final)-1],
               n.ms2.min.intensity*1.4, fu.mz_b_final$mz_b[
                 nrow(fu.mz_b_final)-1]+0.25*(n.ms1.mass-fu.mz_b_final$mz_b[
                   nrow(fu.mz_b_final)-1]), n.ms2.min.intensity*1.4)
      segments(fu.mz_b_final$mz_b[nrow(fu.mz_b_final)-1]+0.75*(
        n.ms1.mass-fu.mz_b_final$mz_b[nrow(fu.mz_b_final)-1]),
        n.ms2.min.intensity*1.4,n.ms1.mass,n.ms2.min.intensity*1.4)


      # 第一次添加黑线25%  75%
      if(i==1){
        segments(0,n.ms2.max.intensity*1.4,0.25*fm.mz_b_final$mz_b[i],
                 n.ms2.max.intensity*1.4)
        segments(0.75*fm.mz_b_final$mz_b[i],n.ms2.max.intensity*1.4,
                 fm.mz_b_final$mz_b[i],n.ms2.max.intensity*1.4)

        segments(0,n.ms2.min.intensity*1.4,0.25*fu.mz_b_final$mz_b[i],
                 n.ms2.min.intensity*1.4)
        segments(0.75*fu.mz_b_final$mz_b[i],n.ms2.min.intensity*1.4,
                 fu.mz_b_final$mz_b[i],n.ms2.min.intensity*1.4)

      }else{
        segments(fm.mz_b_final$mz_b[i-1], n.ms2.max.intensity*1.4,
                 fm.mz_b_final$mz_b[i-1]+0.25*(
                   fm.mz_b_final$mz_b[i]-fm.mz_b_final$mz_b[i-1]),
                 n.ms2.max.intensity*1.4)
        segments(fm.mz_b_final$mz_b[i-1]+0.75*(
          fm.mz_b_final$mz_b[i]-fm.mz_b_final$mz_b[i-1]),
        n.ms2.max.intensity*1.4, fm.mz_b_final$mz_b[i], n.ms2.max.intensity*1.4)

        segments(fu.mz_b_final$mz_b[i-1], n.ms2.min.intensity*1.4,
                 fu.mz_b_final$mz_b[i-1]+0.25*(
                   fu.mz_b_final$mz_b[i]-fu.mz_b_final$mz_b[i-1]),
                 n.ms2.min.intensity*1.4)
        segments(fu.mz_b_final$mz_b[i-1]+0.75*(
          fu.mz_b_final$mz_b[i]-fu.mz_b_final$mz_b[i-1]),
        n.ms2.min.intensity*1.4, fu.mz_b_final$mz_b[i], n.ms2.min.intensity*1.4)

      }




    }
    if(i==1){
      shadowtext(x=fm.mz_b_final$mz_b[i]/2, y=n.ms2.max.intensity*1.4,
                 labels=fm.mz_b_final$AA[i], col="black", bg="white")
      shadowtext(x=fu.mz_b_final$mz_b[i]/2, y=n.ms2.min.intensity*1.4,
                 labels=fu.mz_b_final$AA[i], col="black", bg="white")

    }else{
      shadowtext(x=fm.mz_b_final$mz_b[i-1]+(
        fm.mz_b_final$mz_b[i]-fm.mz_b_final$mz_b[i-1])/2,
        y=n.ms2.max.intensity*1.4, labels=fm.mz_b_final$AA[i],
        col="black", bg="white")
      shadowtext(x=fu.mz_b_final$mz_b[i-1]+(
        fu.mz_b_final$mz_b[i]-fu.mz_b_final$mz_b[i-1])/2,
        y=n.ms2.min.intensity*1.4, labels=fu.mz_b_final$AA[i],
        col="black", bg="white")
    }

  }

  #140%匹配标红
  for(i in 1:nrow(fm.mz_b_final)){
    # 为了防止 需要TRUE/FALSE值的地方不可以用缺少值错误，从下标2开始遍历

    #匹配上色
    for (j in 11:length(fm.mz_b_final)) {
      if(fm.mz_b_final[i,j] != 0){
        segments(x0=fm.mz_b_final$mz_b[i], y0=n.ms2.max.intensity*1.35,
                 x1=fm.mz_b_final$mz_b[i], y1=n.ms2.max.intensity*1.45,
                 col = "green",lwd = 3)
        tag1[i] = fm.mz_b_final$mz_b[i]
      }
    }
    for (j in 11:length(fu.mz_b_final)) {
      if(fu.mz_b_final[i,j] != 0){
        segments(x0=fu.mz_b_final$mz_b[i], y0=n.ms2.min.intensity*1.35,
                 x1=fu.mz_b_final$mz_b[i], y1=n.ms2.min.intensity*1.45,
                 col = "green",lwd = 3)
        tag2[i] = fu.mz_b_final$mz_b[i]
      }
    }
    for (j in 2:length(tag1)) {
      if(tag1[j] > 0 && tag1[j-1] != 0 && fm.mz_b_final$mz_b[i-1] == tag1[j-1]
         && fm.mz_b_final$mz_b[i]== tag1[j]){
        segments(tag1[i-1], n.ms2.max.intensity*1.4, tag1[i-1]+0.25*(
          tag1[i]-tag1[i-1]), n.ms2.max.intensity*1.4,col = "green")
        segments(tag1[i-1]+0.75*(tag1[i]-tag1[i-1]), n.ms2.max.intensity*1.4,
                 tag1[i], n.ms2.max.intensity*1.4,col = "green")

        shadowtext(x=fm.mz_b_final$mz_b[i-1]+(
          fm.mz_b_final$mz_b[i]-fm.mz_b_final$mz_b[i-1])/2,
          y=n.ms2.max.intensity*1.4, labels=fm.mz_b_final$AA[i],
          col="green", bg="white")
      }
    }
    #140%匹配标红
    for (j in 2:length(tag2)) {
      if(tag2[j] > 0 && tag2[j-1] != 0 && fu.mz_b_final$mz_b[i-1] == tag2[j-1]
         && fu.mz_b_final$mz_b[i]== tag2[j]){
        segments(tag2[i-1], n.ms2.min.intensity*1.4, tag2[i-1]+
                0.25*(tag2[i]-tag2[i-1]), n.ms2.min.intensity*1.4,col = "green")
        segments(tag2[i-1]+0.75*(tag2[i]-tag2[i-1]), n.ms2.min.intensity*1.4,
                 tag2[i], n.ms2.min.intensity*1.4,col = "green")
        shadowtext(x=fu.mz_b_final$mz_b[i-1]+(
          fu.mz_b_final$mz_b[i]-fu.mz_b_final$mz_b[i-1])/2,
          y=n.ms2.min.intensity*1.4, labels=fu.mz_b_final$AA[i],
          col="green", bg="white")
      }
    }
  }


  tag1 <- array(0)
  for(k in 1:nrow(fm.mz_y_final)){
    tag1[k] <- 0
  }
  tag2 <- array(0)
  for(k in 1:nrow(fu.mz_y_final)){
    tag2[k] <- 0
  }

  for(i in 1:nrow(fm.mz_y_final)){  # 125%
    par(new = TRUE)
    if(i<nrow(fm.mz_y_final)){
      segments(x0=fm.mz_y_final$mz_y[i], y0=n.ms2.max.intensity*1.2,
               x1=fm.mz_y_final$mz_y[i], y1=n.ms2.max.intensity*1.3)
      segments(x0=fu.mz_y_final$mz_y[i], y0=n.ms2.min.intensity*1.2,
               x1=fu.mz_y_final$mz_y[i], y1=n.ms2.min.intensity*1.3)

      segments(fm.mz_y_final$mz_y[nrow(fm.mz_y_final)-1],
            n.ms2.max.intensity*1.25, fm.mz_y_final$mz_y[nrow(fm.mz_y_final)-1]
               +0.25*(n.ms1.mass-fm.mz_y_final$mz_y[nrow(fm.mz_y_final)-1]),
            n.ms2.max.intensity*1.25)
      segments(fm.mz_y_final$mz_y[nrow(fm.mz_y_final)-1]+0.75*(
        n.ms1.mass-fm.mz_y_final$mz_y[nrow(fm.mz_y_final)-1]),
        n.ms2.max.intensity*1.25,n.ms1.mass,n.ms2.max.intensity*1.25)

      segments(fu.mz_y_final$mz_y[nrow(fu.mz_y_final)-1],
            n.ms2.min.intensity*1.25, fu.mz_y_final$mz_y[nrow(fu.mz_y_final)-1]
               +0.25*(n.ms1.mass-fu.mz_y_final$mz_y[nrow(fu.mz_y_final)-1]),
            n.ms2.min.intensity*1.25)
      segments(fu.mz_y_final$mz_y[nrow(fu.mz_y_final)-1]+0.75*(
        n.ms1.mass-fu.mz_y_final$mz_y[nrow(fu.mz_y_final)-1]),
        n.ms2.min.intensity*1.25,n.ms1.mass,n.ms2.min.intensity*1.25)


      # 25% 75% 黑线
      if(i==1){
        segments(0,n.ms2.max.intensity*1.25,fm.mz_y_final$mz_y[i]*0.25,
                 n.ms2.max.intensity*1.25)
        segments(fm.mz_y_final$mz_y[i]*0.75,n.ms2.max.intensity*1.25,
                 fm.mz_y_final$mz_y[i],n.ms2.max.intensity*1.25)

        segments(0,n.ms2.min.intensity*1.25,fu.mz_y_final$mz_y[i]*0.25,
                 n.ms2.min.intensity*1.25)
        segments(fu.mz_y_final$mz_y[i]*0.75,n.ms2.min.intensity*1.25,
                 fu.mz_y_final$mz_y[i],n.ms2.min.intensity*1.25)
      }else{
        segments(fm.mz_y_final$mz_y[i-1], n.ms2.max.intensity*1.25,
                 fm.mz_y_final$mz_y[i-1]+(
                   fm.mz_y_final$mz_y[i]-fm.mz_y_final$mz_y[i-1])*0.25,
                 n.ms2.max.intensity*1.25)
        segments(fm.mz_y_final$mz_y[i-1]+(
          fm.mz_y_final$mz_y[i]- fm.mz_y_final$mz_y[i-1])*0.75,
      n.ms2.max.intensity*1.25, fm.mz_y_final$mz_y[i], n.ms2.max.intensity*1.25)

        segments(fu.mz_y_final$mz_y[i-1], n.ms2.min.intensity*1.25,
                 fu.mz_y_final$mz_y[i-1]+(
                   fu.mz_y_final$mz_y[i]-fu.mz_y_final$mz_y[i-1])*0.25,
                 n.ms2.min.intensity*1.25)
        segments(fu.mz_y_final$mz_y[i-1]+(
          fu.mz_y_final$mz_y[i]-fu.mz_y_final$mz_y[i-1])*0.75,
      n.ms2.min.intensity*1.25, fu.mz_y_final$mz_y[i], n.ms2.min.intensity*1.25)

      }



    }
    if(i==1){
      shadowtext(x=fm.mz_y_final$mz_y[i]/2, y=n.ms2.max.intensity*1.25,
                 labels=fm.mz_y_final$AA[i], col="black", bg="white")
      shadowtext(x=fu.mz_y_final$mz_y[i]/2, y=n.ms2.min.intensity*1.25,
                 labels=fu.mz_y_final$AA[i], col="black", bg="white")
    }else{
      shadowtext(x=fm.mz_y_final$mz_y[i-1]+(
        fm.mz_y_final$mz_y[i]-fm.mz_y_final$mz_y[i-1])/2,
        y=n.ms2.max.intensity*1.25, labels=fm.mz_y_final$AA[i], col="black",
        bg="white")
      shadowtext(x=fu.mz_y_final$mz_y[i-1]+(
        fu.mz_y_final$mz_y[i]-fu.mz_y_final$mz_y[i-1])/2,
        y=n.ms2.min.intensity*1.25, labels=fu.mz_y_final$AA[i], col="black",
        bg="white")
    }

  }


  for(i in 1:nrow(fm.mz_y_final)){

    # 匹配上色
    for (j in 11:length(fm.mz_y_final)) {
      if(fm.mz_y_final[i,j] != 0){
        segments(x0=fm.mz_y_final$mz_y[i], y0=n.ms2.max.intensity*1.2,
                 x1=fm.mz_y_final$mz_y[i], y1=n.ms2.max.intensity*1.3,
                 col = "orange",lwd = 3)
        tag1[i] = fm.mz_y_final$mz_y[i]
      }
    }
    for (j in 11:length(fu.mz_y_final)) {
      if(fu.mz_y_final[i,j] != 0){
        segments(x0=fu.mz_y_final$mz_y[i], y0=n.ms2.min.intensity*1.2,
                 x1=fu.mz_y_final$mz_y[i], y1=n.ms2.min.intensity*1.3,
                 col = "orange",lwd = 3)
        tag2[i] = fu.mz_y_final$mz_y[i]
      }
    }
    for (j in 2:length(tag1)) {
      if(tag1[j] > 0 && tag1[j-1] != 0 && fm.mz_y_final$mz_y[i-1] == tag1[j-1]
         && fm.mz_y_final$mz_y[i]== tag1[j]){
        segments(tag1[i-1], n.ms2.max.intensity*1.25,
                 tag1[i-1]+0.25*(tag1[i]-tag1[i-1]), n.ms2.max.intensity*1.25,
                 col = "orange")
        segments(tag1[i-1]+0.75*(tag1[i]-tag1[i-1]), n.ms2.max.intensity*1.25,
                 tag1[i], n.ms2.max.intensity*1.25,col = "orange")

        shadowtext(x=fm.mz_y_final$mz_y[i-1]+(
          fm.mz_y_final$mz_y[i]-fm.mz_y_final$mz_y[i-1])/2,
          y=n.ms2.max.intensity*1.25, labels=fm.mz_y_final$AA[i], col="orange",
          bg="white")
      }
    }
    for (j2 in 2:length(tag2)) {
      if(tag2[j2] > 0 && tag2[j2-1] != 0 && fu.mz_y_final$mz_y[i-1] ==tag2[j2-1]
         && fu.mz_y_final$mz_y[i]== tag2[j2]){
        segments(tag2[i-1], n.ms2.min.intensity*1.25, tag2[i-1]+
              0.25*(tag2[i]-tag2[i-1]), n.ms2.min.intensity*1.25,col = "orange")
        segments(tag2[i-1]+0.75*(tag2[i]-tag2[i-1]), n.ms2.min.intensity*1.25,
                 tag2[i], n.ms2.min.intensity*1.25,col = "orange")
        shadowtext(x=fu.mz_y_final$mz_y[i-1]+(
          fu.mz_y_final$mz_y[i]-fu.mz_y_final$mz_y[i-1])/2,
          y=n.ms2.min.intensity*1.25, labels=fu.mz_y_final$AA[i], col="orange",
          bg="white")
      }
    }
  }

  rawfile=unlist(strsplit(f.msms_m[1,1], "\\\\|/|;|=", fixed = FALSE))[
    length(unlist(strsplit(f.msms_m[1,1], "\\\\|/|;|=", fixed = FALSE)))]
  shadowtext(n.ms1.mass/2,n.ms2.max.intensity*1.6,labels =
               paste("File:",rawfile,"    Scan:",f.msms_m$`Scan number`,
                     "    RT(min):",round(f.msms_m$`Retention time`,1),
                     "    m/z:",f.msms_m$`m/z`,
                     paste("    Charge:",f.msms_m$Charge,"+",sep = ""),
                     "    Gene name(s):",f.msms_m$`Gene Names`),col = "black",
             bg = "white")
  shadowtext(n.ms1.mass/2,n.ms2.min.intensity*1.6,labels =
               paste("File:",rawfile,"    Scan:",f.msms_u$`Scan number`,
                     "    RT(min):",round(f.msms_u$`Retention time`,1),
                     "    m/z:",f.msms_u$`m/z`,
                     paste("    Charge:",f.msms_u$Charge,"+",sep = ""),
                     "    Gene name(s):",f.msms_u$`Gene Names`),col = "black",
             bg = "white")


  dev.off()
}
