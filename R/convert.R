#' @title convert
#' @author Oyvind Bleka (in collaboration with Rebecca Just)
#' @description Sequence to LUS format converter
#' @details Details This script enables users to provide their own UAS sequence data and convert the sequences directly to different formats. The formats supported are RU (Repeat Unit/CE), LUS (Longest Uninterupted Stretch), LUS+ (extended version of LUS), or STRseq (The short format provided by STRseq)
#' The lookup table was created by Rebecca Just. The converter  
#' @param dat A data table with format [SampleName,Marker,Allele,Height] where Alleles are sequences. Height column is optional (indicating Evid if present and Ref if not present).
#' @param lookupf The destination/path of the UAS Lookup file (use build-in by default). Must be an xlsx file with same format.
#' @param type The convertion type (RU, LUS, LUS+)
#' @param toFile Boolean of whether the converted data table should be saved to workspace 
#' @param outf The file name of the converted data table
#' @param printMissing A boolean of whether alleles not found in the look-up table should be printed to a seperate file
#' @param missf The file name of the missing alleles (not found in look-up table)
#' @param LUSsymb The separator format to use.
#' @return Returning a table with format [SampleName,Marker,Allele,Height]
#' @export

#rm(list=ls())
#setwd("C:\\Users\\oyvbl\\Dropbox\\Forensic\\MixtureProj\\myDev\\seq2lus\\R")
#type = "LUS"
#lookupf = NULL #file to lookup
#toFile = FALSE
#dat=bigdat

convert = function(dat, lookupf=NULL,type="LUS",toFile=FALSE,outf="converted.csv",printMissing=TRUE, missf = "MissingSequences.csv",LUSsymb="_") {
 #LUSsymb = "_" 
 if( require(readxl) ){
 } else {
   stop("The R-package readxl must be installed to continue!")
 }
 pgkPath <- path.package("seq2lus", quiet = FALSE) # Get package path.
 if(is.null(lookupf)) lookupf <- paste(pgkPath ,"LookupTable.xlsx",sep=.Platform$file.sep ) #default path to freq-files
 sheets = excel_sheets(lookupf) #get sheets to read (name of loci)
 sheets = sheets[-1] #ignore legend 
 degn = c("RU","LUS","LUS+","STRseq") #there are different formats
# degn2 = c("LUS+","STRseq") #there are different formats

 locdata = toupper(dat[,2]) #get loci in dataset
 #locdata #BE CAREFUL THAT PENTA IS GIVEN WITH SPACING!
 locdata[toupper(locdata)=="PENTAD"] = "Penta D" 
 locdata[toupper(locdata)=="PENTAE"] = "Penta E" 
 evids = unique(dat[,1]) #get evidences

 avind = grep("ALLELE",toupper(colnames(dat))) #one or more columns? (give same output also)
 if(length(avind)==0) {
  print("NO ALLELES FOUND IN DATASET!")
  return(NULL)
 }
 hvind = grep("HEIG",toupper(colnames(dat)))
 if(length(hvind )==0) hvind = grep("COVE",toupper(colnames(dat)))

 #create an internal list with the data first[[sampleName]][[locus]]
 outL = list() #table to store values SampleName,Marker,Allele,Height,(we don't know max width yet)
 for(evid in evids) outL[[evid]] = list() #init list

 getVal = function(x) { #helpfunction to get alleles/heights
  tmp = unlist(strsplit(as.character(unlist(x)),"/"))
  return(tmp[!is.na(tmp)])
 }

 missingtab = numeric() #list with NEW SEQUENCES
 maxA = 0 #number of alleles
 for(sh in sheets) {
  #sh= sheets[17] #"D16S539" #
  lookdat = read_xlsx(lookupf,sheet=sh)
  cn = colnames(lookdat)

  typ = degn[4]
  if(type%in%degn[1:3]) typ = "LUS+ ALLELE" #WE ONLY USE LUS+ values!

  colSEQ = which("UAS SEQUENCE"==toupper(cn)) #get col index of sequence
  colUSE = which(toupper(typ)==toupper(cn)) #get col index of sequence

  fromDat = unlist(lookdat[,colSEQ])  #from data
  toDat = unlist(lookdat[,colUSE])  #todata
  if( any(is.na(toDat)) ) { #CHECK IF ALL DATA IS THERE
	stop(paste0("Missing data in LookUpTable for locus ",sh))
  }

  if( type=="RU") {
   toDat = sapply(strsplit(toDat,LUSsymb ),function(x) x[1]) #get RU
  } else if( type=="LUS") {
   toDat = sapply(strsplit(toDat,LUSsymb ),function(x) paste0(x[1:2],collapse=LUSsymb) ) #get LUS
  }

  #Create a small lookup table based on sequences in data
  subdat = dat[toupper(locdata)==toupper(sh),,drop=FALSE]
  seqfrom = unique(getVal(subdat[,avind ])) #get unique seqs in data
 
  indmatch = match(seqfrom,fromDat) #get matching indices
  #all(fromDat[indmatch]==seqfrom) #check
  seqto = as.character(toDat[indmatch]) #data to substitute
 #fromDat[indmatch]

  if( any(is.na(seqto)) ){ #if any new sequences found
   newind = !seqfrom%in%fromDat
   newSEQ = seqfrom[newind] #get new sequence
   print( paste0("New sequences found at marker ",sh,". Printing to file..."))
   missingtab = rbind(missingtab,cbind(sh,newSEQ) )

   seqto[newind] = seqfrom[newind] #COPYING SEQUENCE IF MISSING
  }

  #ACCUMULATE AND REPLACE ALLELES  
  for(evid in evids) { #for each evidence
   indevid = which(subdat[,1]==evid)
   if(length(indevid)==0) {
 	print(paste0("No data found for ",evid," at marker ",sh))
     outL[[evid]][[sh]] = list(adata=character(),hdata=character()) #add empty
 	next;
   }
   #if(length(indevid)>1) {
   # print(paste0("Duplicated data found for ",evid," at marker ",sh,". Keeping first!"))
   # indevid = indevid[1] 
   #}
   eviddat = subdat[indevid,,drop=FALSE] #extract data
   av = getVal(eviddat[,avind])
   if(length(av)==0) next #Skip if empty
   hv = as.numeric(getVal(eviddat[,hvind])) #convert back to numbers

   av2 = seqto[match(av,seqfrom )] #new alleles

   if( length(hv)>0 && length(hv)==length(av2) ) {
    agg = aggregate(hv, by=list(av2),sum) #additivity
    av2 = agg[,1] #Updated alleles
    hv2 = agg[,2] #updated coverage
    outL[[evid]][[sh]] = list(adata=av2,hdata=hv2)
   } else {
    outL[[evid]][[sh]] = list(adata=av2)
   }
   if(maxA<length(av2)) maxA = length(av2) #update max
  } #end for each evid
 } #end for each loci

 #STORE TABLE:
 if(length(missingtab)>0) {
  colnames(missingtab) = c("Marker","Sequence")
  write.table(missingtab ,file=missf,sep="\t",row.names=FALSE,quote=FALSE)
 }

 #RETURN RESULTS IN SAME FORMAT AS INPUT:
 outtab = matrix(nrow=0,ncol=maxA) #output table
 outtab2 = matrix(nrow=0,ncol=1) #output table
 colnames(outtab) = paste0("Allele ",1:maxA)
 colnames(outtab2) = paste0("Allele")
 if(length(hvind)>0) {
  outtab = matrix(nrow=0,ncol=2*maxA) #output table
  outtab2 = matrix(nrow=0,ncol=2) #output table
  colnames(outtab) = c(paste0("Allele ",1:maxA),paste0("Height ",1:maxA))
  colnames(outtab2) = c("Allele","Height")
 }
 NameList<- NameList2 <- matrix(nrow=0,ncol=2)
 colnames(NameList) <- colnames(NameList2) <- c("SampleName","Marker")
 for(evid in names(outL)) {#for each evid
  for(loc in names(outL[[evid]])) {#for each locs
   av = outL[[evid]][[loc]]$adata
   nA = length(av)
   nmiss = maxA-nA
   newrow = c(av,rep(NA,nmiss))
   newrow2 = cbind(av) #paste0(av,collapse="/")
   if(length(hvind)>0) {
    hv = outL[[evid]][[loc]]$hdata
    newrow = c(newrow, c(hv,rep(NA,nmiss)))
    newrow2 = cbind(newrow2 , hv) #paste0(hv,collapse="/"))
   }
   outtab = rbind(outtab, newrow )
   NameList = rbind(NameList,c(evid,loc) )
   if(nA>0) {
	outtab2 = rbind(outtab2, newrow2 )   
	NameList2 = rbind(NameList2,cbind(evid,rep(loc,nA)))
   }
  }
 }

 #RETURN SAME FORMAT
 if( length(avind)==1) { 
  outtab = cbind(NameList2,outtab2)
 } else {
  outtab = cbind(NameList,outtab)
 }

 if(toFile) write.table(outtab,file=outf,row.names=FALSE,quote=FALSE,sep="\t") #if print to file

 return(outtab)

} #end function


