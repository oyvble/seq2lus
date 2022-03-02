#THIS FUNCTION READS STR info from the system
importData = function(ff) { #function to read data
 require(readxl)
 sheet0 = "Autosomal STR Coverage" #name of spreadsheet to use
 dat = read_xlsx(ff,sheet=sheet0 ,col_names = FALSE) #no header
 headsize = 9 #header size (include column names)
 dat = as.matrix( dat[-(1:headsize),c(1,3,6,5)]) #remove header and select columns
 colnames(dat) = c("SampleName","Marker","Allele","Height")

 return( dat)
}