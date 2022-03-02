#THIS FUNCTION READS SNPinfo from the system

importData = function(ff) { #function to read data
  AT=10
  typedOnly=TRUE
  hasFlanks=FALSE
  
  sheet0="Autosomal STRs" #this is UAS report sheet to extract data from
  suppressMessages(
    dat <- readxl::read_excel(ff, sheet = sheet0, col_names = FALSE, progress=FALSE) #no header
  )
  sn = unlist(dat[3,2]) #get sample name
  headrow = which(dat[,1]=="Coverage Information")+1 #recognize header
  header = dat[headrow,] #extract header
  dat = dat[-seq_len(headrow),,drop=FALSE]
  colnames(dat) = header
  
  #TRAVERSE EACH LOCUS AND EXTRACT RELEVANT SEQUENCES
  locs = unique(dat$Locus)
  addws = "PENTA" #adding whitespace after this string (PentaD/PentaE)
  df = NULL
  for(loc in locs) {
    # loc = locs[27]
    loc2 = toupper(loc)
    if( grepl("AM",loc2) ) next #skip AMEL
    if( grepl(addws,loc2) ) loc2 = gsub(addws,paste0(addws," "),loc2) 
    
    #Obtain AT for marker and Update table regarding AT
    AT0 = AT[loc2]
    if(is.na(AT0)) AT0 = AT
    sub = subset(dat,dat$Locus==loc & as.numeric(dat$Reads)>=AT0) #restrict the coverages here!
    
    if(typedOnly) {
      coluse = grep("typed",tolower(colnames(sub)))
      if(length(coluse)==0) {
        print("Typed column not found. Ignoring!")
      } else if(length(coluse)==1) {
        sub = subset(sub,tolower( sub[[coluse]] )=="yes")
      } else {
        print("Multiple typed column found. Ignoring!")
      }
    }
    
    df_new = cbind(SampleName=sn,Marker=loc2,Allele = sub$`Repeat Sequence`,Height=sub$Reads)
    df = rbind(df, df_new)
  }
  return(df)
}