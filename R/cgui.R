#' @title cgui
#' @author Oyvind Bleka
#' @description cgui is a GUI wrapper for converting betweem MPS data formats (sequencies) to LUS
#' @details The function produces a graphical window for letting the user define the settings 
#' @param envirfile A file is given as input (type="file") or data object: possibly type={"list","table"}
#' @export

#envirfile=NULL
cgui = function(envirfile=NULL) {
 require(seq2lus)
 .sep <- .Platform$file.sep 
 pgkPath <- path.package("seq2lus", quiet = FALSE) # Get package path.
 txtformat = c("RU","LUS","LUS+","FWbrack")

 #size of main window
 mwH <- 800
 mwW <- 1000

 #Required for GUI:
 require(gWidgets2tcltk) #requires only gWidgets2.

 #type of gWidgets2-kit
 options(guiToolkit="tcltk")

 #version:
 version =  packageVersion("seq2lus") #follows same version as package number

 #software name:
 softname <- paste0("seq2lus v",version)

 #Spacing between widgets
 spc <- 5
 emptyName = "none (default will be used)" #Text Indicate that nothing is selected
 longspace = "                                                                                 "
 shortspace = "                          "

 ###############
 #HELPFUNCTIONS#
 ###############
 #helpfunction to get environment and file name for different data types
 getEnvirFileNames = function(type) {
   if(type=="DATFOLDS") { #can point to folders
    envirvar = "datfolds"
    fname = DATFOLDS_File
   } else if(type=="DATFILES") {
    envirvar = "datfiles"
    fname = DATFILES_File 
   } else if(type=="SETUP") {
    envirvar = "setup"
    fname = SETUP_File 
   } else {
    envirvar = ""
    fname = ""
   }
   return(c(envirvar,fname ))
 }

 #helpfunction to get folder/IDs vector from list (taken from list envir list only)
  getFolds = function(type) { #get vector of folders from list
   envirvar = getEnvirFileNames(type)[1]
   return( unlist(get(envirvar,envir=mmTK)) ) #return vector of folds
  }

  #helpfunction to set folders/IDs to environment and file (both envir list and file is updated with new info)
  addFold = function(foldadd,type="DATFOLDS") {
   tmp =  getEnvirFileNames(type)
   envirvar = tmp[1]
   fname = tmp[2]
 
   #Update environment variable:
   X = get(envirvar,envir=mmTK) #get list
   X[[length(X)+1]] = foldadd #add file to envir list
   assign(envirvar,X,envir=mmTK) #Store setup values

   #Store to file (update file):
   X = unlist(X)
   write(X,file=fname )    #save to file in installation folder 
   return(X) #return vector of folderes
  }

  saveFolds = function(folds,type="DATFOLDS") {
   tmp =  getEnvirFileNames(type)
   envirvar = tmp[1]
   fname = tmp[2]

   #Store to file (update file):
   write(folds,file=fname )    #save to file in installation folder 

   #Update environment variable:
   assign(envirvar,as.list(folds),envir=mmTK) #Store setup values
  }

  #Helpfunction to store/load values in optList to file
  saveSetup = function(opt) { 
   write(unlist(opt),file=SETUP_File)
  }
  openSetup = function() { 
   #vars = names(optL)
   vars = c("workdir","lookupfile","importDataFile","outputfold","outputname","missingfile","format")
   vart = rep("s",length(vars)) #variable types (s=string,b=boolean,d=double,i=integer)
   dat = scan(file=SETUP_File,what=character(),quiet=TRUE,sep="\n")
   opt = list() #init list
   for(i in 1:length(vars)) { #
      x = dat[i] #string is standard
      if(vart[i]=="b") {
        x = as.logical(x)
      } else if(vart[i]=="d") {
	   x = as.numeric(x)
      } else if(vart[i]=="i") {
	   x = as.integer(x)
      }
      opt[[vars[i]]] = x  #insert correct type of variable for each list element
   }
   return(opt)
  }

 errorMessage = function(msg) {  #Helpfunction to throw error message to user + stop running 
  gWidgets2::gmessage(msg,title="Error",icon="error")
  stop(msg)
 }

 NullIfEmpty = function(x) {
   if(length(x)==0) {
    return(NULL)
   } else {
    return(x)
   }
 }

############################################################
#FUNCTION TO RUN ANALYSIS (running with specified settings)#
############################################################
  runConvertion = function(h,...) { #CALL CONVERTION
    #Step 1: Read search setup from envir vars
    #Step 2: Run str2LUS convert with setup

  	datfoldList = NullIfEmpty(get("datfolds",envir=mmTK))
  	datfileList = NullIfEmpty(get("datfiles",envir=mmTK))
  	opt = get("setup",envir=mmTK)  #receive settings from envir (preassigned or from file)
  
    #prechecks:
  	if(opt$lookupfile==emptyName) {
      print("No lookup file was selected. Using default from package...")
  	  lookupfile <- NULL
  	} else {
   	  lookupfile = opt$lookupfile #use selected fone
  	}
  
  	#READ DATA AND MERGE DATA FROM SELECTED FOLDERS/FILES AND PUT INTO CONVERTER	
    outf = paste0(opt$outputfold,.sep,opt$outputname,".csv")
    
    if(opt$importDataFile!=emptyName) { #import data function was not selected
      source(opt$importDataFile) #substute importData with the one found in a selected R-file 
    } else if( require(euroformix) ) {
      importData = euroformix::tableReader
    } else {
      importData = function(X) read.table(X,header=TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE) #default data reader
    }
  
  	#Gather path to all files to read (add with those inside selected folders)
  	datFiles = datfileList
  	for(fold in datfoldList) datFiles = c(datFiles, list.files(fold,full.names=TRUE,recursive=FALSE)) #add files in selected folders
  	datFiles = unique(unlist(datFiles)) #get unique files to import
  
  	bigdat = numeric() #assume same number of columns
    for(datfile in datFiles) {
  	tryCatch( {
    	dat = importData(datfile)
  	  bigdat = rbind(bigdat , dat)
  	 }, error = function(e) print(e)) #skip files giving errors when importing
    }
  
  	if(nrow(bigdat)==0) { #handle the situation of no data
    	gWidgets2::gmessage("No data to convert was found!")		
    	return()
  	}
  	
  	#"Missing-allele" file (can be named in GUI):
  	missf =  paste0(opt$outputfold,.sep,opt$missingfile,".csv")
    tab <- seq2lus::convert(bigdat,lookupf=lookupfile,type=opt$format,toFile=TRUE,outf=outf,printMissing=TRUE, missf = missf)#,LUSsymb=opt$LUSsymb)
    #write.table(tab,file=outf,row.names=FALSE,sep="\t",quote=FALSE) # SHOULD USER SELECT SEPARATOR?
  
    gWidgets2::gmessage("Conversion successfully completed!")
  } #end run function

 #####################
 #create environment #
 #####################

 #STORING INFO IN BOTH Files and Environment (easy user access)
 #File variable - envir variable:
 #DATFOLDS_File (configDATFOLDS) - DATFOLDS
 #DATFILES_File (configDATFILES) - DATFILES
 #SETUP_File (configSETUP) - SETUP


 DATFOLDS_File <- paste(pgkPath,"configDATFOLDS",sep=.sep) #Setting file for data folders
 DATFILES_File <- paste(pgkPath,"configDATFILES",sep=.sep) #Setting file for data files
 SETUP_File <- paste(pgkPath,"configSETUP",sep=.sep) #Create a file with all settings (not lists)

 if(is.null(envirfile)) {
  mmTK = new.env( parent = emptyenv() ) #create new environment object

  DATFOLDS_List <- DATFILES_List <- list() #This is default (none selected)
  if(file.exists(DATFOLDS_File)) DATFOLDS_List <- scan(file=DATFOLDS_File,what=character(),quiet=TRUE,sep="\n")
  if(file.exists(DATFILES_File)) DATFILES_List <- scan(file=DATFILES_File,what=character(),quiet=TRUE,sep="\n")

  #Default set (empty) of folds:
  assign("datfolds",DATFOLDS_List,envir=mmTK)
  assign("datfiles",DATFILES_List,envir=mmTK)

  #Default Settings set if setupfile not found:
  if(file.exists(SETUP_File)) {
   opt <- openSetup() #get setup list
   assign("setup",opt,envir=mmTK) #Store setup values 
  } else { #DEFAULT SETUP VALUES:
   opt = list() #list of options
   opt$workdir = getwd() #default is work directory
   opt$lookupfile = emptyName #name of freq file (should be full path?). Use file selector 
   opt$importDataFile = emptyName 
   opt$outputfold = getwd() #use work directory as default
   opt$outputname = "outputdata" 
   opt$missingfile = "MissingLookupAlleles"
   opt$format = "LUS" #this is default format
   #opt$LUSsymb = "_" #this is to-format (not used)
   assign("setup",opt,envir=mmTK) #Store setup values 
  } 
 }else {
  load(envirfile) #loading environment
 }
 #optL = get("setup",envir=mmTK)  #receive settings from envir (preassigned or from file)

###################################################################
###########################GUI#####################################
###################################################################

 #Menu bar file-lists:
 f_setwd = function(h,...) {
  dirsel = gWidgets2::gfile(text="Select folder",type="selectdir")
  if(!is.na(dirsel)) {
   setwd(dirsel)
   opt = get("setup",envir=mmTK) #get
   opt$workdir = dirsel
   assign("setup",opt,envir=mmTK) #set
   saveSetup(opt) #save to setup file
  }
 }
 f_openproj = function(h,...) {
  projfile = gWidgets2::gfile(text="Open settings",type="open")
  if(!is.na(projfile)) {
   gWidgets2::dispose(mainwin)
   seq2lus::cgui(projfile) #send environment into program
  }
 }
 f_saveproj = function(h,...) {
  projfile = gWidgets2::gfile(text="Save settings",type="save")
  if(!is.na(projfile) && length(projfile)>0) {
   save(mmTK,file=projfile) #save environment
   print(paste("Settings saved in ",projfile,sep=""))
  }
 }
 f_quitproj = function(h,...) {
  ubool <- gWidgets2::gconfirm("Do you want to save project?",title="Quit Program",icon="info")
  if(ubool) {
   f_saveproj()
  } else { 
   print("Program terminated without saving")
  }
  gWidgets2::dispose(mainwin) #remove window!
 }


 #helpfunction for adding Folder/File when clicking button
  f_addFolder = function(h,...) {
    if(h$action=="DATFOLDS") {
     dirfile = gWidgets2::gfile(text="Select folder",type="selectdir")
     if(is.na(dirfile)) return() 
     fv = addFold(foldadd=dirfile,h$action) #Add folder to environment and file, h=list(action="EVID")
     tab1a[1,2][] = fv #update combolist
     gWidgets2::enabled(tab1a[2,2]) = TRUE
    }

    if(h$action=="DATFILES") {
     dirfile = gWidgets2::gfile(text="Select file",type="open")
     if(is.na(dirfile)) return() 
     fv = addFold(foldadd=dirfile,h$action) #Add folder to environment and file, h=list(action="EVID")
     tab1b[1,2][] = fv #update combolist
     gWidgets2::enabled(tab1b[2,2]) = TRUE
    }
  }
 
  #helpfunction for deleting marked folder when clicking button
  f_delFolder = function(h,...) {
      if(h$action=="DATFOLDS") {
        sel = gWidgets2::svalue(tab1a[1,2])
        vals = tab1a[1,2][]
        folds = setdiff(vals,sel)
        tab1a[1,2][] = folds  #folds #update combolist
        if(length(folds)==0) {
          tab1a[1,2][] = longspace #folds #update combolist
          gWidgets2::svalue(tab1a[1,2]) = longspace 
          gWidgets2::enabled(tab1a[2,2]) = FALSE
        } else {
          gWidgets2::svalue(tab1a[1,2]) = folds[1]
        }
      }
      if(h$action=="DATFILES") {
        sel = gWidgets2::svalue(tab1b[1,2])
        vals = tab1b[1,2][]
        folds = setdiff(vals,sel)
        tab1b[1,2][] = folds #update combolist
        if(length(folds)==0) {
          tab1b[1,2][] = longspace #folds #update combolist
          gWidgets2::svalue(tab1b[1,2]) = longspace 
          gWidgets2::enabled(tab1b[2,2]) = FALSE
        } else {
          gWidgets2::svalue(tab1b[1,2]) = folds[1]
        }
      }
      saveFolds(folds,h$action) #get  h=list(action="EVID")
  } #End delFolder function


##################################################################################################
########### Program starts #######################################################################
##################################################################################################

 ###############
 #start layout:#
 ###############
 mblst = list( #NOTICE THE NEW CODE IN gWidgets2
  File=list(  
    gWidgets2::gaction('Set directory',handler=f_setwd),
    gWidgets2::gaction('Open project', handler=f_openproj),
    gWidgets2::gaction('Save project', handler=f_saveproj),
    gWidgets2::gaction('Quit', handler=f_quitproj,icon="close")
  )
 )

 #change working directory to the one stored in mmTK-environment
 wd=get("setup",envir=mmTK)$workdir #assign working directory to mmTK-environment
 if(!is.null(wd)) {
   tryCatch( { setwd(wd) }, error=function(e) print("Folder not found. Using existing") )
 }
 
 #Main window:
 mainwin <- gWidgets2::gwindow(softname, visible=FALSE, width=mwW,height=mwH)
 gWidgets2::addHandlerUnrealize( mainwin, handler = function(h,...) {
	bool = gWidgets2::gconfirm("Are you sure you want to quit?") 
     if(bool) {
       gWidgets2::dispose(mainwin) #remove window!
     } else {
	  return(TRUE)
     }
 }  ) #call quit function 
 gWidgets2::gmenu(mblst,container=mainwin)
 nb = gWidgets2::gnotebook(container=mainwin)
 tabconvert = gWidgets2::ggroup(expand=TRUE,spacing=spc,container=nb,label="Convert") #tab1: (select project and file storage)
 gWidgets2::svalue(nb) <- 1 #initial start in first tab


#####################################################
###############Tab 1: Convert :######################
#####################################################

  tab1 <- gWidgets2::glayout(spacing=spc,container=tabconvert ) 


  tab1a = gWidgets2::glayout(spacing=spc,container=(tab1[1,1] <-gWidgets2::gframe("Selected folders with data",container=tab1))) 
  tab1a[1,1] <- gWidgets2::glabel("Selected folders:",container=tab1a)
  tab1a[2,1] <- gWidgets2::gbutton("Add a folder",container=tab1a,handler=f_addFolder,action="DATFOLDS")
  tab1a[2,2] <- gWidgets2::gbutton("Remove marked folder",container=tab1a,handler=f_delFolder,action="DATFOLDS")
  folds = getFolds("DATFOLDS") 
  if(is.null(folds) || length(folds)==0)  folds = longspace #numeric()
  tab1a[1,2] <- gWidgets2::gcombobox(items=folds,container=tab1a)
  gWidgets2::size(tab1a[1,2]) = nchar(longspace)
  if(folds[1] == longspace) gWidgets2::enabled(tab1a[2,2]) = FALSE

  tab1b = gWidgets2::glayout(spacing=spc,container=(tab1[2,1] <-gWidgets2::gframe("Selected data files",container=tab1))) 
  tab1b[1,1] <- gWidgets2::glabel("Selected files:",container=tab1b)
  tab1b[2,1] <- gWidgets2::gbutton("Add a file",container=tab1b,handler=f_addFolder,action="DATFILES")
  tab1b[2,2] <- gWidgets2::gbutton("Remove marked file",container=tab1b,handler=f_delFolder,action="DATFILES")
  folds = getFolds("DATFILES") 
  if(is.null(folds) || length(folds)==0)  folds = longspace #numeric()
  tab1b[1,2] <- gWidgets2::gcombobox(items=folds,container=tab1b)
  gWidgets2::size(tab1b[1,2]) = nchar(longspace)
  if(folds[1] == longspace) gWidgets2::enabled(tab1b[2,2]) = FALSE

  tab1c = gWidgets2::glayout(spacing=spc,container=(tab1[3,1] <-gWidgets2::gframe("Import Data file (customizing data format)",container=tab1))) 
  tab1c[1,1] <- gWidgets2::gbutton("Selected importData file:",container=tab1c,handler = 
	function(h,...) { 
      fsel = gWidgets2::gfile(text="Select importData file",type="open")
      if(!is.na(fsel)) {
       opt = get("setup",envir=mmTK) #get
       opt$importDataFile = fsel
       assign("setup",opt,envir=mmTK) #set to envir
       saveSetup(opt) #Save to file
       gWidgets2::svalue(tab1c[1,2]) = fsel
      }
  })
  tab1c[1,2] <- gWidgets2::glabel(get("setup",envir=mmTK)$importDataFile,container=tab1c)
  tab1c[2,1] <- gWidgets2::gbutton("Set back to default",container=tab1c,handler = 
	function(h,...) { 
       opt = get("setup",envir=mmTK) #get
       opt$importDataFile = emptyName
       assign("setup",opt,envir=mmTK) #set to envir
       saveSetup(opt) #Save to file
       gWidgets2::svalue(tab1c[1,2]) = emptyName
     })


  tab1d = gWidgets2::glayout(spacing=spc,container=(tab1[4,1] <-gWidgets2::gframe("Lookup table file selection",container=tab1))) 
  tab1d[1,1] <- gWidgets2::gbutton("Selected LookupTable file:",container=tab1d,handler = 
	function(h,...) { 
      fsel = gWidgets2::gfile(text="Select file (xlsx)",type="open")
      if(!is.na(fsel)) {
       opt = get("setup",envir=mmTK) #get
       opt$lookupfile = fsel
       assign("setup",opt,envir=mmTK) #set to envir
       saveSetup(opt) #Save to file
       gWidgets2::svalue(tab1d[1,2]) = fsel
      }
  })
  tab1d[1,2] <- gWidgets2::glabel(get("setup",envir=mmTK)$lookupfile ,container=tab1d)
  tab1d[2,1] <- gWidgets2::gbutton("Set back to default",container=tab1d,handler = 
	function(h,...) { 
       opt = get("setup",envir=mmTK) #get
       opt$lookupfile = emptyName 
       assign("setup",opt,envir=mmTK) #set to envir
       saveSetup(opt) #Save to file
       gWidgets2::svalue(tab1d[1,2]) = emptyName
     })

  tab1e = gWidgets2::glayout(spacing=spc,container=(tab1[5,1] <-gWidgets2::gframe("Output",container=tab1))) 
  tab1e[1,1] <- gWidgets2::gbutton("Selected output folder:",container=tab1e,handler = 
	function(h,...) { 
      fsel = gWidgets2::gfile(text="Select folder",type="selectdir")
      if(!is.na(fsel)) {
       opt = get("setup",envir=mmTK) #get
       opt$outputfold = fsel
       assign("setup",opt,envir=mmTK) #set to envir
       saveSetup(opt) #Save to file
       gWidgets2::svalue(tab1e[1,2]) = fsel
      }
  })
  tab1e[1,2] <- gWidgets2::glabel(get("setup",envir=mmTK)$outputfold,container=tab1e)

  tab1e[2,1] <- gWidgets2::gbutton("Selected output filename:",container=tab1e,handler = 
	function(h,...) { 
      fsel = gWidgets2::ginput("Set name of converted file",get("setup",envir=mmTK)$outputname )
      if(!is.na(fsel)) {
       opt = get("setup",envir=mmTK) #get
       opt$outputname = fsel
       assign("setup",opt,envir=mmTK) #set to envir
       saveSetup(opt) #Save to file
       gWidgets2::svalue(tab1e[2,2]) = fsel
      }
  })
  tab1e[2,2] <- gWidgets2::glabel(get("setup",envir=mmTK)$outputname,container=tab1e)

  tab1e[3,1] <- gWidgets2::gbutton("Selected missing alleles filename:",container=tab1e,handler = 
    function(h,...) { 
        fsel = gWidgets2::ginput("Set name of file with missing alleles \n(not found in lookup table)",get("setup",envir=mmTK)$missingfile )
        if(!is.na(fsel)) {
          opt = get("setup",envir=mmTK) #get
          opt$missingfile = fsel
          assign("setup",opt,envir=mmTK) #set to envir
          saveSetup(opt) #Save to file
          gWidgets2::svalue(tab1e[3,2]) = fsel
        }
    })
  tab1e[3,2] <- gWidgets2::glabel(get("setup",envir=mmTK)$missingfile,container=tab1e)

  tab1f = gWidgets2::glayout(spacing=spc,container=(tab1[6,1] <-gWidgets2::gframe("Convert",container=tab1))) 
  formatsel = which(get("setup",envir=mmTK)$format==txtformat)
  if(length(formatsel)==0) formatsel = 2 #this is default
  tab1f[1,1] <- gWidgets2::gradio(items=txtformat,container=tab1f ,horizontal = TRUE, selected=formatsel, handler = 
    function(h,...) {
      opt = get("setup",envir=mmTK) #get
      opt$format = gWidgets2::svalue(tab1f[1,1])
      assign("setup",opt,envir=mmTK) #set to envir
      saveSetup(opt) #Save to file
    })
  tab1f[2,1] <- gWidgets2::gbutton("PERFORM CONVERSION",container=tab1f,handler = runConvertion)
  gWidgets2::visible(mainwin) = TRUE

} #end outer function
