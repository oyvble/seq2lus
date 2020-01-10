importData = function(ff) { #function to read data
 require(euroformix)
 return( euroformix::tableReader(ff) )
}