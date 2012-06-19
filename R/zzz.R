.First.lib <-function(lib,pkg)
{
ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
     ver <- as.character(ver)	

options(digits=12)
text = paste("GENEAread", ver, "loaded\n")
packageStartupMessage(text)

##compile bapply and so on if we have a compiler
#require(compiler)
#if (exists("cmpfun")){
#bapply.basic <- cmpfun(bapply.basic)
#bapply <- cmpfun(bapply)
#expand <- cmpfun(expand)
#}

}
