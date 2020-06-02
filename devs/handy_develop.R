#=========================================================================
#handy scripts to be ran during development 
#do not expect any organization or documentation in this file
#=========================================================================



detach("package:msfmethods", unload=TRUE)
system("Rcmd.exe INSTALL --preclean --no-multiarch --with-keep.source C:/Users/worc/Documents/CYER/msfmethods/msfmethods_0.0.0.9000.tar.gz")


devtools::document()

#equivalent to ctrl + b in Rstudio 
devtools::build()



#render Rmd

rmarkdown::render("msf_study.Rmd")

