.onAttach <- function(...) {
  welcome <- paste(""                                              ,
                   "----------------------------------------------",
                   "  'quadrupen' package version 0.4-0           ",
                   ""                                              ,
                   " Still under development... feedback welcome  ",
                   "----------------------------------------------",
                   sep = "\n")
  packageStartupMessage(welcome)
}

