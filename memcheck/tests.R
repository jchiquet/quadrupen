# Code to be run with
#   R -d "valgrind --tool=memcheck --leak-check=full --error-exitcode=1" --vanilla < memcheck/tests.R
library("quadrupen")
devtools::test()
