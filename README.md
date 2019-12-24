# unfold_maxlike
Spectrum unfolding with maximum likelihood function 


The program unfolds or inverts spectrum from measured counts.
## Inputs required:
### Essential
- response function 
- measured counts
###  Optional
-theoretical estimate of flux
- guess flux


The file names containing these data are specified in a single `input file`, which is passed as command line argument to the program.

The format of the `input file` is as follows:
-title
-counts,name of the counts data file 
-response,name of the response data file
-theory,name of the theoretically calculated flux data file
-guess, name of the guess flux data file



The format of the response data file is
title
rows,cols,column number of the energy values (1 or last col),string of energy unit
comma or tab or space separated data




The format of the counts data file is
title
rows,cols     (cols = 1 or 2, so that col 1 could be serial number if required)  
comma or tab or space separated data

The format of the theoretical flux data file is
title
rows,cols
comma or tab or space separated data







