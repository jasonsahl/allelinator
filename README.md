# allelinator  

This script takes in a FASTA file and reports the alleles

-The script assumes that the input file will be named in this format  

AQ01_01.fasta (haven't tested case where no underscores are present)  

-The sequences within the file also need to be named to this convention:  
RHIMI-ARA13-xx-DN-020-JB_S86_17441  

This assumes 3 fields, all separated by "_" with the last field being the read count  

-Test files are included with the repository. First you would run the script on an initial file:  
`python allele_driver.py -f TEST_V1.fasta`  

-Then you would run the new file from the same locus, using the allele file produced in the first script  
`python allele_driver.py -f TEST_V2.fasta -a TEST_V1_alleles.tsv`  




