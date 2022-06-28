# eDNAFlow Pipeline Species Identification.
# https://github.com/mahsa-mousavi/eDNAFlow#using-lca-script-for-an-otuzotuasvetc-file-that-was-not-generated-by-ednaflow

'''
This script is not limited to eDNAFlow created files. 
It can be used for assigning taxonomy of any OTU, ZOTU, and or ASV files, as long as:

1. User provide a tab-delimited blast result file for their OTU or ASV, 
where blastn is performed using the following format:
-outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp"

'''


