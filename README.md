# hamburger

HMmer Based UndeRstandinG of opERons 

A tool to extract and analyse operons or contiguous sets of genes in bacterial genomes, given a concatenated set of protein hidden markov models (HMMs)


=================================================
USAGE:
Extract and plot operons based on hmm profiles
--------------------HaMBURGER--------------------

-------HMmer Based UndeRstandinG of opERons------

optional arguments:
  -h, --help            show this help message and exit
  -i HMM, --hmm HMM     Hmm profile input <required> Set flag
  -g GFF [GFF ...], --gff GFF [GFF ...]
                        Gff file(s) to search <required> Set flag
  -m MIN_GENES, --min_genes MIN_GENES
                        Minimum number of genes in operon, default = 4
  -n GENES_GAP, --genes_gap GENES_GAP
                        Maximum number of genes gap between hits, default = 10
  -u UPSTREAM, --upstream UPSTREAM
                        Number of nucleotides to include upstream of operon,
                        default = 0
  -d DOWNSTREAM, --downstream DOWNSTREAM
                        Number of nucleotides to include downstream of operon,
                        default = 0
  -o OUTPUT, --output OUTPUT
                        Output directory, default = Current date and time

