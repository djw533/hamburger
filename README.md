# hamburger

usage: hamburger.py [-h] -i HMM -g GFF [GFF ...] [-m MIN_GENES] [-n GENES_GAP]
                    [-u UPSTREAM] [-d DOWNSTREAM] [-o OUTPUT]

Extract and plot operons based on hmm profiles
--------------------HaMBURGER--------------------

-------HMmer Based UndeRstandinG of opERons------

              _....----'''----...._
           .-'  o    o    o    o   '-.
          /  o    o     o    o   o    \  	
       __/__o___o_ _ o___ _o _ o_ _ _o_\__
      /                                   \ 	
      \___________________________________/
        \~`-`.__.`-~`._.~`-`~.-~.__.~`-`/
         \                             /
          `-._______________________.-'

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

