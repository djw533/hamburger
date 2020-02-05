#Got these reference files from :

#/home/djwilliams/Documents/thesis/chapter_3/1_creation_of_T6SS_models/testing_secret_6_with_macsyfinder/T6SS_hamburge/protein_alignments_and_stats_with_Db11_and_SM39/

#Then made the names equal by removing the gene number from these, using:


#grep '>' tssB.fasta | while read line; do to=$( echo "${line}" | rev | cut -d '_' -f 2- | rev ) ; echo "${to}" ; sed -i "s/${line}/${to}/g" tssB.fasta;  done

### replace tssB.fasta with tssC.fasta as appropriate
