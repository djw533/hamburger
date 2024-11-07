
import os
import sys

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import translate
from Bio.Seq import reverse_complement
from Bio.SeqUtils import GC




def gff2faa(gff_file, fasta_file, strain,output_dir):
    '''convert a gff file with the appended FASTA to
    gff_file = input gff file
    fasta_file = output file
    output: a protein coding FASTA file '''


    out = open(output_dir+"/"+strain+"/tmp_fasta.fa", "w")
    contigs = {}
    list_of_names = {}
    genes_and_contig = {}
    gene_names_altered = 'no'
    filename = gff_file.split('.gff')[0]
    with open(gff_file) as f:
        fasta = False
        counter = 1 # setting hte counter here so that the script has an acutal integer for the gene - and keeping the list so that the gene identifier can be called back later
        for line in f:
            if fasta:
                out.write(line)
                continue
            if line.startswith("##FASTA"):
                fasta = True
                continue
            if line.startswith("#"):
                continue
            toks = line.strip().split("\t")
            if toks[2] != "CDS":
                continue
            #some gffs have had odd symbols - so if the start or stop can't be made into an integer, then continue
            try:
                int(toks[3])
            except:
                continue
            try:
                int(toks[4])
            except:
                continue

            #parse the details in the info column:
            various_details = toks[8]

            #check to see if ID there:
            if "ID=" in toks[8]:
                true_name = toks[8].split("ID=")[1].split(";")[0]
            elif "locus_tag=" in toks[8]:
                true_name = toks[8].split("locus_tag=")[1].split(";")[0]
            else:
                true_name = "No_name_set_{number}".format(number = str(counter))

            #move on to set the hamburger/artificially set gene number (as integer) etc
            name = strain+'_'+str(counter)
            genes_and_contig[name] = toks[0]
            counter += 1
            list_of_names[name] = true_name
            if toks[0] not in contigs:
                contigs[toks[0]] = []

            try:
                contigs[toks[0]].append({"name": name, "start": int(toks[3])-1, "stop": int(toks[4]), "strand": toks[6]})
            except Exception as e:
                raise

    out.close()

    ## read the contigs and save the final fasta file
    out = open(fasta_file, "w")
    with open(output_dir+"/"+strain+"/tmp_fasta.fa") as handle:
        for values in SimpleFastaParser(handle):
            curr_contig = values[0].split()[0]
            #print curr_contig
            if curr_contig not in contigs: # no CDSs in this contig
                #print 'no contig here'
                continue
            for cds in contigs[curr_contig]:
                out.write(">" + cds["name"] + "\n")
                seq = values[1][cds["start"]:cds["stop"]]
                if cds["strand"] == "-":
                    seq = reverse_complement(seq)
                out.write(translate(seq) + "\n")
    out.close()
    os.remove(output_dir+"/"+strain+"/tmp_fasta.fa")
    return gene_names_altered,list_of_names,genes_and_contig



def gff_splitter(gff_file, output_dir):

    gff_data = open(gff_file)
    gff_lines = gff_data.readlines()
    gff_data.close()


    if "##FASTA\n" not in gff_lines:
        print("Couldn't find the fasta sequence in {gff_file} Please use gff3 format otherwise Hamburger will not work".format(gff_file=gff_file))
    fasta_target = gff_lines.index('##FASTA\n')

    annotation_lines = gff_lines[:fasta_target]
    nuc_lines = gff_lines[fasta_target+1:]


    #### get GC content of the whole assembly:

    whole_genome_string = ''.join([x for x in nuc_lines if not x.startswith('>')]).replace("\n","")
    GC_genome=GC(whole_genome_string)


    #write out

    # with open("{output_dir}/annotation.gff".format(output_dir = output_dir),"w") as f:
    #     for line in annotation_lines:
    #         f.write(line)

    with open("{output_dir}/contigs.fna".format(output_dir = output_dir),"w") as f:
        for line in nuc_lines:
            f.write(line)

    return(annotation_lines,GC_genome)


def concat_gff_and_fasta(fasta_file,gff_file,output_gff_file):

    gff_data = open(gff_file)
    gff_lines = gff_data.readlines()
    gff_data.close()

    #check that the annotation isn't already in here...
    if "##FASTA\n" in gff_lines:
        sys.exit("Annotation in the gff file. Can't concatenate a gff with sequence with another fasta file")

    fasta_data = open(fasta_file)
    fasta_lines = fasta_data.readlines()
    fasta_data.close()

    with open(output_gff_file,"w") as output:
        for line in gff_lines:
            output.write(line)
        output.write("##FASTA\n")
        for line in fasta_lines:
            output.write(line)
