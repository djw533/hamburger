#!/usr/bin/env python3

import traceback, os, math, sys

import os.path


from datetime import datetime


### import the bio requirements
import Bio
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import translate
from Bio.Seq import reverse_complement
from Bio import SeqIO
from Bio.SeqUtils import GC
import time
from shutil import which


### multiprocessing
import multiprocessing as mp
from multiprocessing import Process
import multiprocessing
import tqdm

#print("Number of processors: ", mp.cpu_count())

#from multiprocessing import Process

#import threading, time

###### modify the path/name for hmmsearch below:

hmmsearch='hmmsearch'

# ### need the directory where hamburger is installed:
hamburger_base_directory = os.path.abspath(__file__).split("/scripts/hamburger.py")[0]


### and for the supplied models for standard pipeline & parameters

T6SS_core = "{hamburger_base_directory}/models/T6SS/T6SS_core.hmm".format(hamburger_base_directory=hamburger_base_directory)
T6SS_accessory = "{hamburger_base_directory}/models/T6SS/T6SS_accessory.hmm".format(hamburger_base_directory=hamburger_base_directory)


###### parse arguments

def parseArgs():
    """Parse command line arguments"""

    import argparse

    try:
        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter ,
            description="""Extract and plot gene_clusters based on hmm profiles
--------------------HaMBURGER--------------------\n
\n


-------HMmer Based UndeRstandinG of gene clustERs------\n
\n
              _....----'''----...._
           .-'  o    o    o    o   '-.
          /  o    o     o    o   o    \  \t
       __/__o___o_ _ o___ _o _ o_ _ _o_\__
      /                                   \ \t
      \___________________________________/
        \~`-`.__.`-~`._.~`-`~.-~.__.~`-`/
         \                             /
          `-._______________________.-'
""")
        parser.add_argument('-i',
                '--mandatory',
                action='store',
            #    required=True,
                help='Mandatory hmm profile input <required if not using -t flag> ')
        parser.add_argument('-a',
                '--accessory',
                action='store',
                required=False,
                help='Accessory hmm profile input')
        parser.add_argument('-g',
			    '--gff',
			    action='store',
                nargs='+',
                #help='<Required> Set flag',
                required=True,
			    help='Gff file(s) to search <required>')
        parser.add_argument('-m',
			    '--min_genes',
			    action='store',
                default=4,
			    help='Minimum number of genes in gene cluster, default = 4')
        parser.add_argument('-l',
			    '--genes_gap',
			    action='store',
                default=10,
			    help='Maximum number of genes gap between hits, default = 10')
        parser.add_argument('-u',
			    '--upstream',
                type = int,
			    action='store',
                default=0,
			    help='Number of nucleotides to include upstream/"right" of gene cluster, default = 0')
        parser.add_argument('-d',
			    '--downstream',
                type = int,
			    action='store',
                default=0,
			    help='Number of nucleotides to include downstream/"left" of gene cluster, default = 0')
        parser.add_argument('-c',
			    '--cutoff',
                type = float,
			    action='store',
                default=20,
			    help='Cutoff HMMER score for each hit, default = 20')
        parser.add_argument('-t',
			    '--t6ss',
			    action='store_true',
			    help='Automatic searching for T6SSs, uses min_genes = 8, genes_gap = 12, mandatory hmm profile of all 13 tss genes')
        parser.add_argument('-n',
			    '--num_threads',
			    action='store',
                default = 1,
			    help='Number of threads to use, default = 1')
        parser.add_argument('-k',
			    '--keep_files',
			    action='store_true',
			    help='Keep all intermediate files produced, default = False')
#        parser.add_argument('-p',
#			    '--pfam',
#			    action='store',
#                default = False
#			    help='Text file of pfam domains to use for hmms - one per line')
        parser.add_argument('-o',
			    '--output',
			    action='store',
                default="Hamburger_output",
                #default="HaMBURGER_output_{time}".format(time='_'.join('_'.join(str(datetime.now()).split('.')[0].split(':')).split())),
			    help='Output directory, default =Hamburger output. Will not write over a previously existing output folder!')
        parser.add_argument('-q',
                '--itol',
                action='store_true',
                help='Create itol output for number of T6SSs and subtypes per strain')

    except:
        print("An exception occurred with argument parsing. Check your provided options.")
        traceback.print_exc()

    return parser.parse_args()


def __is_tool__(name):
    """Check whether `name` is on PATH and marked as executable."""

    # from whichcraft import which

    return which(name) is not None



def __read_hmms__(input_hmm_file):
    ### open the hmm file
    hmm_data = open(input_hmm_file)
    hmm_lines = hmm_data.readlines()
    hmm_data.close()

    hmm_names = []

    ### get all the names of the entries:
    for l in hmm_lines:
        if l.startswith("NAME"):
            name = l.split()[1].strip() # take the name and add it to the list
            hmm_names.append(name)

    ### return the list of the names:
    return hmm_names

def __gff2faa__(gff_file, fasta_file, strain,output_dir):
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
            true_name = toks[8].split("ID=")[1].split(";")[0]
            name = strain+'_'+str(counter)
            genes_and_contig[name] = toks[0]
            counter += 1
            list_of_names[name] = true_name
            if toks[0] not in contigs:
                contigs[toks[0]] = []
            contigs[toks[0]].append({"name": name, "start": int(toks[3])-1, "stop": int(toks[4]), "strand": toks[6]})
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

def __run_hmmer__(hmmsearch_command, hmm_queries, input, hmmer_cutoff, output_dir, name):
    # run hmmer, apply cutoff to only get the results that are wanted, and output the hmmer results
    os.system("{hmmsearch} -o {output_dir}/{name}_log --tblout {output_dir}/{name}_hmmer_scores.tbl {hmm_queries} {input}".format(hmmsearch = hmmsearch_command, hmm_queries = hmm_queries, input = input, output_dir = output_dir, name = name))

    ## Sometimes this doesn't seem to be working - so if the file doesn't exist - then wait and repeat until it does:
    while not os.path.exists("{output_dir}/{name}_hmmer_scores.tbl".format(name = name, output_dir = output_dir)):
        os.system("{hmmsearch} -o {output_dir}/{name}_log --tblout {output_dir}/{name}_hmmer_scores.tbl {hmm_queries} {input}".format(hmmsearch = hmmsearch_command, hmm_queries = hmm_queries, input = input, output_dir = output_dir, name = name))




    ### create list for the hmmer output
    hmmer_output = []

    ##create list for the gene names for clustering later:
    gene_names = []

    ### read in the hmmer_output
    hmmer_output_data = open("{output_dir}/{name}_hmmer_scores.tbl".format(name = name, output_dir = output_dir))
    hmmer_lines = hmmer_output_data.readlines()
    hmmer_output_data.close()

    ### filter the hmmer_output

    for line in hmmer_lines:
        if line.startswith("#"):
            pass
        else:
            toks = line.strip().split()
            hmmer_score = float(toks[5])
            gene_number = int(toks[0].split('_')[-1])

            ####only take the hit if the score is above the threshold:
            if hmmer_score >= hmmer_cutoff:
                hmmer_output.append(line)
                gene_names.append(gene_number)

    #now return the filtered hmmer scores:

    return(hmmer_output,gene_names)

def __clustering_func__(numbers_list,max_gap,groupsize_min):
    ### Setting maximum gap between genes in the gene cluster

    maxgap = int(max_gap)

    ### clustering CDSs based on position and maxgap value

    numbers_list.sort()
    #print numbers_list
    groups = [[numbers_list[0]]]
    for x in numbers_list[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])

    ### filter groups based on the min_genes value

    filtered_groups = []


    for group in groups:
        if len(group) >= int(groupsize_min):
            filtered_groups.append(group)
    return filtered_groups

def __gff_splitter__(gff_file, output_dir):

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

def __multifasta_to_singlefasta__(input_fasta, input_dir, output_dir):

    os.makedirs(output_dir)

    data = open("{input_dir}/{input_fasta}".format(input_dir = input_dir, input_fasta = input_fasta))

    lines = data.readlines()

    lines_string = ''.join(lines)

    entries_list = (''.join(lines)).split('>')

    for entry in entries_list:
         name = entry.split('\n')[0]
         ## going to try to just take the first part to the first break as the output name for the singlefasta.. (if there is more than one):
         if len(name.split()) > 1:
             name = name.split()[0]
         with open("{output_dir}/{name}.fna".format(output_dir = output_dir, name = name), "w") as output:
              output.write('>'+entry)

def __same_contigs_check__(gene_numbers, strain, genes_and_contig):

    contigs = []

    for num in gene_numbers:
        gene = "{strain}_{num}".format(strain = strain, num = str(num)) # covert the number back into the format in the annotation

        ### get the contig:
        contig = genes_and_contig[gene]

        contigs.append(contig)

    ### check if all the contigs are the same:
    if all(x == contigs[0] for x in contigs) == True:
        return(True, contig)

    else:
        return(False, "No contig")

def __gene_names_in_cluster__(gene_num_in_cluster, hmmer_output, strain):

    ###first get a dictionary with the gene number and the hmmer gene that was hit:

    ### also think - what if two hmmer queries match the same protein at some level?

    names_to_genes = {}

    for line in hmmer_output:
        gene_num = int(line.split()[0].split(strain+"_")[1]) # take the number after the strain, e.g. {strain}_0234 or {strain}_05422
        hmmer_query = line.split()[2]
        #put in dict:
        names_to_genes[gene_num] = hmmer_query


    ## cycle through the gene_numbers in the gene cluster and count number of each hmmer query type

    gene_names = {}


    for gene in gene_num_in_cluster:
        hmmer_query_name = names_to_genes[gene]

        ## add to the names_to_genes, or create key if not in the dict:
        if hmmer_query_name not in gene_names:
            gene_names[hmmer_query_name] = []

        gene_names[hmmer_query_name].append(gene)


    ### now report the number of unique gene hits found in the cluster:
    return(gene_names)

def extract_a2b(start_gene, stop_gene, annotation, strain, cluster_num, upstream, downstream, contig, output_dir, orientation):
    """ Take sequence and annotation from gff3 files between genomic point A and B, given both base indices
    Orientation should be either 'same', 'forward', or 'reverse' """



    #check the contig size is long enough:
    record = SeqIO.read("{output_dir}/{strain}/singlefastas/{contig}.fna".format(output_dir = output_dir, strain = strain, contig = contig), "fasta")
    contig_length = len(record)



    ### creating list for the extracted annotation for later use:
    cluster_annotation = []

    start_got = False
    stop_got = False

    ### check in the start and the stop are the same...

    if orientation == "same":
        for gff_line in annotation:
            if gff_line.startswith('#'):
                continue # ignore the hashed lines

            id = gff_line.split("\t")[8].split(";")[0].split("=")[1]

            if id == start_gene:
                feature_start = gff_line.split('\t')[3]
                feature_stop = gff_line.split('\t')[4]
                start_index = min(feature_start, feature_stop)
                stop_index = max(feature_start, feature_stop)
                start_got = True
                stop_got == True


    else:
        ## find start and stop
        for gff_line in annotation:
            if gff_line.startswith('#'):
                continue

            id = gff_line.split("\t")[8].split(";")[0].split("=")[1]

            if id == start_gene:
                feature_start = gff_line.split('\t')[3]
                feature_stop = gff_line.split('\t')[4]
                if orientation == "forward":
                    start_index = min(feature_start, feature_stop)
                elif orientation == "reverse":
                    start_index = max(feature_start, feature_stop)
                start_got = True

            elif id == stop_gene:
                feature_start = gff_line.split('\t')[3]
                feature_stop = gff_line.split('\t')[4]
                stop_index = min(feature_start, feature_stop)
                if orientation == "forward":
                    stop_index = max(feature_start, feature_stop)
                elif orientation == "reverse":
                    stop_index = min(feature_start, feature_stop)
                stop_got = True

            elif stop_got == True and start_got == True:
                continue


    #print 'Checking to see if hit is on the positive or negative strand: -----'
    ### N.B this should always be positive - this is a bit that has been taken over from a previous script that equated blast hits (which had been reversed by blast)
    ### to gff subsections from the original annotated genome
    result = int(stop_index) - int(start_index)
    if result > 0:
        #print 'Hit is on the positive strand, proceeding : ------'
        start = (int(start_index) - 1) - int(upstream)
        stop = int(stop_index) + int(downstream)
        gc_start = (int(start_index) - 1)
        gc_stop = int(stop_index)

        #check if the contig is long enough:
        if start < 0 :
            start = 0
        if stop > contig_length:
            stop = contig_length



    elif result < 0:
        #print 'Hit is on the negative strand, redefining the start and stop indices'
        start = int(stop_index) - int(downstream)
        stop = (int(start_index) -1 ) + int(upstream)
        gc_start = int(stop_index)
        gc_stop = (int(start_index) -1 )


        #again check if the contig is long enough: (and correct if not):
        if start < 0 :
            start = 0
        if stop > contig_length:
            stop = contig_length



    elif result == 0:
        print('Start and Stop are the same, break')


    ## -- Next line of gff file written detailing the sequence region, how long the region is, and what the name of the chromosome is in the associated fasta sequence
    with open("{output_dir}/{strain}/{strain}_cluster_{cluster_num}.gff".format(output_dir = output_dir, strain = strain, cluster_num = cluster_num), "w") as output:
        output.write("##sequence-region {contig} 1 {end_pos}\n".format(contig = contig, end_pos = str(stop-(start+1))))

    for gff_line in annotation:
        if gff_line.startswith('#'):
            continue
        else:
            feature_start = gff_line.split('\t')[3]
            feature_stop = gff_line.split('\t')[4]

        if int(feature_start) >= start and int(feature_stop) <= stop and contig in gff_line:
            new_list = gff_line.split("\t")[0:3] + gff_line.split("\t")[5:8]
            new_list.append(" ".join(gff_line.strip().split("\t")[8:]))
            new_list.insert(3, str(int(feature_start) - start))
            new_list.insert(4, str(int(feature_stop) - start))
            linegff3 = "\t".join(new_list)
            with open("{output_dir}/{strain}/{strain}_cluster_{cluster_num}.gff".format(output_dir = output_dir, strain = strain, cluster_num = cluster_num), "a") as output:
                output.write(linegff3)
                cluster_annotation.append(linegff3) ### append only the cluster annotation for later gggenes input csv creation
                output.write("\n") ####### have removed the newline character "\n"

        else:
            pass

    ## -- Append ##fasta header to the annotation to indicate that the fasta sequence will be below

    with open("{output_dir}/{strain}/{strain}_cluster_{cluster_num}.gff".format(output_dir = output_dir, strain = strain, cluster_num = cluster_num), "a") as output:
        output.write('##FASTA\n')


    ## -- Append fasta sequence to the end of the gff file

    record = SeqIO.read("{output_dir}/{strain}/singlefastas/{contig}.fna".format(output_dir = output_dir, strain = strain, contig = contig), "fasta")
    with open("{output_dir}/{strain}/{strain}_cluster_{cluster_num}.gff".format(output_dir = output_dir, strain = strain, cluster_num = cluster_num), "a") as out:
        SeqIO.write(record[start:stop], out, "fasta")

    with open("{output_dir}/{strain}/{strain}_cluster_{cluster_num}.fna".format(output_dir = output_dir, strain = strain, cluster_num = cluster_num), "w") as out:
        SeqIO.write(record[start:stop], out, "fasta")




    ### get the GC content of the operon from the hit to the rest of the hit as well
    with open("{output_dir}/{strain}/temp_gc_string.fna".format(output_dir = output_dir, strain = strain, cluster_num = cluster_num), "w") as out:
        SeqIO.write(record[gc_start:gc_stop], out, "fasta")

    cluster_nuc_lines = open("{output_dir}/{strain}/temp_gc_string.fna".format(output_dir = output_dir, strain = strain, cluster_num = cluster_num)).readlines()
    cluster_sequence_string = ''.join([x for x in cluster_nuc_lines if not x.startswith('>')]).replace("\n","")
    GC_cluster=GC(cluster_sequence_string)

    return(GC_cluster, gc_start, gc_stop, start, stop, cluster_annotation)

def Merge(dict1, dict2):
    return(dict2.update(dict1))

def __extract_gggenes_info__(cluster_annotation, named_genes,list_of_genes,start,stop,strain,output_dir,cluster_num):

    ###match the gene names back to what they are in the annotation

    old_named_genes = {}

    for model, new_numbers in named_genes.items():
        for num in new_numbers:
            old_named_genes[list_of_genes[strain+"_"+str(num)]] = model # create new dict with each gene numebr rin the original annotation having the model assigned to it

    #### now go through the list of the cluster annotation and put all into the gggenes csv
    gg_number = 0
    gggenes_output = []

    for line in cluster_annotation:
        if line.split("\t")[2] == "CDS":


            ## update counter:
            gg_number += 1

            gg_start = str(line.split("\t")[3])
            gg_stop = str(line.split("\t")[4])
            strand = str(line.split("\t")[6])
            id = line.split("\t")[8].split(';')[0].split("ID=")[1]
            ## now get the new name assigned by hamburger:
            for new_name, old_name in list_of_genes.items():    # for name, age in dictionary.iteritems():  (for Python 2.x)
                if old_name == id:
                    hamburger_assigned_id = new_name
                    continue

            if id in old_named_genes:
                model_name = old_named_genes[id] # assigning the model name - if it's there;
            else:
                model_name = "Non-model"

            if strand == "+":
                gg_strand = "forward"
                gg_direction = "1"
            if strand == "-":
                gg_strand = "reverse"
                gg_direction = "-1"

            gg_gene_number = id


            #output to a list:
            gggenes_output.append("{cluster},{number},{start},{end},{gene},{strand},{direction},{strain},{id},{new_name}".format(cluster = strain+"_cluster_"+cluster_num,number = gg_number, gene = model_name, start = gg_start, end = gg_stop, strand = gg_strand, direction = gg_direction, strain = strain, id = id, new_name  = hamburger_assigned_id))


    ## return the list

    return(gggenes_output)

def __search_single_genome__(gff_file):

    args = parseArgs()

    ##### repeat the name changing here....
    ### if t6ss flag requested then set these using the t6ss default setting
    if args.t6ss == True:
        min_genes_num = 8
        genes_gap_num = 12
        mandatory_models = T6SS_core#"/home/djwilliams/github/hamburger/models/T6SS/T6SS_core.hmm"
        accessory_models = T6SS_accessory#"/home/djwilliams/github/hamburger/models/T6SS/T6SS_accessory.hmm"
    elif args.t6ss == False:
        if args.mandatory == False:
            print("Need to provide input mandatory hmm profile(s): --mandatory")
            sys.exit()
        if args.min_genes == False:
            print("Need to set the minimum number of hits: --min_genes)")
            sys.exit()
        if args.genes_gap == False:
            print("Need to set the maximum genes gap: --genes_gap)")
            sys.exit()
        # now change the variable names from args.xxx
        min_genes_num = args.min_genes
        genes_gap_num = args.genes_gap
        mandatory_models = args.mandatory
        if args.accessory is not None: # if accessory models were given change variable name
            accessory_models = args.accessory



    mandatory_names = __read_hmms__(mandatory_models)

    if args.accessory is not None or args.t6ss == True:
        accessory_names = __read_hmms__(accessory_models)

    number_clusters_in_strain = 0
    number_rejected_clusters_in_strain = 0
    number_contig_break_clusters_in_strain = 0


    output_dir = args.output

    #print(gff_file)

    ### take the strain name (last field from "_" separator and remove the .gff from the end:)

    strain = gff_file.split('/')[-1].replace(".gff","")

    ### now made directory for this strain:
    strain_dir = "{output_dir}/{strain}".format(output_dir = output_dir, strain = strain)
    os.makedirs(strain_dir)

#3 - extract protein sequences from the gff file

### extract the protein sequences
    prot_seqs = "{output_dir}/{strain}/{strain}.faa".format(output_dir = output_dir, strain = strain)
    gff2faa_output = __gff2faa__(gff_file,prot_seqs,strain,output_dir)
    #altered_gene_names = gff2faa_output[0] # a "yes" or "no" - is only no because this is currently set - look at later??
    list_of_genes = gff2faa_output[1] # list of the new names to the names in the gff file
    genes_and_contig = gff2faa_output[2]

#4 - create TEMPORARY singlefasta files, fasta of the whole sequence, and separate annotation files

    gff_split_output = __gff_splitter__(gff_file, strain_dir)

    annotation = gff_split_output[0]
    GC_genome = gff_split_output[1]

    ##set singlefasta directory
    singlefasta_dir = "{strain_dir}/singlefastas".format(strain_dir = strain_dir)
    __multifasta_to_singlefasta__("contigs.fna", strain_dir, singlefasta_dir)

#5 - Run HMMER and filter according to the cutoff
#- for both mandatory and accessory

    mandatory_hmmer_tuple = __run_hmmer__(hmmsearch, mandatory_models, prot_seqs, args.cutoff, strain_dir, "mandatory")
    mandatory_hmmer_output = mandatory_hmmer_tuple[0]
    mandatory_hmmer_genes = mandatory_hmmer_tuple[1]


    ###### then combine the accessory and mandatory results if required.... (think about adding any genes that are not allowed in gene clusters?)
    if accessory_models:
        accessory_hmmer_tuple = __run_hmmer__(hmmsearch, accessory_models, prot_seqs, args.cutoff, strain_dir, "accessory")
        accessory_hmmer_output = accessory_hmmer_tuple[0]
        accessory_hmmer_genes = accessory_hmmer_tuple[1]

        ### add accessory to the mandatory



        total_hmmer_output = mandatory_hmmer_output + accessory_hmmer_output
        total_hmmer_genes = mandatory_hmmer_genes + accessory_hmmer_genes

    else:
        total_hmmer_output = mandatory_hmmer_output
        total_hmmer_genes = mandatory_hmmer_genes

#5 - Cluster the HMMER hits according to the parameters


    if len(mandatory_hmmer_genes) == 0:
        ##no hits
        with open(output_dir+"/log_file.txt", "a") as output:
            output.write("No hmmer hits in {strain}".format(strain =strain))
        ### and write to strain_statistics.csv
        with open(output_dir + "/"+strain+"/strain_statistics.csv","a") as output:
            output.write("{strain},{clusters},{rejected_clusters}\n".format(strain=strain,clusters = "0",rejected_clusters = "0"))
        return

    filtered_groups = __clustering_func__(total_hmmer_genes,genes_gap_num,min_genes_num)

#6 - Check that clusters are on the the same contigs


    # print(annotation) # lines of the annotation from the gff file
    # print(filtered_groups) # list of lists for each gene cluster from the (simple) algorithm
    # print(strain) # name of the strain
    # print(list_of_genes) # dictionary with the "new" gene number/identifier as the key and the "old" gene number/identifier as the value
    # print(genes_and_contig) #dictionary with the "new" gene number/identifier as the key and the contig as the value

    ##set counter for the number of the gene cluster
    counter = 0

    #### check to see if there are no gene clusters that are reported:
    #if len(filtered_groups) == 0:
        #print("No gene clusters found")
        #continue

    for group in filtered_groups: # now getting into working with each gene cluster/gene cluster
        #print(group)
        ## make name

        contig_check =  __same_contigs_check__(group, strain, genes_and_contig)
        #print(contig_check)

        if contig_check[0] == False:
            ### the hits are spread across more than one contigs

            with open(output_dir+"/log_file.txt", "a") as output:
                output.write("Contig break over gene cluster in {strain} \n\n".format(strain = strain))
                number_contig_break_clusters_in_strain += 1

            ### need to output this in some manner

            continue

        elif contig_check[0] == True: # if true, hits are on a single contig


            contig = contig_check[1] # store the name of the contig
            #update counter

            ### now check for the presence / absence of mandatory / accessory genes:

            if accessory_models: # if accessory was used - remove these hits (if present) from the file:

                mandatory_in_cluster =  [x for x in group if x not in accessory_hmmer_genes]

                ### now associate each hit with the hmmer output to get the number of queries that were satisfied from the mandatory, concatenated hmm model
                mandatory_genes_in_cluster = __gene_names_in_cluster__(mandatory_in_cluster, mandatory_hmmer_output, strain)


                #get the names of the accessory genes for later:
                accessory_in_cluster =  [x for x in group if x not in mandatory_hmmer_genes]

                ### now associate each hit with the hmmer output to get the number of queries that were satisfied from the mandatory, concatenated hmm model
                accessory_genes_in_cluster = __gene_names_in_cluster__(accessory_in_cluster, accessory_hmmer_output, strain)

                ### cancel the search if the number of unique genes is not less than or equal to the min number set in the search:
                if len(mandatory_genes_in_cluster) < int(min_genes_num):
                #    print("Not enough unique mandatory genes in the cluster found, not reporting gene cluster")
                    number_rejected_clusters_in_strain += 1

                    continue

                    #### could also still extract the region, but report that t didn't pass threshold??

            else: ### if accessory models weren't used - but still need to query whether enough unique genes were found..

                query_types = __gene_names_in_cluster__(group, mandatory_hmmer_output, strain)

                ## do the same again - are there enough genes:
                #print(number_query_types)
                if len(query_types) < int(min_genes_num):
                #    print("Not enough unique mandatory genes in the cluster found, not reporting gene cluster")

                    continue


#7 - Extract the gff subsequence - with and without the extra sequences

        #### now set a name for the gene cluster:
        counter += 1 # update counter now that the identified gene cluster has passed filtering
        number_clusters_in_strain += 1 # and same for counting clusters / strain

        cluster_name = "{strain}_cluster_{counter}".format(strain=strain, counter = str(counter))

        #### find numbers for the beginning and the end of the gene cluster:
        min_gene = list_of_genes[strain+"_"+str(min(group))]
        max_gene = list_of_genes[strain+"_"+str(max(group))]

        ## get orientation, or if single gene:
        if min_gene == max_gene:
            orientation = "same"
        elif min_gene < max_gene:
            orientation = "forward"
        elif min_gene > max_gene:
            orientation = "reverse"


        extraction_details = extract_a2b(min_gene, max_gene, annotation, strain, str(counter), args.upstream, args.downstream, contig, output_dir,orientation)

        GC_cluster = extraction_details[0]
        cluster_start = extraction_details[1]
        cluster_stop = extraction_details[2]
        extended_start = extraction_details[3]
        extended_stop = extraction_details[4]
        cluster_annotation = extraction_details[5]


        length_of_operon = cluster_start - cluster_stop
        if length_of_operon < 0:
            length_of_operon = length_of_operon * -1  ### change if the direction isn't left to write (i.s. on the reverse strand if you were to look at it in artemis)



#8 Write out the stats


        #### get lists of the numbers of genes for each hmmer query, for the genes that are present in the operon
        list_of_mandatory  = []
        for gene in mandatory_names:
            if gene not in mandatory_genes_in_cluster:
                list_of_mandatory.append("0")
            else:
                list_of_mandatory.append(str(len(mandatory_genes_in_cluster[gene])))

        #### and now for accessory
        list_of_accessory  = []
        for gene in accessory_names:
            if gene not in accessory_genes_in_cluster:
                list_of_accessory.append("0")
            else:
                list_of_accessory.append(str(len(accessory_genes_in_cluster[gene])))

        #### start creating the output:

        with open(output_dir+"/"+strain+"/cluster_stats.csv", "a") as output:
            output.write("""{strain}_cluster_{counter},{strain},{contig},{start_on_contig},{stop_on_contig},{length},{number_of_mandatory_genes},{found_number_of_mandatory_genes},{percent_of_mandatory_genes_in_query},{number_of_accessory_genes},{found_number_of_accessory_genes},{percent_of_accessory_genes_in_query},{hmm_genes},{GC_cluster},{GC_genome},{GCoperonbyGCgenome}\n""".format(strain=strain,counter=str(counter),contig=contig,start_on_contig=str(cluster_start),stop_on_contig=str(cluster_stop),length=abs(int(length_of_operon)),number_of_mandatory_genes=len(mandatory_names),found_number_of_mandatory_genes=len(mandatory_genes_in_cluster),percent_of_mandatory_genes_in_query=(float(len(mandatory_genes_in_cluster))/float(len(mandatory_names)))*100,number_of_accessory_genes=len(accessory_names),found_number_of_accessory_genes=len(accessory_genes_in_cluster),percent_of_accessory_genes_in_query=(float(len(accessory_genes_in_cluster))/float(len(accessory_names)))*100,hmm_genes=','.join(list_of_mandatory+list_of_accessory),GC_cluster=GC_cluster,GC_genome=GC_genome,GCoperonbyGCgenome=str(GC_cluster/GC_genome)))


        ####now work out and write out the input for gggenes:

        ### first, if using the accessory flag, need to create a new dictionary:

        if accessory_models:
            all_genes_in_cluster = dict(mandatory_genes_in_cluster)
            all_genes_in_cluster.update(accessory_genes_in_cluster)

        else:
            all_genes_in_cluster = mandatory_genes_in_cluster

        ### now run a function to extract the gggenes information

        gggenes_input = __extract_gggenes_info__(cluster_annotation, all_genes_in_cluster,list_of_genes, extended_start, extended_stop, strain, output_dir, str(counter))



        with open(output_dir+"/"+strain+"/gggenes_input.csv","a") as output:
            for line in gggenes_input:
                output.write(line+"\n")




    ####### record number of clusters / strain:
    with open(output_dir + "/"+strain+"/strain_statistics.csv","a") as output:
        output.write("{strain},{clusters},{rejected_clusters},{contig_break_clusters}\n".format(strain=strain,clusters = number_clusters_in_strain,rejected_clusters = number_rejected_clusters_in_strain,contig_break_clusters = number_contig_break_clusters_in_strain))


    ### remove singlefastas folder:

    os.system("rm -r {output_dir}/{strain}/singlefastas".format(output_dir=output_dir,strain=strain))
    os.system("rm  {output_dir}/{strain}/contigs.fna".format(output_dir=output_dir,strain=strain))

def __clean_up_files__(gff_file):

    args = parseArgs()
    output_dir = args.output


    strain = gff_file.split('/')[-1].replace(".gff","")

    ### now made directory for this strain:
    strain_dir = "{output_dir}/{strain}".format(output_dir = output_dir, strain = strain)

    ## now remove the files that aren't:
    ## faa file of the genome (used to find T6SS genes)
    ## gff files
    ## gggenes input
    for root, dirs, files in os.walk(strain_dir):
        for file in files:
            if file.endswith(".gff"):
                continue
            elif file == "gggenes_input.csv":
                continue
            elif file == "{strain}.faa".format(strain = strain):
                continue
            elif file == "strain_statistics.csv".format(strain = strain):
                continue

            else:
                os.remove("{output_dir}/{strain}/{file}".format(output_dir = output_dir, strain = strain, file = file))



def __search_single_genome_no_accessory__(gff_file):

    args = parseArgs()

    ##### repeat the name changing here....
    ### if t6ss flag requested then set these using the t6ss default setting
    if args.t6ss == True:
        min_genes_num = 8
        genes_gap_num = 12
        mandatory_models = T6SS_core#"/home/djwilliams/github/hamburger/models/T6SS/T6SS_core.hmm"
    elif args.t6ss == False:
        if args.mandatory == False:
            print("Need to provide input mandatory hmm profile(s): --mandatory")
            sys.exit()
        if args.min_genes == False:
            print("Need to set the minimum number of hits: --min_genes)")
            sys.exit()
        if args.genes_gap == False:
            print("Need to set the maximum genes gap: --genes_gap)")
            sys.exit()
        # now change the variable names from args.xxx
        min_genes_num = args.min_genes
        genes_gap_num = args.genes_gap
        mandatory_models = args.mandatory



    mandatory_names = __read_hmms__(mandatory_models)


    number_clusters_in_strain = 0
    number_rejected_clusters_in_strain = 0
    number_contig_break_clusters_in_strain = 0


    output_dir = args.output

    #print(gff_file)

    ### take the strain name (last field from "_" separator and remove the .gff from the end:)

    strain = gff_file.split('/')[-1].replace(".gff","")

    ### now made directory for this strain:
    strain_dir = "{output_dir}/{strain}".format(output_dir = output_dir, strain = strain)
    os.makedirs(strain_dir)

#3 - extract protein sequences from the gff file

### extract the protein sequences
    prot_seqs = "{output_dir}/{strain}/{strain}.faa".format(output_dir = output_dir, strain = strain)
    gff2faa_output = __gff2faa__(gff_file,prot_seqs,strain,output_dir)
    #altered_gene_names = gff2faa_output[0] # a "yes" or "no" - is only no because this is currently set - look at later??
    list_of_genes = gff2faa_output[1] # list of the new names to the names in the gff file
    genes_and_contig = gff2faa_output[2]

#4 - create TEMPORARY singlefasta files, fasta of the whole sequence, and separate annotation files

    gff_split_output = __gff_splitter__(gff_file, strain_dir)

    annotation = gff_split_output[0]
    GC_genome = gff_split_output[1]

    ##set singlefasta directory
    singlefasta_dir = "{strain_dir}/singlefastas".format(strain_dir = strain_dir)
    __multifasta_to_singlefasta__("contigs.fna", strain_dir, singlefasta_dir)

#5 - Run HMMER and filter according to the cutoff
#- for both mandatory and accessory

    mandatory_hmmer_tuple = __run_hmmer__(hmmsearch, mandatory_models, prot_seqs, args.cutoff, strain_dir, "mandatory")
    mandatory_hmmer_output = mandatory_hmmer_tuple[0]
    mandatory_hmmer_genes = mandatory_hmmer_tuple[1]


    ###### then combine the accessory and mandatory results if required.... (think about adding any genes that are not allowed in gene clusters?)

    total_hmmer_output = mandatory_hmmer_output
    total_hmmer_genes = mandatory_hmmer_genes

#5 - Cluster the HMMER hits according to the parameters

    ### introduce break if there are no hmmer hits at all!:

    if len(mandatory_hmmer_genes) == 0:
        ##no hits
        with open(output_dir+"/log_file.txt", "a") as output:
            output.write("No hmmer hits in {strain}".format(strain =strain))
        ### and write to strain_statistics.csv
        with open(output_dir + "/"+strain+"/strain_statistics.csv","a") as output:
            output.write("{strain},{clusters},{rejected_clusters},{contig_break_clusters}\n".format(strain=strain,clusters = "0",rejected_clusters = "0",contig_break_clusters = "0"))
        return

    filtered_groups = __clustering_func__(total_hmmer_genes,genes_gap_num,min_genes_num)

#6 - Check that clusters are on the the same contigs


    # print(annotation) # lines of the annotation from the gff file
    # print(filtered_groups) # list of lists for each gene cluster from the (simple) algorithm
    # print(strain) # name of the strain
    # print(list_of_genes) # dictionary with the "new" gene number/identifier as the key and the "old" gene number/identifier as the value
    # print(genes_and_contig) #dictionary with the "new" gene number/identifier as the key and the contig as the value

    ##set counter for the number of the gene cluster
    counter = 0

    #### check to see if there are no gene clusters that are reported:
    #if len(filtered_groups) == 0:
        #print("No gene clusters found")
        #continue

    for group in filtered_groups: # now getting into working with each gene cluster/gene cluster
        #print(group)
        ## make name

        contig_check =  __same_contigs_check__(group, strain, genes_and_contig)
        #print(contig_check)

        if contig_check[0] == False:
            ### the hits are spread across more than one contigs

            with open(output_dir+"/log_file.txt", "a") as output:
                output.write("Contig break over gene cluster in {strain} \n\n".format(strain = strain))
            ### need to output this in some manner

            continue

        elif contig_check[0] == True: # if true, hits are on a single contig


            contig = contig_check[1] # store the name of the contig
            #update counter

            ### now check for the presence / absence of mandatory / accessory genes:




            query_types = __gene_names_in_cluster__(group, mandatory_hmmer_output, strain)

            ## do the same again - are there enough genes:
            #print(number_query_types)
            if len(query_types) < int(min_genes_num):
            #    print("Not enough unique mandatory genes in the cluster found, not reporting gene cluster")

                continue


#7 - Extract the gff subsequence - with and without the extra sequences

        #### now set a name for the gene cluster:
        counter += 1 # update counter now that the identified gene cluster has passed filtering
        number_clusters_in_strain += 1 # and same for counting clusters / strain

        cluster_name = "{strain}_cluster_{counter}".format(strain=strain, counter = str(counter))

        #### find numbers for the beginning and the end of the gene cluster:
        min_gene = list_of_genes[strain+"_"+str(min(group))]
        max_gene = list_of_genes[strain+"_"+str(max(group))]


        ## get orientation, or if single gene:
        if min_gene == max_gene:
            orientation = "same"
        elif min_gene < max_gene:
            orientation = "forward"
        elif min_gene > max_gene:
            orientation = "reverse"



        extraction_details = extract_a2b(min_gene, max_gene, annotation, strain, str(counter), args.upstream, args.downstream, contig, output_dir,orientation)

        GC_cluster = extraction_details[0]
        cluster_start = extraction_details[1]
        cluster_stop = extraction_details[2]
        extended_start = extraction_details[3]
        extended_stop = extraction_details[4]
        cluster_annotation = extraction_details[5]


        length_of_operon = cluster_start - cluster_stop
        if length_of_operon < 0:
            length_of_operon = length_of_operon * -1  ### change if the direction isn't left to write (i.s. on the reverse strand if you were to look at it in artemis)



#8 Write out the stats


        ### get dictionary for gene names:
        mandatory_genes_in_cluster = __gene_names_in_cluster__(group, mandatory_hmmer_output, strain)


        #### get lists of the numbers of genes for each hmmer query, for the genes that are present in the operon
        list_of_mandatory  = []
        for gene in mandatory_names:
            if gene not in mandatory_genes_in_cluster:
                list_of_mandatory.append("0")
            else:
                list_of_mandatory.append(str(len(mandatory_genes_in_cluster[gene])))

        #### start creating the output:

        with open(output_dir+"/"+strain+"/cluster_stats.csv", "a") as output:
            output.write("""{strain}_cluster_{counter},{strain},{contig},{start_on_contig},{stop_on_contig},{length},{number_of_mandatory_genes},{found_number_of_mandatory_genes},{percent_of_mandatory_genes_in_query},{hmm_genes},{GC_cluster},{GC_genome},{GCoperonbyGCgenome}\n""".format(strain=strain,counter=str(counter),contig=contig,start_on_contig=str(cluster_start),stop_on_contig=str(cluster_stop),length=abs(int(length_of_operon)),number_of_mandatory_genes=len(mandatory_names),found_number_of_mandatory_genes=len(mandatory_genes_in_cluster),percent_of_mandatory_genes_in_query=(float(len(mandatory_genes_in_cluster))/float(len(mandatory_names)))*100,hmm_genes=','.join(list_of_mandatory),GC_cluster=GC_cluster,GC_genome=GC_genome,GCoperonbyGCgenome=str(GC_cluster/GC_genome)))


        ####now work out and write out the input for gggenes:

        ### first, if using the accessory flag, need to create a new dictionary:


        all_genes_in_cluster = mandatory_genes_in_cluster

        ### now run a function to extract the gggenes information

        gggenes_input = __extract_gggenes_info__(cluster_annotation, all_genes_in_cluster,list_of_genes, extended_start, extended_stop, strain, output_dir, str(counter))



        with open(output_dir+"/"+strain+"/gggenes_input.csv","a") as output:
            for line in gggenes_input:
                output.write(line+"\n")




    ####### record number of clusters / strain:
    with open(output_dir + "/"+strain+"/strain_statistics.csv","a") as output:
        output.write("{strain},{clusters},{rejected_clusters},{contig_break_clusters}\n".format(strain=strain,clusters = number_clusters_in_strain,rejected_clusters = number_rejected_clusters_in_strain,contig_break_clusters = str(number_contig_break_clusters_in_strain)))


    ### remove singlefastas folder:

    os.system("rm -r {output_dir}/{strain}/singlefastas".format(output_dir=output_dir,strain=strain))
    os.system("rm  {output_dir}/{strain}/contigs.fna".format(output_dir=output_dir,strain=strain))

def f(q):
    q.put([42, None, 'hello'])

#################################
def main():



    ## check if hmmsearch is installed:

    if __is_tool__("hmmsearch") == False:
        print("Please install hmmer before running hamburger")
        sys.exit()





    start = time.time()
    args = parseArgs()



    ### if t6ss flag requested then set these using the t6ss default setting
    if args.t6ss == True:
        min_genes_num = 8
        genes_gap_num = 12
        mandatory_models = T6SS_core#"/home/djwilliams/github/hamburger/models/T6SS/T6SS_core.hmm"
        accessory_models = T6SS_accessory#"/home/djwilliams/github/hamburger/models/T6SS/T6SS_accessory.hmm"
    elif args.t6ss == False:
        if args.mandatory == False:
            print("Need to provide input mandatory hmm profile(s): --mandatory")
            sys.exit()
        if args.min_genes == False:
            print("Need to set the minimum number of hits: --min_genes)")
            sys.exit()
        if args.genes_gap == False:
            print("Need to set the maximum genes gap: --genes_gap)")
            sys.exit()
        # now change the variable names from args.xxx
        min_genes_num = args.min_genes
        genes_gap_num = args.genes_gap
        mandatory_models = args.mandatory
        if args.accessory is not None: # if accessory models were given change variable name
            accessory_models = args.accessory



    ### make new output directory first:
    output_dir = args.output
    os.makedirs(output_dir)


##### steps Required

#1 - Read the HMM queries and get the names of each query
    # - for both mandatory and accessory

    mandatory_names = __read_hmms__(mandatory_models)

    if args.accessory is not None or args.t6ss == True:
        accessory_names = __read_hmms__(accessory_models)




        ### make log file and write into it:


##### other outputs :

    with open(output_dir+"/strain_statistics.csv","w") as output:
        output.write("strain,number_of_gene_clusters,number_rejected_clusters,number_contig_break_clusters\n")

    with open(output_dir+"/gggenes_input.csv","w") as output:
        output.write("operon,number,start,end,gene,strand,direction,strain,CDS_identifier,hamburger_CDS_identifier\n")

#2 for each gff file make an ouput folder and work in it:

### cycle through each of the gff files in the provided list of files:



    gff_files = args.gff


    start = time.time()


    pool = multiprocessing.Pool(processes=int(args.num_threads))



    if args.accessory is not None or args.t6ss == True:

        with open(output_dir+"/log_file.txt", "w") as output:
            output.write("""--------------------HaMBURGER--------------------
            \n
            -------HMmer Based UndeRstandinG of opERons------
            \n
            Using the following parameters as input: \n\n
            \tMinimum number of genes for cluster:  -> {min_genes}\n
            \tMaximum gap between genes in cluster: -> {genes_gap}\n
            \tUpstream region length:               -> {upstream}\n
            \tDownstream region length:             -> {downstream}\n
            \tHmmer cutoff:                         -> {cutoff}\n
            \tSearching for the following genes     ->\n
            \tMandatory genes:\n
            \t\t{mandatory_genes}\n\n
            \tAccessory genes:\n
            \t\t{accessory_genes}\n
            """.format(min_genes=min_genes_num, genes_gap=genes_gap_num, cutoff=args.cutoff, upstream=args.upstream, downstream=args.downstream, mandatory_genes='\n\t\t'.join(mandatory_names), accessory_genes='\n\t\t'.join(accessory_names)))



        with open(output_dir+"/cluster_stats.csv", "w") as output:
            output.write("gene_cluster,strain,contig,start,stop,length,number_of_mandatory_genes,found_number_of_mandatory_genes,percent_of_mandatory_genes_in_query,number_of_accessory_genes,found_number_of_accessory_genes,percent_of_accessory_genes_in_query,{hmm_genes},GC_cluster,GC_genome,GCcluster/GCgenome\n".format(hmm_genes=','.join(mandatory_names+accessory_names)))


        # multiprocess
        #result = pool.map(__search_single_genome__, gff_files)

        print("Running Hamburger - using both mandatory and accessory HMMs")
        for _ in tqdm.tqdm(pool.imap_unordered(__search_single_genome__, gff_files), total=len(gff_files)):
            pass


    else:


        with open(output_dir+"/log_file.txt", "w") as output:
            output.write("""--------------------HaMBURGER--------------------
            \n
            -------HMmer Based UndeRstandinG of opERons------
            \n
            Using the following parameters as input: \n\n
            \tMinimum number of genes for cluster:  -> {min_genes}\n
            \tMaximum gap between genes in cluster: -> {genes_gap}\n
            \tUpstream region length:               -> {upstream}\n
            \tDownstream region length:             -> {downstream}\n
            \tHmmer cutoff:                         -> {cutoff}\n
            \tSearching for the following genes     ->\n
            \tMandatory genes:\n
            \t\t{mandatory_genes}\n\n
            """.format(min_genes=min_genes_num, genes_gap=genes_gap_num, cutoff=args.cutoff, upstream=args.upstream, downstream=args.downstream, mandatory_genes='\n\t\t'.join(mandatory_names)))


        with open(output_dir+"/cluster_stats.csv", "w") as output:
            output.write("gene_cluster,strain,contig,start,stop,length,number_of_mandatory_genes,found_number_of_mandatory_genes,percent_of_mandatory_genes_in_query,{hmm_genes},GC_cluster,GC_genome,GCcluster/GCgenome\n".format(hmm_genes=','.join(mandatory_names)))


        # multiprocess
        #result = pool.map(__search_single_genome_no_accessory__, gff_files)

        print("Running Hamburger - no accessory HMMs given")
        for _ in tqdm.tqdm(pool.imap_unordered(__search_single_genome_no_accessory__, gff_files), total=len(gff_files)):
            pass


    end = time.time()


    print(end-start)



    #strain by strain for checking:
    # for file in args.gff:
    #     __search_single_genome__(file)


#10 Now group all of the different statistics files together

    os.system("cat {output_dir}/*/strain_statistics.csv >> {output_dir}/strain_statistics.csv".format(output_dir=output_dir))
    os.system("cat {output_dir}/*/gggenes_input.csv >> {output_dir}/gggenes_input.csv".format(output_dir=output_dir))
    os.system("cat {output_dir}/*/cluster_stats.csv >> {output_dir}/cluster_stats.csv".format(output_dir=output_dir))


#10 - Cleanup:

    if args.keep_files is False:
        for gff_file in gff_files:
            __clean_up_files__(gff_file)


    if args.itol == True:
    ## now convert the T6SS subtypes info for itol output:
        stats_file  = open("{output_dir}/strain_statistics.csv".format(output_dir=output_dir))
        stats_data = stats_file.readlines()
        stats_header = stats_data[0].split(',') # get the different headers required:
        for num in range(1,len(stats_header)): # go through each one and make an itol input file:
            #create list to write to file:
            to_write = []
            to_write.append("DATASET_GRADIENT\nSEPARATOR SPACE\nDATASET_LABEL {subtype}\nCOLOR #ff0000\nDATA\n".format(subtype = stats_header[num]))

            #then cycle through the subtypes file:
            for line in stats_data[1:]:
                toks = line.split(',')
                to_write.append("{strain} {number_observations}\n".format(strain = toks[0], number_observations = toks[num]))

            ## now write out:
            with open("{output_dir}/{subtype}_itol.txt".format(subtype=stats_header[num],output_dir=output_dir), "w") as output:
                for l in to_write:
                    output.write(l)



if __name__ == '__main__':
    main()
