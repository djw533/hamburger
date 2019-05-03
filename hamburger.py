#! /usr/bin/env python

###### modify the path/name for hmmsearch below:

hmmsearch='hmmsearch'


import sys
import os
import traceback
import Bio
import shutil
import re
from Bio import SeqIO

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import translate
from Bio.Seq import reverse_complement
from shutil import copyfile
import numpy as np
from datetime import datetime
from Bio.SeqUtils import GC
import time
from datetime import timedelta

###### parse arguments

def parseArgs():
    """Parse command line arguments"""

    import argparse

    try:
        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter ,
            description="""Extract and plot operons based on hmm profiles
--------------------HaMBURGER--------------------\n
\n


-------HMmer Based UndeRstandinG of opERons------\n
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
                '--hmm',
                action='store',
                required=True,
                help='Hmm profile input <required> Set flag')
        parser.add_argument('-g',
			    '--gff',
			    action='store',
                nargs='+',
                #help='<Required> Set flag',
                required=True,
			    help='Gff file(s) to search <required> Set flag')
        parser.add_argument('-m',
			    '--min_genes',
			    action='store',
                default=4,
			    help='Minimum number of genes in operon, default = 4')
        parser.add_argument('-n',
			    '--genes_gap',
			    action='store',
                default=10,
			    help='Maximum number of genes gap between hits, default = 10')
        parser.add_argument('-u',
			    '--upstream',
			    action='store',
                default=0,
			    help='Number of nucleotides to include upstream of operon, default = 0')
        parser.add_argument('-d',
			    '--downstream',
			    action='store',
                default=0,
			    help='Number of nucleotides to include downstream of operon, default = 0')
#        parser.add_argument('-p',
#			    '--pfam',
#			    action='store',
#                default = False
#			    help='Text file of pfam domains to use for hmms - one per line')
        parser.add_argument('-o',
			    '--output',
			    action='store',
                default="HaMBURGER_output_{time}".format(time='_'.join('_'.join(str(datetime.now()).split('.')[0].split(':')).split())),
			    help='Output directory, default = Current date and time')
    except:
        print "An exception occurred with argument parsing. Check your provided options."
        traceback.print_exc()

    return parser.parse_args()

##### has nummbers function

def hasNumbers(inputString):
    return any(char.isdigit() for char in inputString)


######## gff2faa function


def gff2faa(gff_file, fasta_file):
    '''convert a gff file with the appended FASTA to
    gff_file = input gff file
    fasta_file = output file
    output: a protein coding FASTA file '''
    out = open("tmp_fasta.fa", "w")
    contigs = {}
    list_of_names = {}
    gene_names_altered = 'no'
    filename = gff_file.split('.gff')[0]
    with open(gff_file) as f:
        fasta = False
        counter = 1
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
            name = filename+'_'+str(counter)
            counter += 1
            list_of_names[name] = true_name
            if toks[0] not in contigs:
                contigs[toks[0]] = []
            contigs[toks[0]].append({"name": name, "start": int(toks[3])-1, "stop": int(toks[4]), "strand": toks[6]})
    out.close()

    ## read the contigs and save the final fasta file
    out = open(fasta_file, "w")
    with open("tmp_fasta.fa") as handle:
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
    os.remove("tmp_fasta.fa")
    return gene_names_altered,list_of_names

###################################################################################################################################


#### extracting from 2 indexes function   -  editing as if there is no upstream of downstream, then the start region is always missing one base...
def extract_a2b(start_index, stop_index, gff_input, strain, T6SS_number, contig, upstream, downstream):
    """ Take sequence and annotation from gff3 files between genomic point A and B, given both base indices """

    #### get cwd
    currentdir = os.getcwd()+'/'

    #print 'Checking to see if hit is on the positive or negative strand: -----'
    ### N.B this should always be positive - this is a bit that has been taken over from a previous script that equated blast hits (which had been reversed by blast)
    ### to gff subsections from the original annotated genome
    result = stop_index - start_index
    if result > 0:
        #print 'Hit is on the positive strand, proceeding : ------'
        start = int(start_index) - int(upstream)
        stop = int(stop_index) + int(downstream)
    elif result < 0:
        #print 'Hit is on the negative strand, redefining the start and stop indices'
        start = int(stop_index) - int(downstream)
        stop = int(start_index) + int(upstream)
    elif result == 0:
        print 'Start and Stop are the same, please seek the advice of a medical professional'

    ## -- loop through each gff file in the list made above and start reading from the 2nd line of the gff file
    ## -- create new file and insert header for a new gff file

    #for gff in gff_list:

    gff_file = open(gff_input)
    gff_file_lines = gff_file.readlines()
    gff_file.close

    with open(strain+"_"+T6SS_number+".gff", "w") as output:
        output.write("##gff-version 3")
        output.write("\n")

   # print 'Beginning extraction from ', strain


    ## -- index the line where the fasta sequence begins, and make a list of all gff features above that line
    ## -- i.e. all annotation

    for line in gff_file_lines:
        if not line.startswith("#"):
            line_to_index = line
            break

    #print line_to_index

    fasta_target = gff_file_lines.index('##FASTA\n')
    hash_index = gff_file_lines.index(line_to_index)
    gff_features = gff_file_lines[hash_index:fasta_target]

    #print "fasta line is: ", fasta_target


    ## -- Next line of gff file written detailing the sequence region, how long the region is, and what the name of the chromosome is in the associated fasta sequence
    with open(strain+"_"+T6SS_number+".gff", "a") as output:
        output.write('##sequence-region'+' '+contig+' 1 '+str(stop-(start+1)))
        output.write("\n")

    ## -- Loop through annotation - define each feature start and stop indices
    ## -- If feature start and stop sites are between vgrG1_1 and vgrG1_2 (start and stop),
    ## -- annotation line is appended to new_list as a new item
    ## -- Feature start and stop indices are corrected so as to refer to the new sequence length

    for gff_line in gff_file_lines[hash_index:(fasta_target)]:
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
            with open(strain+"_"+T6SS_number+".gff", "a") as output:
                output.write(linegff3)
                output.write("\n") ####### have removed the newline character "\n"

        else:
            pass

    ## -- Append ##fasta header to the annotation to indicate that the fasta sequence will be below

    with open(strain+"_"+T6SS_number+".gff", "a") as output:
        output.write('##FASTA\n')


    ## -- Append fasta sequence to the end of the gff file

    record = SeqIO.read(currentdir+"only_fasta/singlefastas/"+contig+'.fasta', "fasta")
    with open(strain+"_"+T6SS_number+".gff", "a") as out:
        SeqIO.write(record[start:stop], out, "fasta")

    with open(strain+"_"+T6SS_number+".fasta", "w") as out:
        SeqIO.write(record[start:stop], out, "fasta")


    ## -- Print counter to terminal screen to show progress

   # print 'Sequence and annotation extracted from', strain

###### function to check all the same

def all_same(items):
    return all(x == items[0] for x in items)

##### function to cluster numbers

def clustering_func(numbers_list,max_gap,groupsize_min):
    ### Setting maximum gap between genes in the operon

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

#### function to re-order the input for gggenes so that all operons are in the same direction

def reorder_operons(input_genes,output_genes):

    data = open(input_genes)

    lines = data.readlines()

    operons = []


    for line in lines[1:]:
        if line.split(',')[1] not in operons:
            operons.append(line.split(',')[1])

    new_lines = []

    open(output_genes, "w")

    with open(output_genes, "a") as output:
        output.write("number,operon,start,end,gene,strand,direction\n")

    for operon in operons:
        directions = []
        contig_lines = []
        for line in lines[1:]:
            if operon in line:
                directions.append(int(line.split(',')[6].strip()))
                contig_lines.append(line)

        average = float(float(sum(directions))/float(len(directions)))
        if average < 0.00:
            #### get the last number and use this to reverse everything......
            length = int(contig_lines[-1].split(',')[3])
            for new_gene in contig_lines:
                new_start = length - int(new_gene.split(',')[2])
                new_stop = length - int(new_gene.split(',')[3])
                with open(output_genes, "a") as output:
                    output.write(','.join(new_gene.split(',')[0:2])+','+str(new_start)+','+str(new_stop)+','+','.join(new_gene.split(',')[4:]))
        else:
            for new_gene in contig_lines:
                with open(output_genes, "a") as output:
                    output.write(new_gene)



def pull_out_list_of_sequences(list_file, fasta_file, output_file):
    data = open(fasta_file)
    lines = data.readlines()

    query_data = open(list_file)
    query_lines = query_data.readlines()

    entries_list = (''.join(lines)).split('>')
    hits_list = []
    for entry in entries_list:
         for query in query_lines:
              if query in entry:
                   hits_list.append(entry)

    ##### Now rejoin together :->
    open(output_file, "w")
    for hit in hits_list:
         with open(output_file, "a") as output:
              output.write('>'+hit)



########## end of funtions, start of script

def main():
    printed_columns = False

    startTime = datetime.now()
    starttime = time.time()
    args = parseArgs()

    ## Set counter for the number of files processed

    files_processed = 0

    number_of_files = len(args.gff)
    print "Number of annotated genomes to extract from: "+str(number_of_files)



    ###set up and down variables as these are changed later:
    backup_upstream = args.upstream
    backup_downstream = args.downstream


    originaldir = os.getcwd()+'/'

    gene_names = []
    #make dictionary to store the names of all the genes that fit the hmm model(s)
    hmm_genes_dict_original = {}
    hmm_data = open(args.hmm)
    hmm_data_lines = hmm_data.readlines()
    for hmmer_line in hmm_data_lines:
        if hmmer_line.startswith("NAME"):
            #add hmm name into a list in order to write out in the output list / for presence and absence of genes in the sets.......
            gene_names.append(hmmer_line.split()[1].strip())
            #finally set keys in a dictionary of each of the hmm models in the query - value is an empty list to fit all of the gene names / identifiers into
            hmm_genes_dict_original[hmmer_line.split()[1].strip()] = []


    #print gff_files

    if os.path.isdir(originaldir+args.output):
        sys.exit(args.output+" is already a directory! Exiting script")
    else:
        os.makedirs(originaldir+args.output)

    os.chdir(originaldir+args.output)
    os.makedirs("extracted_fastas")
    os.makedirs("proteins_hit")


    open("log_file.txt","w")
    with open("strain_statistics.csv","w") as output:
        output.write("strain,number_of_operons\n")

    with open("master_GGgenes.csv","w") as output:
        output.write("number,operon,start,end,gene,strand,direction\n")

    with open("operon_stats.csv", "w") as output:
        output.write("operon_name,strain,contig,start,stop,length,number_of_hmm_genes,total_number_of_hmm_genes,percent_of_genes_in_query,{hmm_genes},GC_operon,GC_genome,GCoperon/GCgenome\n".format(hmm_genes=','.join(gene_names)))


    with open("log_file.txt", "a") as output:
        output.write("""--------------------HaMBURGER--------------------
\n
-------HMmer Based UndeRstandinG of opERons------
\n
Using the following parameters as input: \n\n
\tMinimum number of genes for cluster:  -> {min_genes}\n
\tMaximum gap between genes in cluster: -> {genes_gap}\n
\tUpstream region length:               -> {upstream}\n
\tDownstream region length:             -> {downstream}\n
\tSearching for the following genes     ->\n
\t\t{genes}\n
""".format(min_genes=args.min_genes, genes_gap=args.genes_gap, upstream=args.upstream, downstream=args.downstream, genes='\n\t\t'.join(gene_names)))

    master_number = 1
    for new_input_file in args.gff:
        hmm_genes_dict = {}

        for gene in gene_names:
            hmm_genes_dict[gene] = []
            with open('proteins_hit/'+gene+'_gene_names.txt','w') as output:
                output.write("strain,CDS_identifier,parsed_identifier\n")


        base_dir = os.getcwd()+'/'
        parsing_gff_file = new_input_file.split('/')[-1]

        ### osrt out so that if absolute filepath is given then use that, otherwise, construct an absolute filepath
        if new_input_file.startswith('/') or new_input_file.startswith('~'):
            absolute_gff_file = new_input_file
        else:
            absolute_gff_file = originaldir+new_input_file


        ### if protein multifasta doesn't exist, make one

        strain = '.'.join(parsing_gff_file.split(".")[:-1])
        #print strain
        with open("log_file.txt", "a") as output:
            output.write("Extracting operons from "+strain+"\n")


        os.makedirs(base_dir+strain)


        copyfile(absolute_gff_file,base_dir+strain+"/"+parsing_gff_file )

        ###move into dir for each strain

        os.chdir(base_dir+strain)


        currentdir = os.getcwd()+'/'
        #### create dir for the extracted sequences
        os.makedirs(currentdir+"extracted_protein_sequences")

        #### add options to not have to give complete path for the hmm file
        if args.hmm.startswith('/') or args.hmm.startswith('~'):
            pass
        else:
            args.hmm = originaldir+args.hmm

        if os.path.isfile(currentdir+strain+".faa"):
            pass
        else:
            gff2faa_output_tuple = gff2faa(parsing_gff_file, strain+".faa")
            gene_names_altered = gff2faa_output_tuple[0]
            list_of_genes = gff2faa_output_tuple[1] ### dict - original name is the value - new name is the key

        #### perform hmmer
        os.system(hmmsearch+" -o log --tblout hits.tab "+args.hmm+" "+strain+".faa")


        #### should make into a function!!!
        ### extract fasta sequence and annotation from gff file
        fasta_data = open(parsing_gff_file)
        gff_lines_for_fasta = fasta_data.readlines()
        number_for_parsing_from = gff_lines_for_fasta.index("##FASTA\n")
        fasta_lines = gff_lines_for_fasta[number_for_parsing_from+1:]
        annotation_lines = gff_lines_for_fasta[:number_for_parsing_from]

        ##fasta
        if not os.path.exists(currentdir+"only_fasta"):
            os.makedirs(currentdir+"only_fasta")
            open(currentdir+"only_fasta/"+strain+".fasta","w")
        with open(currentdir+"only_fasta/"+strain+".fasta","a") as output:
            for line in fasta_lines:
                output.write(line)

        ##annotation
        if not os.path.exists(currentdir+"only_annotation"):
            os.makedirs(currentdir+"only_annotation")
            open(currentdir+"only_annotation/"+parsing_gff_file,"w")
        with open(currentdir+"only_annotation/"+parsing_gff_file,"a") as output:
            for line in annotation_lines:
                output.write(line)


        ###should also make this a function!!!!!
        ### create singlefastas to be read for each operon
        if not os.path.exists(currentdir+"only_fasta/singlefastas"):
          os.makedirs(currentdir+"only_fasta/singlefastas")

        data = open("only_fasta/"+strain+".fasta")
        lines = data.readlines()

        lines_string = ''.join(lines)
        #print lines_string
        #print GC(lines_string)
        no_entries = [x for x in lines if not x.startswith('>')]
        no_entries_string = ''.join(no_entries).replace('\n','')
        GC_genome = GC(no_entries_string)

        entries_list = (''.join(lines)).split('>')

        ### make dictionary for the contig lengths:
        contig_lengths = {}

        ##make dictionary for protein hits results:
        hits_results = {}

        for entry in entries_list:
             if entry != "":
                 name = entry.split('\n')[0].split()[0]
                 contig_lengths[name] = len(''.join(entry.split('\n')[1:]))
                 #print name
                 with open(currentdir+"only_fasta/singlefastas/"+name+'.fasta', "w") as output:
                      output.write('>'+entry)


        annotation_data = open("only_annotation/"+parsing_gff_file)
        annotation = annotation_data.readlines()

        #### parse hmmer hits

        data = open("hits.tab")

        lines = data.readlines()

        new_list = []
        protein_deets = {}

        for line in lines:
            if line.startswith("#"):
                pass
            else:
                new_list.append(line)
        #print new_list

        it_number = 0
        if len(new_list) == 0:
            print "No hits in "+strain

        else:

            ##### should really make a dictionary of the output that is used here rather than re-reading(!!!) the file thats created lower down!
            open("output.csv", "w")

            with open("output.csv", "a") as output:
                output.write("number,target name,accession,query name,accession,E-value,score,bias,E-value-D,score-D,bias-D,exp,reg,clu,ov,env,dom,rep,inc,description of target,original_name\n")


            numbers = []
            full_genes_dict = {}

            for line in new_list:
                with open("output.csv", "a") as output:
                    last_col = ' '.join(line.split()[18:])
                    last_col.replace(",", "_")
                    gene_hit= line.split()[0]
                    protein_deets[gene_hit] = {}
                    protein_deets[gene_hit]["short_number"] = line.split()[0].split('_')[-1]
                    protein_deets[gene_hit]["real_name"] = list_of_genes.get(gene_hit)
                    protein_deets[gene_hit]["hmm_query"] = line.split()[2]
                    output.write(line.split()[0].split('_')[-1]+','+','.join(line.split()[:18])+","+last_col+","+list_of_genes.get(gene_hit)+"\n")
                numbers.append(int(line.split()[0].split('_')[-1]))
                ###fixing bug for gene numbers losing the preceding 'zeroes'
                full_genes_dict[int(line.split()[0].split('_')[-1])] =  line.split()[0].strip()



            #### cluster genes

            filtered_groups = clustering_func(numbers,args.genes_gap,args.min_genes)

            open("GGgenes_input.csv", "w")
            with open("GGgenes_input.csv", "a") as output:
                output.write("number,operon,start,end,gene,strand,direction\n")


            number = 1
            #print 'filtered_groups are '
            #print filtered_groups
            for group in filtered_groups:

                contigs = []
                full_gene_list = []

                for gene in group:
                    for item in new_list:
                        #### if item.split()[0].strip().endswith(str(gene)): ### removed to use a dictionary to fix a bug
                        #### what are we trying to do here???
                        #### if the name of the hit from hmmer is the same as the long name value of the parsed hmmer hits:
                        #### gene a gene name????
                        if item.split()[0].strip() == full_genes_dict.get(gene):
                            gene = list_of_genes.get(item.split()[0].strip())

                            #####altered for NCBI demon genome files
                            if gene_names_altered == 'yes':
                                m = re.search("\d", gene)
                                gene = gene[:m.start()-1] + gene[m.start():]
                    full_gene_list.append(gene)
                    for line in annotation:
                        if len(line.split("\t")) == 9:
                            if line.split("\t")[8].split("=")[1].startswith(gene):
                                contigs.append(line.split("\t")[0].strip())
                #print 'contigs'
                #print contigs
                if all_same(contigs):
                    it_number += 1
                    #print contigs
                    #print "All on contig: " + contigs[0]
                    with open("../log_file.txt", "a") as output:
                        output.write("Group containing the following CDSs all present on: {contig}: \n{genes}\n".format(genes=','.join(map(str,full_gene_list)), contig=contigs[0]))

                    start_gene = full_gene_list[0]
                    #print 'start_gene is ' + start_gene
                    stop_gene = full_gene_list[-1]
                    #print 'stop_gene is ' + stop_gene
                    start_site = 0
                    stop_site = 0

                    for line in annotation:
                        #print line  ##### so- for some the start and stop sites are not being assigned - need to troubleshoot   - REMEMBER TO SPLIT USING A TAB --- e.g. line.split('\t')
                        if line.startswith('#'):
                            continue
                        if start_gene in line and line.split('\t')[2] == 'CDS':
                            start_site = line.split('\t')[3]
                            #print start_site
                        if stop_gene in line and line.split('\t')[2] == 'CDS':
                            stop_site = line.split('\t')[4]
                            #print stop_site

                    with open("../log_file.txt", "a") as output:
                        output.write("operon hit runs from {sgene} at position {start} to {egene} at position {stop}\n\n\n".format(sgene=start_gene, start=start_site, egene=stop_gene, stop=stop_site))

                    #print 'start is ' + str(start_site)
                    #print 'stop is ' + str(stop_site)

                    ###check contig length is enough to extract the gff correctly:
                    contig_length = contig_lengths.get(contigs[0])



                    ## modify upstream value
                    if int(stop_site) + int(args.upstream) >= contig_length:
                        args.upstream = int((contig_length - int(stop_site)) - 2 )

                    else:
                        #print 'no change upstream'
                        args.upstream = backup_upstream

                    #modify downstream value
                    if int(start_site) - int(args.downstream) <=  0:
                        args.downstream = int(start_site) - 2

                    else:
                        #print 'no change downstram'
                        args.downstream = backup_downstream



                    extract_a2b(int(start_site) - 1 , int(stop_site), parsing_gff_file, strain, "operon_"+str(it_number), contigs[0], args.upstream, args.downstream)

                    shutil.move(strain+"_"+"operon_"+str(it_number)+".fasta",base_dir+"extracted_fastas")

                    T6SS_data = open(strain+"_"+"operon_"+str(it_number)+".gff")
                    T6SS_lines = T6SS_data.readlines()
                    fasta_target = T6SS_lines.index("##FASTA\n")
                    parse_lines = T6SS_lines[2:fasta_target]

                    gene_names_in_operon = {}


                    for line in parse_lines:
                        #print 'parse_lines' ###### problem here - keeps creating entries in the list parse_lines that are just empty  - checkin the T6SS gff file? - yep - has empty lines in it
                        #print parse_lines

                        #print line + ' here'
                        if str(line.split("\t")[2]) == 'CDS':
                            #print 'CDS is here'
                            gg_start = str(line.split("\t")[3])
                            gg_stop = str(line.split("\t")[4])
                            gg_strand = str(line.split("\t")[6])

                            if gg_strand == "+":
                                gg_stran = "forward"
                                gg_direction = "1"
                            if gg_strand == "-":
                                gg_stran = "reverse"
                                gg_direction = "-1"
                            gg_gene_number = str(line.split("\t")[8].split('ID=')[1].split(';')[0])
                            #print 'gg-gene-number is ' +   gg_gene_number

                           # if '_' not in gg_gene_number:

                            #    if hasNumbers(gg_gene_number):
                             #        m = re.search("\d", gg_gene_number)
                              #       gg_gene_number = gg_gene_number[:m.start()] + '_' + gg_gene_number[m.start():]
                               #      gg_gene_number_altered = 'yes'
                                #     #print name
                                     #print new_name
                            translate_data = open("output.csv")
                            translate = translate_data.readlines()
                            #print translate
                            gg_gene = "X-unknown"
                            for info in translate[1:]:
                                #print 'info'
                                #print info
                                ### add in contig to protein_deets dictionary
                                protein_deets[info.split(',')[1]]["contig"] = contigs[0]
                                protein_deets[info.split(',')[1]]["operon"] =  strain+"_"+"operon_"+str(it_number)
                                protein_deets[info.split(',')[1]]["operon_file"] =  strain+"_"+"operon_"+str(it_number)+".gff"
                                if gg_gene_number in info.strip('\n'):
                                    #print 'is in info'

                                    gg_gene = info.split(',')[3]


                                    ### add hmm gene names to dictionary as key, with gene number as value, if value already there, add a second gene number to create a list
                                    if gg_gene in gene_names_in_operon:
                                        gene_names_in_operon[gg_gene].append(gg_gene_number)
                                        hmm_genes_dict[gg_gene].append(gg_gene_number)
                                    else:
                                        gg_gene_number_list = []
                                        gg_gene_number_list.append(gg_gene_number)
                                        gene_names_in_operon[gg_gene] = gg_gene_number_list
                                        hmm_genes_dict[gg_gene].append(gg_gene_number)


                            with open("hit_details.csv", "w") as output:
                                output.write("hamburger_name,short_name,contig,real_name,hmm_query,operon,operon_file\n")
                            with open("hit_details.csv", "a") as output:
                                for key, value in protein_deets.items():
                                    output.write("{hamburger_name},{short_name},{contig},{real_name},{hmm_query},{operon},{operon_file}\n".format(hamburger_name=key,short_name=value.get("short_name"),contig=value.get("contig"),real_name=value.get("real_name"),hmm_query=value.get("hmm_query"),operon=value.get("operon"),operon_file=value.get("operon_file")))





                            with open("GGgenes_input.csv", "a") as output:
                                output.write(str(number)+","+strain+"_operon_"+str(it_number)+","+gg_start+","+gg_stop+","+gg_gene+","+gg_stran+","+gg_direction+"\n")
                            with open("../master_GGgenes.csv","a") as output:
                                output.write(str(master_number)+","+strain+"_operon_"+str(it_number)+","+gg_start+","+gg_stop+","+gg_gene+","+gg_stran+","+gg_direction+"\n")
                            number +=1
                            master_number +=1


                    number_of_hmm_genes=len(gene_names_in_operon)
                    total_number_of_hmm_genes=len(gene_names)
                    percent_of_genes_in_query=str(float(float(len(gene_names_in_operon))/float(len(gene_names)))*100)

                    # record the number of occurences of each gene in the operon that was queried in the hmm input
                    presence_absence_string = []
                    for gene in gene_names:
                        if gene in gene_names_in_operon:
                            presence_absence_string.append(len(gene_names_in_operon[gene]))
                        else:
                            presence_absence_string.append(0)

                    ## extract fasta of operon sequence to get the GC content of the operon
                    data_for_reading_fasta_sequence = open(strain+"_operon_"+str(it_number)+".gff")
                    operon_gff_lines = data_for_reading_fasta_sequence.readlines()
                    target_for_sequence = operon_gff_lines.index("##FASTA\n")
                    operon_sequence = operon_gff_lines[target_for_sequence+1:]
                    [x for x in lines if not x.startswith('>')]
                    operon_sequence_string = ''.join([x for x in operon_sequence if not x.startswith('>')]).replace("\n","")
                    GC_operon=GC(operon_sequence_string)
                    #rint "length of operon is "+str(len(operon_sequence_string))

                    with open("../operon_stats.csv", "a") as output:
                        #output.write("operon_name,strain,contig,start,stop,length,number_of_hmm_genes,percent_of_genes_in_query,total_number_of_hmm_genes,{hmm_genes},GC,GC_genome,GCoperon/GCgenome\n".format(hmm_genes=','.join(gene_names)))
                        output.write("{strain}_operon_{it_number},{strain},{contig},{start_on_contig},{stop_on_contig},{length},{number_of_hmm_genes},{total_number_of_hmm_genes},{percent_of_genes_in_query},{hmm_genes},{GC_operon},{GC_genome},{GCoperonbyGCgenome}\n".format(strain=strain,it_number=str(it_number),contig=contigs[0],start_on_contig=start_site,stop_on_contig=stop_site,length=abs(int(stop_site)-int(start_site)),number_of_hmm_genes=number_of_hmm_genes,percent_of_genes_in_query=percent_of_genes_in_query,total_number_of_hmm_genes=total_number_of_hmm_genes,hmm_genes=','.join(map(str,presence_absence_string)),GC_operon=GC_operon,GC_genome=GC_genome,GCoperonbyGCgenome=str(GC_operon/GC_genome)))

                else:
                    #print "Not on the same contig"
                    with open("../log_file.txt", "a") as output:
                        output.write("{contig} possesses no operon".format(contig=contigs[0]))

                    pass

            ##### Perhaps make a function??

                #print "Number of operons Extracted = "+ str(it_number)

            for key, value in hmm_genes_dict.items():
                with open(strain+'_'+str(key)+'_gene_names.txt','a') as output:
                    for item in value:
                        for new_name, old_name in list_of_genes.items():    # for name, age in dictionary.iteritems():  (for Python 2.x)
                            if old_name == item:
                                #print(new_name)
                                output.write(new_name+'\n')

                pull_out_list_of_sequences(strain+'_'+str(key)+'_gene_names.txt', strain+'.faa', strain+'_'+str(key)+'_sequences.faa')

                #move output
                shutil.move(strain+'_'+str(key)+'_gene_names.txt', currentdir+"extracted_protein_sequences/"+strain+'_'+str(key)+'_gene_names.txt')
                shutil.move(strain+'_'+str(key)+'_sequences.faa', currentdir+"extracted_protein_sequences/"+strain+'_'+str(key)+'_sequences.faa')


            #os.system("Rscript ~/GGgenes_T6SS.R")
        os.remove(base_dir+strain+"/"+parsing_gff_file)
        shutil.rmtree(base_dir+strain+"/only_fasta")
        shutil.rmtree(base_dir+strain+"/only_annotation")
        os.chdir(base_dir)
        #print 'Time elasped: '+str(datetime.now() - startTime)

        with open("strain_statistics.csv", "a") as output:
            output.write("{strain},{number_of_operons}\n".format(strain=strain,number_of_operons=str(it_number)))

        reorder_operons("master_GGgenes.csv","same_direction_gggenes_input.csv")


        ##### print out the dictionary that has the lists of each of the genes that fit each of the queried hmm models

        for key, value in hmm_genes_dict.items():
            #print key
           # print value
            with open('proteins_hit/'+str(key)+'_gene_names.txt','a') as output:
                for item in value:
                    for new_name, old_name in list_of_genes.items():
                        if old_name == item:
                            output.write(strain+','+item+','+new_name+'\n')

            #shutil.move(str(key)+'_gene_names.txt', 'proteins_hit/'+str(key)+'_gene_names.txt')

        hmm_genes_dict.clear()




        ###### update counter and print progress to screen
        files_processed += 1
        percent_done = float(float(files_processed) / float(number_of_files)) * 100
        current_time = time.time()
        time_taken = (current_time - starttime) / 60
        time_left = (float(time_taken) / percent_done) * ( 100.00 - percent_done)

        if printed_columns == False:
            print("{: <20} {: <15} {: <30} {: <30}".format("Percent completed", "Time taken", "Predicted time remanining", "Predicted finish time"))
            print("{: <20} {: <15} {: <30} {: <30}".format(str("{0:.2f}".format(percent_done))+"%", str("{0:.2f}".format(time_taken))+"mins", str("{0:.2f}".format(time_left))+"mins", str(datetime.now() + timedelta(minutes=time_left))))
            printed_columns = True

        if files_processed%25 == 0:
            print("{: <20} {: <15} {: <30} {: <30}".format(str("{0:.2f}".format(percent_done))+"%", str("{0:.2f}".format(time_taken))+"mins", str("{0:.2f}".format(time_left))+"mins", str(datetime.now() + timedelta(minutes=time_left))))

            #print "Completed "+str("{0:.2f}".format(percent_done))+"%"
            #current_time = time.time()
            #time_taken = (current_time - starttime) / 60
            #time_left = (float(time_taken) / percent_done) * ( 100.00 - percent_done)
            #print 'Time taken : ' + str("{0:.2f}".format(time_taken))+" minutes"
            #print 'Predicted time remaining : '+ str("{0:.2f}".format(time_left))+" minutes"
            #print 'Predicted finish time : ' + str(datetime.now() + timedelta(minutes=time_left))



        if percent_done == 100:
            print("{: <20} {: <15} {: <30} {: <30}".format(str("{0:.2f}".format(percent_done))+"%", str("{0:.2f}".format(time_taken))+"mins", str("{0:.2f}".format(time_left))+"mins", str(datetime.now() + timedelta(minutes=time_left))))
            print "Completed 100%"
            current_time = time.time()
            time_taken = (current_time - starttime) / 60
            print 'Time taken : ' + str("{0:.2f}".format(time_taken))+" minutes"

        #os.system("Rscript ~/master_GGgenes_T6SS.R")

if __name__ == '__main__':
    main()
