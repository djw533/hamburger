#! /usr/env python


import os
import traceback
from datetime import datetime
from Bio import SeqIO
import sys
from os import listdir


def parseArgs():
    """Parse command line arguments"""

    import argparse

    try:
        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter ,
            description="""Extract and plot homologous DNA regions with or without an associated phylogeny
Requires and uses blastn to find homologous regions
Use either -q or both -s and -e""")
        #### make mutually exclusive group
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('-q',
                '--query',
                action='store',
                help='Blastn query to search for <use this flag independently of -s and -e>')
        group.add_argument('-s',
                '--start',
                action='store',
                help='Blastn query to search from <use with -e>')
        parser.add_argument('-e',
                '--stop',
                action='store',
                help='Blastn query to search to <use with -s>')
        ### rest of the arguments as normal
        parser.add_argument('-b',
			    '--blastdb',
			    action='store',
			    help='Blastdb if previously defined - will speed up considerably! - only use blast databases created using this script')
        parser.add_argument('-g',
			    '--gff',
			    action='store',
                nargs='+',
                required=True,
			    help='Gff file(s) to search <required> Set flag')
        parser.add_argument('-c',
			    '--query_cover',
			    action='store',
                default=90,
			    help='Length of query covered in target (as percentage) - default = 90')
        parser.add_argument('-p',
			    '--perc_id',
			    action='store',
                default=90,
			    help='Percentage identity of query vs target - default = 90')
        parser.add_argument('-t',
			    '--tree',
			    action='store',
			    help='Associated phylogeny - tip labels must be identical to the gff file prefixes')
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
        parser.add_argument('-o',
			    '--output',
			    action='store',
                default="BLASTBUrger_output_{time}".format(time='_'.join('_'.join(str(datetime.now()).split('.')[0].split(':')).split())),
			    help='Output directory, default = Current date and time')
    except:
        print "An exception occurred with argument parsing. Check your provided options."
        traceback.print_exc()

    return parser.parse_args()

##### FUNCTIONS

##############################################################################################################

def __gff_to_fna__(gff_file,output_fasta,contig_dir):
    ''' contig dir is for the output of the singlefastas'''

    out = open(output_fasta, "w")
    strain = gff_file.split('.gff')[0]
    contigs = []
    list_of_names = {}
    with open(gff_file) as f:
        fasta = False
        for line in f:
            if fasta:
                if line.startswith('>'):
                    contigs.append('>'.join(line.strip().split('>')[1:])+strain) # slightly roundhouse way of creading unique identifiers for each contig so that there are no duplicate identifiers in the blastdb
                    out.write(line.strip('\n')+strain+'\n')
                    continue
                out.write(line)
                continue
            if line.startswith("##FASTA"):
                fasta = True
                continue
            if line.startswith("#"):
                continue
    out.close()

    #### now for each contig file:
    contig_sequences = {}
    for contig in contigs:
        if not os.path.exists(contig_dir):
            os.makedirs(contig_dir)
        contig_sequences[contig] = []

    with open(gff_file) as f:
        fasta = False
        for line in f:
            if fasta:
                if line.startswith('>'):
                    current_contig = '>'.join(line.strip().split('>')[1:])+strain
                contig_sequences[current_contig].append(line)
                continue
            if line.startswith("##FASTA"):
                fasta = True
                continue
            if line.startswith("#"):
                continue

    ##### write out contigs
    for contig, sequence in contig_sequences.items():

        out = open(contig_dir+'/'+strain.join(contig.split(strain)[:-1])+'.fasta', "w") # messy - but see above for why the contig is named weirdly - to stop duplicate id's in the blastdb

        for line in sequence:
            out.write(line)

        out.close()

    return contigs

##############################################################################################################

def __create_blast_db__(fasta_files,blastdb_name,output_dir):

    ###concatenate input fasta files
    os.system("cat {fasta_files} > {dir}/blastdb_input.fasta".format(fasta_files = ' '.join(fasta_files), dir = output_dir))

    ### make blastdb
    os.system("makeblastdb -in {dir}/blastdb_input.fasta -dbtype nucl -parse_seqids -out {dir}/{name} -title {name}".format(dir = output_dir, name = blastdb_name))

##############################################################################################################

def __run_blast__(query,blastdb,output_file):

    ### run blast
    os.system("blastn -task blastn -query {query} -gapopen 0 -gapextend 4 -db {blastdb} -outfmt 6 -out {output}".format(query = query, blastdb = blastdb, output = output_file))

##############################################################################################################

def __filter_blast__(blast_results,perc_id,query_cover,query,filtered_output):

    ### get length of query:
    FastaFile = open(query, 'rU')

    for rec in SeqIO.parse(FastaFile, 'fasta'):
        name = rec.id
        seq = rec.seq
        seqLen = len(rec)

    FastaFile.close()

    ###  filter blast results

    filtered_blast_results = []

    data = open(blast_results)
    lines = data.readlines()

    for line in lines:
        toks = line.strip().split('\t')

        hit_perc_id = float(toks[2]) # can have a decimal if wanted
        hit_length = float(toks[3]) # length of hit

        hit_query_cover = (hit_length / float(seqLen)) * 100  # calculated perc query cover

        if hit_perc_id >= perc_id and hit_query_cover >= query_cover:# only return blast hits that fit the criteria
            filtered_blast_results.append(line.strip())

    print 'query is {query}'.format(query = query)
    for line in filtered_blast_results:
        print line

    with open(filtered_output,"a") as output:
        for line in filtered_blast_results:
            output.write(line)

    return filtered_blast_results

##############################################################################################################

def __check_blast_hits__(blast_results, contigs, permitted_number):
    '''Checking only one blast hit in each strain gff file
    Contig input is a dictionary with each key being the
    strain name and the key being a list of contigs in this.
    Blast results is in blastn output format 6'''

    strains = contigs.keys() # all strains in the contigs dictionary

    ## get all of the contigs that were hits in the blast file
    hit_contigs = []
    for result in blast_results:
        toks = result.strip().split('\t')
        hit_contigs.append(toks[1])

    hits_per_strain = {} # to store number of hits per strain - keys will be strain, value will be int of # of hits
    for strain in strains:
        hits_per_strain[strain] = 0 # set number at 0 now

    for contig in hit_contigs:
        for key, value in contigs.items(): # key is the strain, value is the list of contigs in that strain
            if contig in value:
                hits_per_strain[key] += 1

    # print hits_per_strain

    ### finally check that no strain has more than 1 hit - currently leaving as fine for a strain to not have a hit:

    over_permitted = False # set boolean to stop script if required - perform in main script not here if so!

    over_limit = []

    for strain, number in hits_per_strain.items():
        if number > permitted_number:
            over_permitted = True
            over_limit.append(strain)

        if number == 0:
            print '{strain} has no hits'.format(strain = strain)

    return over_permitted, over_limit

##############################################################################################################

def __get_start_and_stop_from_blast_hits__(blast_results, contigs):

    start_and_stop = {} # dictionary for strain as another dictionary, and then start and stop as keys of that, value is integer of the start site

    ## get all of the contigs that were hits in the blast file

    for result in blast_results:
        toks = result.strip().split('\t')
        target_contig = toks[1]
        start = int(toks[8])
        stop = int(toks[9])

        ##cycle through the contigs to get the strain names:
        for strain, contigs_in_gff in contigs.items():
            if target_contig in contigs_in_gff:
                ### add details to dctionary
                start_and_stop[strain] = {}
                start_and_stop[strain]["start"] = start
                start_and_stop[strain]["stop"] = stop
                start_and_stop[strain]["contig"] = target_contig

    ###return dictionary:
    return start_and_stop

##############################################################################################################

def extract_a2b(start_index, stop_index, gff_input, strain, T6SS_number, contig, upstream, downstream, singlefastas_dir, output_dir):
    """ Take sequence and annotation from gff3 files between genomic point A and B, given both base indices """

    #### get cwd
    # currentdir = os.getcwd()+'/'

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

    with open(output_dir+strain+T6SS_number+".gff", "w") as output:
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
    with open(output_dir+strain+T6SS_number+".gff", "a") as output:
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
            with open(output_dir+strain+T6SS_number+".gff", "a") as output:
                output.write(linegff3)
                output.write("\n") ####### have removed the newline character "\n"

        else:
            pass

    ## -- Append ##fasta header to the annotation to indicate that the fasta sequence will be below

    with open(output_dir+strain+T6SS_number+".gff", "a") as output:
        output.write('##FASTA\n')


    ## -- Append fasta sequence to the end of the gff file

    record = SeqIO.read(singlefastas_dir+contig+'.fasta', "fasta")
    with open(output_dir+strain+T6SS_number+".gff", "a") as out:
        SeqIO.write(record[start:stop], out, "fasta")

    with open(output_dir+strain+T6SS_number+".fasta", "w") as out:
        SeqIO.write(record[start:stop], out, "fasta")

##############################################################################################################
def __gff_to_gggenes__(input_gffs,output_gggenes,containing_folder):


    open(output_gggenes, "w")
    with open(output_gggenes, "a") as output:
        output.write("number,operon,start,end,gene,strand,direction\n")

    for gff in input_gffs:

        data = open(containing_folder+gff)
        parse_lines = data.readlines()


        fasta = False
        strain = gff.split(".gff")[0]
        number = 1

        for line in parse_lines:

            if line.startswith('##FASTA'):
        		fasta = True
            if line.startswith('#'):
            	continue

            	print fasta

            if fasta == True:
            	continue
            if str(line.split("\t")[2]) == 'CDS':

                gg_start = str(line.split("\t")[3])
                gg_stop = str(line.split("\t")[4])
                gg_strand = str(line.split("\t")[6])
                if 'gene=' in line.split("\t")[8]:
                	##### assigning gg_gene for the name that will be printed in the figure
                	gg_gene = '_'.join(line.split("\t")[8].split('gene=')[1].split(';')[0].split(','))   # used join to get rid of any commas in the case there are some.... just so the that the gggenes input is always correct
                else:
                	gg_gene = 'X-Unknown'

                if gg_strand == "+":
                    gg_stran = "forward"
                    gg_direction = "1"
                if gg_strand == "-":
                    gg_stran = "reverse"
                    gg_direction = "-1"

                gg_gene_number = str(line.split("\t")[8].split('ID=')[1].split(';')[0])


                with open(output_gggenes, "a") as output:
                    output.write(str(number)+","+strain+","+gg_start+","+gg_stop+","+gg_gene+","+gg_stran+","+gg_direction+"\n")
                number +=1

################################################################################################
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


####### SCRIPT


def main():

    args = parseArgs()
    os.makedirs(args.output)

    #### check that either start and stop are set, or that query is set

    if args.query and args.stop != None:
        print "Can't have -q and -e set simultaneously"
        print "Exiting script"
        sys.exit()

    if args.start and args.stop == None:
        print 'Please provide stop query with -e flag'
        print "Exiting script"
        sys.exit()

    # if args.query == False and args.start == False and args.stop == False:
    #     print 'Need to use fasta query with flag "-q" or start and stop with flags "-s" and "-e" respectively'
    #     print 'Exiting script'
    #     sys.exit()
    #
    # if args.query != False and args.start


    ### set perc_id and query_cover as integers

    args.query_cover = int(args.query_cover)
    args.perc_id = int(args.perc_id)
    args.upstream = int(args.upstream)
    args.downstream = int(args.downstream)


    ## set db :
    if args.blastdb:
        blastdb_location = args.blastdb
    else:
        blastdb_location = args.output+'/blast_burger_db'


    #### 1. Get fasta sequences from gff file

    fasta_files = []
    contigs = {} # create dictionary of contigs to check later that eeach genome only has one hit

    for gff in args.gff:
        strain = '.'.join(gff.strip().split('.')[:-1])
        contigs[strain] = __gff_to_fna__(gff,args.output+'/'+strain+'.fna',args.output+'/'+strain)   ### extract sequences:
        fasta_files.append(args.output+'/'+strain+'.fna') # list of fasta files to pass to create blastdb

    # print contigs

    #### 2. create blastdb
    if args.blastdb == None:
        __create_blast_db__(fasta_files,"blast_burger_db",args.output)

    #### 3. run blast

    if args.query:
        __run_blast__(args.query,blastdb_location,args.output+"/single_query_blast_results.tbl")

        #### 4. Filter blast

        filtered_blast_results = __filter_blast__(args.output+"/single_query_blast_results.tbl", args.perc_id, args.query_cover, args.query, 'single_query_filtered_blast_results')

        # for line in filtered_blast_results:
        #     print line

        #### 5. Check that there is onlt one hit in each genome:

        blast_check  = __check_blast_hits__(filtered_blast_results, contigs, 1)
        over_permitted = blast_check[0]
        strains_over_limit = blast_check[1]
        if over_permitted == True:
            print '{strains} have more than one hit within the criteria.\nExiting script now.'.format(strains = ','.join(strains_over_limit))
            sys.exit()


    #### 6. Get start and stop from the filtered_blast_results

        start_and_stop = __get_start_and_stop_from_blast_hits__(filtered_blast_results, contigs)


    #### 7. Extract gff from a to b

        os.makedirs(args.output+'/extracted_gffs_and_sequences')

        for strain , details  in start_and_stop.items():

            start = details["start"]
            stop = details["stop"]
            contig = strain.join(details["contig"].split(strain)[:-1]) # remove strain name from the end of the contig

            extract_a2b(start, stop, strain+'.gff', strain, "", contig, args.upstream, args.downstream, args.output+'/'+strain+'/',args.output+'/extracted_gffs_and_sequences/')


    ##### repeat for start and stop - i.e. 2 blast queries:

    if args.start and args.stop:
        __run_blast__(args.start,blastdb_location,args.output+"/start_blast_results.tbl")
        __run_blast__(args.stop,blastdb_location,args.output+"/stop_blast_results.tbl")

        #### 4. Filter blast

        start_filtered_blast_results = __filter_blast__(args.output+"/start_blast_results.tbl", args.perc_id, args.query_cover, args.start, "start_filtered_blast_results.tbl")
        stop_filtered_blast_results = __filter_blast__(args.output+"/stop_blast_results.tbl", args.perc_id, args.query_cover, args.stop, "stop_filtered_blast_results.tbl")

        # for line in filtered_blast_results:
        #     print line

        #### 5. Check that there is onlt one hit in each genome:

        ###for start
        start_blast_check  = __check_blast_hits__(start_filtered_blast_results, contigs, 1)
        over_permitted = start_blast_check[0]
        strains_over_limit = start_blast_check[1]
        if over_permitted == True:
            print '{strains} have more than one hit within the criteria.\nExiting script now.'.format(strains = ','.join(strains_over_limit))
            sys.exit()

        ##for stop
        stop_blast_check  = __check_blast_hits__(stop_filtered_blast_results, contigs, 1)
        over_permitted = stop_blast_check[0]
        strains_over_limit = stop_blast_check[1]
        if over_permitted == True:
            print '{strains} have more than one hit within the criteria.\nExiting script now.'.format(strains = ','.join(strains_over_limit))
            sys.exit()


    #### 6. Get start and stop from the filtered_blast_results

        double_queries_dicts = {}
        double_queries_dicts["start"] = __get_start_and_stop_from_blast_hits__(start_filtered_blast_results, contigs)
        double_queries_dicts["stop"] = __get_start_and_stop_from_blast_hits__(stop_filtered_blast_results, contigs)


    #### 7. Extract gff from a to b

        os.makedirs(args.output+'/extracted_gffs_and_sequences')

        ## check that the same strains are in both dictionaries
        start_strains = double_queries_dicts["start"].keys()
        start_strains.sort()
        stop_strains = double_queries_dicts["stop"].keys()
        stop_strains.sort()

        if start_strains == stop_strains:
            print 'strain hit lists are the same'

        for strain in start_strains:# is start or stop - as created above
            print strain
            start = double_queries_dicts["start"][strain]["start"] # N.B - was concerned that this wouldn't work - but as long as the queries are selected from the same strand, then the start and stops should result in extracted regions being wholly inclusive of the two fasta queries / blast matches of them
            stop = double_queries_dicts["stop"][strain]["stop"]
            #### check that the two hits are on the same contig and only extractgff if this is the case:
            start_contig = double_queries_dicts["start"][strain]["contig"]
            stop_contig = double_queries_dicts["start"][strain]["contig"]
            if start_contig != stop_contig:
                print "hits for {strain} are not on the same contig".format(strain = strain)
                print "Not extracting any operon"
            elif start_contig == stop_contig:
                contig = strain.join(start_contig.split(strain)[:-1])
                extract_a2b(start, stop, strain+'.gff', strain, "", contig, args.upstream, args.downstream, args.output+'/'+strain+'/',args.output+'/extracted_gffs_and_sequences/')

        #####now rejoin the main script at point 8 for both a single query or an end-end query:

    #### 8. gff to gggenes input

    ##list of gffs in the extracted folder to work with:
    gffs_to_plot = []
    for file in listdir(args.output+'/extracted_gffs_and_sequences/'):
        if file.endswith('.gff'):
            gffs_to_plot.append(file)


    __gff_to_gggenes__(gffs_to_plot,args.output+"/output_gggenes.csv",args.output+'/extracted_gffs_and_sequences/')

    ###### 9. reorder_operons

    reorder_operons(args.output+"/output_gggenes.csv",args.output+"/same_direction_gggenes.csv")

    #### 9. Rscript for the gggenes input alongside the phylogeny - if given
    os.chdir(args.output)
    if args.tree:
        os.system("Rscript ~/github/hamburger/blastburger_plots.R {tree}".format(tree = args.tree))
    else:
        os.system("Rscript ~/github/hamburger/blastburger_plots_no_tree.R ")

    ##ORDER OF FUNCTIONS
    #
    # 1. Create blastdb of the fasta sequences from the gffs:
    #
    # 2. run blast on the sequences given - either provide a single sequence - or provide a start and stop sequence........
    #
    # 3. Filter the blast results given the following cutoffs:
    #     a. length that the query should match
    #     b. perc id
    #     c. print error if there is more than one hit in any given genome - i.e. report back if there is a duplication / common sequence - suggest that the user uses a longer query
    #
    # 4. Retrieve this region from a gff file
    #
    # 5. (optional?) - extract protein sequences and cluster using phmmer?
    #
    # 6. convert this to gggenes input
    #
    # 7. create plot (also with phylogeny if user provides tree (gff suffix and the tip labels must match))






if __name__ == '__main__':
    main()

### idea of how to function and order
