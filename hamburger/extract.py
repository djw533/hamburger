
from Bio import SeqIO
from Bio.SeqUtils import GC

from Bio.Seq import translate
from Bio.Seq import reverse_complement


def extract_gggenes_info(cluster_annotation, named_genes,list_of_genes,start,stop,strain,output_dir,cluster_num,protein = True):

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
            prot_seq = str(line.split("\t")[9])


            if "ID=" in line.split("\t")[8]:
                id = line.split("\t")[8].split("ID=")[1].split(";")[0]
            elif "locus_tag=" in line.split("\t")[8]:
                id = line.split("\t")[8].split("locus_tag=")[1].split(";")[0]
            else:
                id = "No_name_set_{number}".format(number = str(counter))

            #id = line.split("\t")[8].split(';')[0].split("ID=")[1]
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
            if protein == False:
                append_line = "{cluster},{number},{start},{end},{gene},{strand},{direction},{strain},{id},{new_name}".format(
                    cluster = strain+"_cluster_"+cluster_num,
                    number = gg_number,
                    gene = model_name,
                    start = gg_start,
                    end = gg_stop,
                    strand = gg_strand,
                    direction = gg_direction,
                    strain = strain,
                    id = id,
                    new_name  = hamburger_assigned_id
                 )

            elif protein == True:
                append_line = "{cluster},{number},{start},{end},{gene},{strand},{direction},{strain},{id},{new_name},{protein}".format(
                    cluster = strain+"_cluster_"+cluster_num,
                    number = gg_number,
                    gene = model_name,
                    start = gg_start,
                    end = gg_stop,
                    strand = gg_strand,
                    direction = gg_direction,
                    strain = strain,
                    id = id,
                    new_name  = hamburger_assigned_id,
                    protein = prot_seq.strip()
                 )


            gggenes_output.append(append_line)

    ## return the list

    return(gggenes_output)



def extract_a2b(start_gene, stop_gene, annotation, strain, cluster_num, upstream, downstream, contig, output_dir, orientation):
    """ Take sequence and annotation from gff3 files between genomic point A and B, given both base indices
    Orientation should be either 'same', 'forward', or 'reverse' """

    print(start_gene)
    print(stop_gene)


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
                feature_start = int(gff_line.split('\t')[3])
                feature_stop = int(gff_line.split('\t')[4])
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
                feature_start = int(gff_line.split('\t')[3])
                feature_stop = int(gff_line.split('\t')[4])
                if orientation == "forward":
                    start_index = min(feature_start, feature_stop)
                elif orientation == "reverse":
                    start_index = max(feature_start, feature_stop)
                start_got = True

            elif id == stop_gene:
                feature_start = int(gff_line.split('\t')[3])
                feature_stop = int(gff_line.split('\t')[4])
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

            #extract the protein sequence as well:
            strand = gff_line.split("\t")[6]
            if strand == "+":
                prot_seq = str(translate(record.seq[int(feature_start) - 1:int(feature_stop)]))
            elif strand == "-":
                prot_seq = str(translate(reverse_complement(record.seq[int(feature_start) - 1 :int(feature_stop)])))



            with open("{output_dir}/{strain}/{strain}_cluster_{cluster_num}.gff".format(output_dir = output_dir, strain = strain, cluster_num = cluster_num), "a") as output:
                output.write(linegff3)
                cluster_annotation.append(linegff3.strip("\n") + "\t" + prot_seq + "\n") ### append only the cluster annotation for later gggenes input csv creation
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
