
import os


def read_hmms(input_hmm_file):
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



def run_hmmer(hmmsearch_command, hmm_queries, input, hmmer_cutoff, output_dir, name):
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
