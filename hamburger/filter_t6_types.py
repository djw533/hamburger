

from hamburger import tool_check
from hamburger import istarmap

import time, os, sys


import subprocess
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

#import multiprocessing
from multiprocessing import Pool
import tqdm
from functools import partial

#import Rscripts
from pathlib import Path
import importlib.resources as pkg_resources

def resource_path(*parts: str) -> Path:
    """
    Return a Path inside the installed hamburger package.
    e.g. resource_path("r_scripts", "check_R_dependencies.R")
    """
    return pkg_resources.files("hamburger").joinpath(*parts)



#hamburger_base_directory = os.path.abspath(__file__).split("/hamburger/filter_t6_types.py")[0]


def get_tssB_tssC(hamburger_dir):

    ## read in gggenes:
    data = pd.read_csv("{d}/gggenes_input.csv".format(
        d = hamburger_dir
    ))


    #create output_files for tssB and tssC sequences

    open(hamburger_dir + "/all_observed_tssB.faa","w")
    open(hamburger_dir + "/all_observed_tssC.faa","w")


    #get unique operons
    clusters = set(data["operon"])

    for cluster in clusters:

        #filter the pandas dataset
        temp_df = data[data["operon"] == cluster]

        #find tssB and tssC
        tssB = temp_df[temp_df["gene"] == "TssB"]
        tssC = temp_df[temp_df["gene"] == "TssC"]

        #check length
        if len(tssB) == 1 and len(tssC) == 1:
            #genes found - write out sequences
            with open(hamburger_dir + "/all_observed_tssB.faa","a") as output:
                output.write(">{id}\n{seq}\n".format(
                    id = tssB["operon"].item(),
                    seq = tssB["prot_seq"].item()))
            with open(hamburger_dir + "/all_observed_tssC.faa","a") as output:
                output.write(">{id}\n{seq}\n".format(
                    id = tssC["operon"].item(),
                    seq = tssC["prot_seq"].item()))

        elif len(tssB) == 0 or len(tssC) == 0:
            print("Not the genes required in {cluster}".format(cluster = cluster))
            continue


        else:
            #find the tssB and tssC that are adjacent
            tssB["gene_nums"] = [int(g.split("_")[-1]) for g in tssB["hamburger_CDS_identifier"]]
            tssC["gene_nums"] = [int(g.split("_")[-1]) for g in tssC["hamburger_CDS_identifier"]]

            # now look for the two closest numbers: (i.e. tssB and tssC genes that are adjacent)
            closest_genes = None
            chosen_TssB = None
            chosen_TssC = None
            for tssB_num in tssB["gene_nums"]:
                for tssC_num in tssC["gene_nums"]:
                    distance = abs(tssB_num - tssC_num)
                    if closest_genes == None:
                        closest_genes = distance
                        chosen_TssB = tssB_num
                        chosen_TssC = tssC_num
                    elif distance < closest_genes:
                        closest_genes = distance
                        chosen_TssB = tssB_num
                        chosen_TssC = tssC_num

            # get the selected sequences
            tssB = tssB[tssB["gene_nums"] == chosen_TssB]
            tssC = tssC[tssC["gene_nums"] == chosen_TssC]

            #write out
            with open(hamburger_dir + "/all_observed_tssB.faa","a") as output:
                output.write(">{id}\n{seq}\n".format(
                    id = tssB["operon"].item(),
                    seq = tssB["prot_seq"].item()))
            with open(hamburger_dir + "/all_observed_tssC.faa","a") as output:
                output.write(">{id}\n{seq}\n".format(
                    id = tssC["operon"].item(),
                    seq = tssC["prot_seq"].item()))


    return([hamburger_dir + "/all_observed_tssB.faa",hamburger_dir + "/all_observed_tssC.faa"])



def concatenate_alignments(aln_1, aln_2, concatenated_output_aln): # must have the same names for this to work

    alignment_1 = open(aln_1).readlines()

    alignment_2 = open(aln_2).readlines()



    sequences = {}


    for line in alignment_1:
        new_sequence = False

        if line.startswith(">"):
            new_sequence = True

        if new_sequence == True:
            identifier = line.strip().split('>')[1]
            if identifier not in sequences:
                sequences[identifier] = []

        if new_sequence == False:
            sequences[identifier].append(line)

    ## now do same with second file:

    for line in alignment_2:
        new_sequence = False

        if line.startswith(">"):
            new_sequence = True

        if new_sequence == True:
            identifier = line.strip().split('>')[1]
            if identifier not in sequences:
                sequences[identifier] = []

        if new_sequence == False:
            sequences[identifier].append(line)


    ### now write out:sys.argv[3],"w"

    open(concatenated_output_aln,"w")

    for identifier, sequence in sequences.items():
        with open(concatenated_output_aln,"a") as output:
            #write header
            output.write('>'+identifier+'\n')
            #write sequence:
            output.write(''.join(sequence))


def cat_files(file1, file2, output):

    os.system("cat " + str(file1) + " " + str(file2) + " > " + output)


def get_muscle_version() -> str:
    """
    Run `muscle -version` (or `--version`) and return the version string.
    """
    for flag in ("-version", "--version"):
        try:
            # Capture both stdout/stderr since some versions print to stderr
            completed = subprocess.run(
                ["muscle", flag],
                capture_output=True,
                text=True,
                check=True
            )
            output = completed.stdout or completed.stderr
            # The first line usually contains the version
            return output.strip().splitlines()[0].split()[1].strip("v").split(".")[0]
        except (subprocess.CalledProcessError, FileNotFoundError):
            continue
    raise RuntimeError("Could not invoke 'muscle' to get version; is it installed and on your PATH?")


def align_sequences(fasta_file, output_file): # cat tssB/tssC sequences in reference and those found by hamburger

    #check version of muscle:

    muscle_vers = get_muscle_version()

    #os.system("muscle -in {output_dir}/all_observed_{gene_name}.faa -out {output_dir}/all_observed_{gene_name}_aligned.fasta".format(output_dir=args.input_dir, gene_name=gene_name))

    if muscle_vers == "3":

        subprocess.run(
            args = [
            "muscle",
            "-in",
            fasta_file,
            "-out",
            output_file#,
            #"-quiet"
            ],
            #shell=True,
            #check=True,
            #stdout=subprocess.DEVNULL,
            #stderr=subprocess.DEVNULL
        )

    if muscle_vers == "5":

        subprocess.run(
            args = [
            "muscle",
            "-align",
            fasta_file,
            "-output",
            output_file#,
            #"-quiet"
            ],
            #shell=True,
            #check=True,
            #stdout=subprocess.DEVNULL,
            #stderr=subprocess.DEVNULL
        )





def draw_tree(input_alignment, output_tree):

    subprocess.run(
        args = [
        "fasttree",
        "-out",
        output_tree,
        #"-quiet",
        "-nopr",
        input_alignment
        ]#,
        # shell=True,
        # check=True,
        # stdout=subprocess.DEVNULL,
        # stderr=subprocess.STDOUT
    )


def run_filter(hamburger_dir,threads = 2, itol = False, keep_files = False):

    start = time.time()


    ## check that programs are installed:
    if tool_check.is_tool("muscle") == False:
        print("Please install muscle before running hamburger. Exiting script.")
        sys.exit()

    if tool_check.is_tool("fasttree") == False and tool_check.is_tool("FastTree") == False:
        print("Please install fasttree before running hamburger. Exiting script.")
        sys.exit()

    elif tool_check.is_tool("fasttree") == True:
        tree_software = "fasttree"

    elif tool_check.is_tool("FastTree") == True:
        tree_software = "FastTree"

    # get the Path to the script inside the installed package
    check_script = resource_path("r_scripts", 'check_R_dependencies.R')

    # call it, passing base_dir as the first and only argument
    subprocess.run(
        ["Rscript", str(check_script)],
        check=True
    )

    ## check the dependencies_check
    dependencies_check = open("dependencies_check.csv")
    d = dependencies_check.readlines()
    dependencies_check.close()
    if len(d) > 1: # then some R packages aren't installed:
        print("The following R libraries are not installed, please install them before running this script. Dependencies can be installed using the install_R_packages.R script in the hamburger/scripts directory.")
        for l in d[1:]:
            print(l.strip())
        print("Exiting script.")
        os.remove("dependencies_check.csv")
        sys.exit()
    ## also remove here:
    os.remove("dependencies_check.csv")


    # check to see if no T6SSs were found:
    if os.path.isfile("{dir}/strain_statistics.csv".format(dir = hamburger_dir)):
        strain_stats = open("{dir}/strain_statistics.csv".format(dir = hamburger_dir))
        strain_data = strain_stats.readlines()
        strain_stats.close()
    else:
        sys.exit("No strain_statistics.csv file in {dir}".format(dir = hamburger_dir))

    #now check:
    num_T6SSs_found = 0
    for line in strain_data[1:]:
        num_T6SSs_found += int(line.split(',')[1]) # add the number of T6SSs found in each strain statistics file

    if num_T6SSs_found == 0:
        sys.exit("No T6SSs were found in the previous hamburger run. Therefore any T6SS typing is not possible. Exiting.")

    #set threads
    if int(threads) > 1:
        num_threads_for_muscle = 2

    else:
        num_threads_for_muscle = 1


    #get the TssB and TssC sequences
    sequence_files = get_tssB_tssC(hamburger_dir)

    #concatenate and align the reference sequences
    print("Adding reference TssB and TssC sequences and aligning")
    #align them - for tssB
    cat_files(
        sequence_files[0],
        resource_path("t6ss_reference_set", "tssB.fasta"),
        "{output_dir}/all_tssB.fasta".format(output_dir = hamburger_dir)
    )
    align_sequences(
        fasta_file = "{output_dir}/all_tssB.fasta".format(output_dir = hamburger_dir),
        output_file = "{output_dir}/all_tssB.aln".format(output_dir = hamburger_dir)
    )

    #and for tssC
    cat_files(
        sequence_files[1],
        resource_path("t6ss_reference_set", "tssC.fasta"),
        "{output_dir}/all_tssC.fasta".format(output_dir = hamburger_dir)
    )
    align_sequences(
        fasta_file = "{output_dir}/all_tssC.fasta".format(output_dir = hamburger_dir),
        output_file = "{output_dir}/all_tssC.aln".format(output_dir = hamburger_dir)
    )

    #cat alignments
    concatenate_alignments(
        "{output_dir}/all_tssB.aln".format(output_dir = hamburger_dir),
        "{output_dir}/all_tssC.aln".format(output_dir = hamburger_dir),
        "{output_dir}/tssBC_alignment.aln".format(output_dir = hamburger_dir)
    )

    # draw tree
    print("Drawing TssBC tree")
    draw_tree(
        "{output_dir}/tssBC_alignment.aln".format(output_dir = hamburger_dir),
        "{output_dir}/tssBC_alignment.fasta.treefile".format(output_dir = hamburger_dir)
    )

    #run the Rscript
    os.chdir(hamburger_dir)
    print("Filtering T6SS subtypes according to the TssBC phylogeny")

    filter_script = resource_path("r_scripts", "filter_t6_types.R")
    # If filter_t6_types needs extra args (e.g. node_to_root), append them here
    subprocess.run(
    ["Rscript", str(filter_script), pkg_resources.files("hamburger")],
    check=True
    )


    #os.system("Rscript {hamburger_base_directory}/scripts/filter_t6_types.R {hamburger_base_directory}".format(hamburger_base_directory=hamburger_base_directory))


    if itol == True:
    ## now convert the T6SS subtypes info for itol output:
        subtypes_file  = open("T6SS_cluster_types.csv")
        subtypes_data = subtypes_file.readlines()
        subtypes_header = subtypes_data[0].split(',') # get the different headers required:
        for num in range(1,len(subtypes_header)): # go through each one and make an itol input file:
            #create list to write to file:
            to_write = []
            to_write.append("DATASET_GRADIENT\nSEPARATOR SPACE\nDATASET_LABEL {subtype}\nCOLOR #ff0000\nDATA\n".format(subtype = subtypes_header[num]))

            #then cycle through the subtypes file:
            for line in subtypes_data[1:]:
                toks = line.split(',')
                to_write.append("{strain} {number_observations}\n".format(strain = toks[0], number_observations = toks[num]))

            ## now write out:
            with open("{subtype}_itol.txt".format(subtype=subtypes_header[num]), "w") as output:
                for l in to_write:
                    output.write(l)


    # if keep_files is False:
    #     os.system("rm */*_tssC.faa */*_tssB.faa all_observed_tssC.faa all_observed_tssB.faa all_observed_tssB_aligned.fasta all_observed_tssC_aligned.fasta tssBC_alignment.fasta tssBC_alignment.fasta.treefile")




    end = time.time()


    print(end-start)
