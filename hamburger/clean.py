
import os



def clean_up_files(directory,strain):

    # args = parseArgs()
    # output_dir = args.output


    # strain = gff_file.split('/')[-1].replace(".gff","")

    ### now made directory for this strain:
    #strain_dir = "{output_dir}/{strain}".format(output_dir = output_dir, strain = strain)

    ## now remove the files that aren't:
    ## faa file of the genome (used to find T6SS genes)
    ## gff files
    ## gggenes input

    for root, dirs, files in os.walk(directory):
        #remove the singlefastas directory:
        for file in files:
            if file.endswith(".gff") and "cluster" in file:
                continue
            elif file == "cluster_stats.csv":
                continue
            elif file == "gggenes_input.csv":
                continue
            elif file == "{strain}.faa".format(strain = strain):
                continue
            elif file == "strain_statistics.csv".format(strain = strain):
                continue

            else:
                os.remove("{root}/{file}".format(root = root, file = file))
