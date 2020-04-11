#! /usr/bin/python3


## Script to rename the paths for the models - need to run from the directory that currently in



import os



cwd = os.getcwd()
print(cwd)
files = os.listdir(cwd)

## only run if in the hamburger directory

if cwd.endswith("/hamburger"):# and "hamburger.py" in files:
    ### now alter the file with the address for the T6SS models

    # T6SS_core = cwd+"/models/T6SS/T6SS_core.hmm"
    # T6SS_accessory = cwd+"/models/T6SS/T6SS_accessory.hmm"

    #now alter the hamburger_multithread.py file:

    data = open("scripts/hamburger.py")
    lines = data.readlines()

    # and the t6ss search file:
    data_post_search = open("scripts/post_hamburger_t6SS_search.py")
    data_post_search_lines = data_post_search.readlines()

    Rscript_file = open("scripts/filter_t6_types.R")
    Rscript_file_lines = data_post_search.readlines()




    # insert the lines as second and third element:

    lines.insert(1, 'hamburger_base_directory = "{address}"\n'.format(address = cwd))
    data_post_search_lines.insert(1, 'hamburger_base_directory = "{address}"\n'.format(address = cwd))
    Rscript_file_lines.insert(1, 'hamburger_base_directory <- "{address}"\n'.format(address = cwd))

    ##now write back out:

    with open("scripts/hamburger_check.py","w") as output:
        output.write(''.join(lines))

    with open("scripts/post_hamburger_t6SS_search_check.py","w") as output:
        output.write(''.join(lines))

    with open("scripts/filter_t6_types_check.py","w") as output:
        output.write(''.join(lines))


    print("Hamburger base directory set to scripts")
