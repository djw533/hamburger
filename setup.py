#! /usr/bin/python3


## Script to rename the paths for the models - need to run from the directory that currently in



import os



cwd = os.getcwd()

files = os.listdir(cwd)

## only run if in the hamburger directory

if cwd.endswith("hamburger") and "hamburger_multithread.py" in files:
    ### now alter the file with the address for the T6SS models

    T6SS_core = cwd+"/models/T6SS/T6SS_core.hmm"
    T6SS_accessory = cwd+"/models/T6SS/T6SS_accessory.hmm"

    #now alter the hamburger_multithread.py file:

    data = open("hamburger_multithread.py")
    lines = data.readlines()


    # insert the lines as second and third element:

    lines.insert(1, 'T6SS_core = "{address}"\n'.format(address = T6SS_core))
    lines.insert(2, 'T6SS_accessory = "{address}"\n'.format(address = T6SS_accessory))


    print("Found T6SS models")

    ##now write back out:

    with open("hamburger_multithread_w_names.py","w") as output:
        output.write(''.join(lines))

    print("T6SS models known to hamburger; inserted into top of script")
