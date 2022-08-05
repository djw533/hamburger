import subprocess



def run_prodigal(input_fasta,output_gff):
    #run prodigal and pipe output into a list, then write out into a gff file with the initial fasta sequence concatenated

    gff_data = subprocess.run(["prodigal", "-i", input_fasta, "-f", "gff", "-q"], stdout=subprocess.PIPE)
    gff_lines = gff_data.stdout.decode('utf-8').splitlines()

    gff_lines.append("##FASTA")

    fasta_data = open(input_fasta)
    fasta_lines = fasta_data.readlines()
    fasta_data.close()

    for line in fasta_lines:
        gff_lines.append(line.strip())

    #write out:
    with open(output_gff,"w") as output:
        for line in gff_lines:
            output.write(line+"\n")
