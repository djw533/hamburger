
import os
from Bio import SeqIO




def multifasta_to_singlefasta(input_fasta, input_dir, output_dir):

    os.makedirs(output_dir)

    data = open("{input_dir}/{input_fasta}".format(input_dir = input_dir, input_fasta = input_fasta))

    for record in SeqIO.parse("{input_dir}/{input_fasta}".format(input_dir = input_dir, input_fasta = input_fasta), "fasta"):

        #print(record.seq)

        SeqIO.write(record, "{output_dir}/{name}.fna".format(output_dir = output_dir, name = record.id), "fasta")

    # lines = data.readlines()
    #
    # lines_string = ''.join(lines)
    #
    # entries_list = (''.join(lines)).split('>')
    # print(len(entries_list))
    #
    # for entry in entries_list:
    #      name = entry.split('\n')[0]
    #      ## going to try to just take the first part to the first break as the output name for the singlefasta.. (if there is more than one):
    #      if len(name.split()) > 1:
    #          name = name.split()[0]
    #      with open("{output_dir}/{name}.fna".format(output_dir = output_dir, name = name), "w") as output:
    #           output.write('>'+entry)
