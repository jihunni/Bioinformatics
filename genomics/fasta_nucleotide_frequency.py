# Ref : https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/fasta_n/
# Requirement : Biopython

#Open a FASTA input file of gene nucleotides sequences:
input_file = open('/home/data/ref_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa', 'r')

#Open an output file to record the counts in.
#tsv is short for "Tab Separated Variables",
#also known as "Tab Delimited Format".
#
#This is a universal format, you can read it
#with any text editor - Microsoft Excel is
#also a good choice.
output_file = open('nucleotide_counts.tsv','w')

#We will now write a header line to our output file.
#
#We must write \t to mean a tab, and \n to mean
#an end of line (new line) character.
#
#i.e.
#Gene (tab) A (tab) C (tab) G (tab) T (tab) Length (tab) CG%
output_file.write('Gene\tA\tC\tG\tT\tLength\tCG%\n')

#We are going to need BioPython's SeqIO library, so we
#must tell Python to load this ready for us:
from Bio import SeqIO

#Get SeqIO to read this file in "fasta" format,
#and use it to see each record in the file one-by-one
for cur_record in SeqIO.parse(input_file, "fasta") :
    #Because we used the Bio.SeqIO parser, each record
    #is SeqRecord object which includes name and seq
    #properties.
    gene_name = cur_record.name

    #Just like a string in python, a Biopython sequence
    #object has a 'count' method we can use:
    A_count = cur_record.seq.count('A')
    C_count = cur_record.seq.count('C')
    G_count = cur_record.seq.count('G')
    T_count = cur_record.seq.count('T')

    #We would also like to know the number of nucleotides
    #in this gene (which should add up to the four
    #base counts, if there are no unknown bases, N)
    length = len(cur_record.seq)

    #Now work out the CG percentage for this gene.
    #We must switch from integers into floating point
    #(non-integer) because integer division will
    #just give 0 or 1 as the answer
    cg_percentage = float(C_count + G_count) / length

    #Finally, we are going to save this information
    #as a single tab separated line in our output file.
    #
    #As before (when we wrote the header line), we must
    #write \t to mean a tab, and \n to mean an end of line
    #(new line) character.
    #
    #We are using the string formatting (or interpolation)
    # operator % so the %s means insert a string,
    #while %i means insert an integer and %f a floating
    #point (non-integer).
    #
    #The \ character means this command continues on the
    #next line
    output_line = '%s\t%i\t%i\t%i\t%i\t%i\t%f\n' % \
    (gene_name, A_count, C_count, G_count, T_count, length, cg_percentage)

    # 
    output_file.write(output_line)

#Now we have finished all the genes, we can close the output file:
output_file.close()

#and close the input file:
input_file.close()
