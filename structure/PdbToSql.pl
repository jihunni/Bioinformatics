#!/usr/bin/perl -w

#############################################################
# PdbToSql.pl
# Translates PDB files to FASTA format sequence files
# Minimal error-checking
# Handles multi-chain files
# Author: Stefan Larson
#
# This is a simple Perl script which is designed to read in PDB files and then output their sequence in one letter code FASTA format files. 
# Simply type in "./PdbToSql.pl" at the prompt on a Unix machine, and it will ask for the name of the pdb file, outputting the sequence. 
# Users of Mac OS 9 and below may need to visit MacPerl.com to use this, or simply install it in your home directory on any unix machine 
# (like the Stanford servers).
#############################################################


# set up a hash to translate 3-letter 
# amino acid to 1-letter code
%aa=qw(
	ALA 	A
	CYS 	C
	ASP 	D
	GLU 	E 
	PHE 	F
	GLY	G
	HIS	H
	ILE	I
	LYS	K
	LEU	L
	MET	M
	ASN	N
	PRO	P
	GLN	Q
	ARG	R
	SER	S
	THR	T
	VAL	V
	TRP	W
	TYR	Y);

# open input PDB file and output sequence file
print "Which .pdb file?\n";
chomp($file=<STDIN>);
open (FILE, $file);
open(OUT, ">$file.sql");

$oldresno=-1;	# ATOMS belong to existing residue?
$chain=1;		# in case PDB holds multiple protein chains
@1=();		# numbered arrays hold the chain sequence

while (<FILE>)		# read each line of the PDB file
{
	if (/^ATOM/) 	# look for lines that start with ATOM
	{
		@line=split; 				# parse fields into @line array
		chomp($res=$line[3]); 			# residue is 4th field

		if ($line[4]=~/[A-Z]/) 			# contains chain identifier as 5th field
		{
			chomp($resno=$line[5]);		# $resno is the residue number
		}
		else						# does not contain chain identifier
		{
			chomp($resno=$line[4]);		# as above
		}
		chomp($residue=$aa{$res});		# translate 3-letter to 1-letter
		
		if ($resno>$oldresno)	
		{						# if it's a new residue ...
			push(@$chain, $residue);
		}						# ... add it to the chain sequence
		$oldresno=$resno;
	}
	if (/^TER/)		
	{							# TER signals end of chain ...
		$chain++;
		@$chain=();
		$oldresno=-1;
	}							# ... store chain sequence and start a new one
} # end of reading loop

# finished reading PDB; output all chain sequences
for $i (1..$chain-1)	
{
	print OUT ">$file.$i\n";
	foreach $residue (@$i)
	{
		print OUT "$residue";
	}
	print OUT "\n";
}			
