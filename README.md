# Over_presented_sequences_From_Plant_Genome
-We can understand gene sequences as the combination of four letters: A, T, C, G. 
-A, T, C, G are the four nucleotides found in DNA. 
-Certain nucleotide sequences (cis-element) in the promoter region of certain gene could suggest the regulation of this gene's expression. 
-To better understand genes' behavior, this project is to identify all over-presented nucleotide sequences in the promoter region
 of a group of genes involved in the same biological processes. 
-Sample genome used here is maize genome. 

#Files included in the project: 
	-Over_presented_sequences.py: the python code for this project. 
	-place.fasta: file that contains all the already known nucleotide sequences. (Downloaded from PLACE database)
	-example.txt: target genes with its gene accession ID in maize. 
	-Zea_mays.sample.fa: sample genome sequences file.(Downloaded from maizeGDB, too big to upload to Github)
	-Zmays_284_6a.gene: maize gene annotation file.(Downloaded from maizeGDB, too big to upload to Github)

#Steps to use the program
	1.Download the required files to your working directory
	2.Open the Over_presented_sequences.py and change the path to your working directory at 12th line.
	3.run the program:
		•Select target species: type ‘1’ for maize.
		•Input length of promoter for analysis: I use ‘500’ as input, but other length is available.
		•Input list name that contains your target genes: type ‘example.txt’.
		•Input your significant threshold for test: I use ‘0.05’ as significant threshold value.
	4.Several options are provided sequentially after the analysis:
		•show cis-element distribution plot.
		•write file to local disk.
		•start a new analysis



