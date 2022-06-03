samiRa : semi automated secondary structure prediction of miRNAs

samiRa is a semi-automation for miRNA(microRNA) secondary structure prediction. A conventional method to find conserved miRNAs in any organism is to use blsatN based homology search, where the query sequence is known plant mature miRNAs and the database is the reference genome/transcriptome/EST/GSS/STS dataset of an organism. Further steps require the matching hits to be extended by 100-200 flanking and searching that region into the secondary structure prediction standalone-tools/webservers, which is a tedious job. "samiRa" perform the same task using automation and provied the predicted structures, which just requires a manual selection - whether to choose that structure or not.

Prior requirements to run samiRa
1.	query mature miRNA(known miRNAs) fasta file
2.	Reference fasta file
	This is a file, which is used to be searched by query miRNAs.
	It may be a Genome/Denovo-assembled transcriptome/EST/GSS/STS datafile.
	It must be in single-line fasta format.
	Multi-line to Single-line fasta format can be converted using a provided perl program "Multiple_to_OneLine_FASTA.pl"
		use: perl Multiple_to_OneLine_FASTA.pl <Input-fasta-filename>	<Output-fasta-filename>
3.	Linux OS
	This program is tested on Ubuntu
4.	Stand-alone mfold program
	can be downloaded from: http://www.unafold.org/mfold/software/download-mfold.php
5.	Stand-alone VARNA program
	can be downloaded from: http://varna.lri.fr/index.php?lang=en&page=downloads&css=varna
	after installation, the respective path must be set at "$VARNA_path" variable
6.	Stand-alone blast
	can be downloaded from: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/
	after installation, the respective path must be set at "$blsatN_path" variable
7.	Use makeblastdb utility of blastn to index the reference dataset
	usage: makeblastdb -in <input_fasta_1line-format> -input_type fasta -dbtype nucl
8.	samiRa usage:
	perl samiRa.pl <query-mature-miRNAs-fastafile> <Reference_file_1line-format>	<Process_keyword>
	Here, the "<Process_keyword>" is an alphabet or word used as a prefix in generated outfile.
