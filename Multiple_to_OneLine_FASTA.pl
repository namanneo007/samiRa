#	CREATED BY NAMAN MANGUKIA
#	THIS PROGRAM CONVERT MULTILINE FASTA INTO SINGLE LINE FASTA
#
#	EX.
#	seq1.fasta
#	
#	>s1
#	ACGATGC
#	CA
#	CAATA
#	>s2
#	CAGA
#	AC
#	>s3
#	GCATA
#	CAACACA
#	
#	perl Multiple_to_OneLine_FASTA.pl seq1.fasta out1.fasta
#	
#	out1.fasta
#	
#	>s1
#	ACGATGCCACAATA
#	>s2
#	CAGAAC
#	>s3
#	GCATACAACACA


$arg1=$#ARGV+1;
if($arg1!=2)
{
	print "USAGE : perl	$0	<input-FASTA-file>	<output-FileName>\n";
	exit 0;
}

open(IN,"$ARGV[0]");
open(OUT,">$ARGV[1]");
while (<IN>)
{
	if (/>/)
	{	chomp;$_=~s/[\r]//g;
		if ($. == 1)
		{
			print OUT "$_\n";
		}
		else {print OUT "\n$_\n";}
	}
	else
	{
		chomp;$_=~s/[ \t\n\r]//g;
		print OUT $_;
	}
}
close IN;


