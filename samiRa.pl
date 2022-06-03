#	SCRIPT BY NAMAN MANGUKIA

#	HERE, we are assuming that, MAKEBLASTDB for Reference/database is completed

#	WE ARE ALSO ASSUMING THAT , 1-LINE FASTA IS ALREADY CREATED


$arg1=$#ARGV+1;

if($arg1!=3)
{
	print "USAGE : perl	$0	<mature-miRNA-fasta>	<Genome-fasta-1LineFormat>	<process-counter>\n";
	exit 0;
}

$t1=time();

($query1,$genome1)='';

$query1=$ARGV[0];
chomp($query1);
$genome1=$ARGV[1];
chomp($genome1);


##############	HASH_GENERATION FOR REFERENCE-GENOME_1L FASTA

%hash_refgenome1L=();
open(IN0,"$genome1");
while(<IN0>)
{
	if(/^>/)
	{
		chomp;
		$head='';
		$head=(split ' ',$_)[0];
		$head=~s/\>//;
	}
	else
	{
		chomp;
		$hash_refgenome1L{$head}=$_;
	}
}
close IN0;

$process_counter='';
$process_counter=$ARGV[2];
chomp($process_counter);

($out_dir1,$out_dir2,$out_dir3,$out_dir4,$out_dir5,$out_dir6)='';
$out_dir1="G".$process_counter."_FASTA";
$out_dir2="G".$process_counter."_CT";
$out_dir3="G".$process_counter."_mfold";
$out_dir4="G".$process_counter."_VARNA";
$out_dir5="G".$process_counter."_info";
$out_dir6="G".$process_counter."_blastN_Sorted";

#	Set your VARNA executable path here
$VARNA_path='';
$VARNA_path="/opt/VARNA/VARNAv3-93.jar";

system("mkdir $out_dir1;mkdir $out_dir2;mkdir $out_dir3;mkdir $out_dir4;mkdir $out_dir5;mkdir $out_dir6;");

print "6 directory generated\n\n";

$out1_fnm='';
$out1_fnm="out1_blastN_on_Genome".$process_counter.".txt";

#	Set your blastN executable path here
$blastN_path='';
$blsatN_path="/opt/blast/blast_2p11p0plus/bin/blastn";


$cmd1='';
$cmd1=$blsatN_path." -query ".$query1." -db ".$genome1." -outfmt \"6 qaccver saccver length qlen slen pident nident mismatch gaps evalue bitscore qstart qend sstart send qseq sseq frames\" -evalue 1000 -max_target_seqs 1000 -word_size 18 -out ".$out1_fnm;

system("$cmd1");		#********** UNHASH-THIS! BEFORE FINAL RUN

print "blsatN completed\n\n";

$out2_fnm='';
$out2_fnm="out1_y1_Genome".$process_counter."_BothStrand_blastN_sort.txt";

$cmd2='';
$cmd2="sort \-t \'\-\' \-Vk2\,2 ".$out1_fnm." \>".$out2_fnm;

system("$cmd2");

$entry_cnt=0;

$new_info='';
$new_info="G".$process_counter."_info.txt";
open(OUT4,">$new_info");

open(IN2,"$out2_fnm");
while(<IN2>)
{
	chomp;
	($G_id,$aln_len,$G_len,$G_st,$G_end,$frame_val)='';

	$G_id=(split '\t',$_)[1];
	$aln_len=(split '\t',$_)[2];
	$G_len=(split '\t',$_)[4];
	$G_st=(split '\t',$_)[13];
	$G_end=(split '\t',$_)[14];
	$frame_val=(split '\t',$_)[17];

	$temp_seq='';
	$temp_seq=$hash_refgenome1L{$G_id};

	if($frame_val=~/\-/)
	{
		$frame_name='';
		$frame_name="neg";
		$entry_cnt++;

		($seq,$rc_seq)='';
		$rc_seq=reverse($temp_seq);
		$rc_seq=~tr/atgcATGC/tacgTACG/;

		($rc_st,$rc_end)='';
		$rc_st=$G_len-$G_st+1;
		$rc_end=$rc_st+$aln_len-1;

		#	Left: 10-x-190
		$frame_mode="1L";

		($frame_st,$frame_end,$frame_range1,$frame_range2,$frame_seq,$frame_cmd1,$frame_color_st,$frame_color_end)='';
	
		if(($rc_st-10)>=1){$frame_st=($rc_st-10);}
		else{$frame_st=1;}

		if(($rc_end+190)<=$G_len){$frame_end=($rc_end+190);}
		else{$frame_end=$G_len;}

		$frame_range1=($frame_st-1);
		$frame_range2=($frame_end-$frame_st+1);

		$frame_seq=substr($rc_seq,$frame_range1,$frame_range2);
		chomp($frame_seq);

		$frame_color_st=($rc_st-$frame_st)+1;
		$frame_color_end=$frame_color_st+$aln_len-1;

		$temp_prefix='';
		$temp_prefix="temp_".$t1;

		($temp_fasta,$temp_ct,$temp_mfold)='';
		$temp_fasta=$temp_prefix.".fasta";
		$temp_ct=$temp_prefix.".ct";
		$temp_mfold=$temp_prefix."_1.png";

		open(OUT3,">$temp_fasta");
		print OUT3 ">s".$entry_cnt."_".$frame_mode."_".$frame_name."".$rc_st."\n".$frame_seq."\n";
		close OUT3;

		$new_fnm_prefix='';
		$new_fnm_prefix="s".$entry_cnt."_".$frame_mode."_".$frame_name."_".$rc_st."-".$rc_end;

		($new_fasta,$new_ct,$new_mfold,$new_VARNA)='';
		$new_fasta=$new_fnm_prefix.".fasta";
		$new_ct=$new_fnm_prefix.".ct";
		$new_mfold=$new_fnm_prefix."_mfold.png";
		$new_VARNA=$new_fnm_prefix."_VARNA.PNG";

		print OUT4 "seq".$entry_cnt."\tFrame:".$frame_mode."\t".$frame_name.":".$G_st."-".$G_end."\tSeqHitPos:".$rc_st."-".$rc_end."\tFramePos:".$frame_st."-".$frame_end."\tFrameColorPos:".$frame_color_st."-".$frame_color_end."\n";

		$frame_cmd1="mfold SEQ=".$temp_fasta." MAX\=1\ >/dev/null 2>/dev/null </dev/null;java \-cp ".$VARNA_path." fr\.orsay\.lri\.varna\.applications\.VARNAcmd -i ".$temp_ct." -o ".$new_VARNA." \-highlightRegion \"".$frame_color_st."-".$frame_color_end."\:fill\=\#FF0000\" \-flat False >/dev/null 2>/dev/null </dev/null; mv ".$temp_fasta." ".$new_fasta.";mv ".$temp_ct." ".$new_ct.";mv ".$temp_mfold." ".$new_mfold.";mv ".$new_fasta." ".$out_dir1.";mv ".$new_ct." ".$out_dir2.";mv ".$new_mfold." ".$out_dir3.";mv ".$new_VARNA." ".$out_dir4.";rm ".$temp_prefix."*";
		system("$frame_cmd1");

		#	MIDDLE: 100-x-100
		$frame_mode="2M";

		($frame_st,$frame_end,$frame_range1,$frame_range2,$frame_seq,$frame_cmd1,$frame_color_st,$frame_color_end)='';
	
		if(($rc_st-100)>=1){$frame_st=($rc_st-100);}
		else{$frame_st=1;}

		if(($rc_end+100)<=$G_len){$frame_end=($rc_end+100);}
		else{$frame_end=$G_len;}

		$frame_range1=($frame_st-1);
		$frame_range2=($frame_end-$frame_st+1);

		$frame_seq=substr($rc_seq,$frame_range1,$frame_range2);
		chomp($frame_seq);

		$frame_color_st=($rc_st-$frame_st)+1;
		$frame_color_end=$frame_color_st+$aln_len-1;

		$temp_prefix='';
		$temp_prefix="temp_".$t1;

		($temp_fasta,$temp_ct,$temp_mfold)='';
		$temp_fasta=$temp_prefix.".fasta";
		$temp_ct=$temp_prefix.".ct";
		$temp_mfold=$temp_prefix."_1.png";

		open(OUT3,">$temp_fasta");
		print OUT3 ">s".$entry_cnt."_".$frame_mode."_".$frame_name."".$rc_st."\n".$frame_seq."\n";
		close OUT3;

		$new_fnm_prefix='';
		$new_fnm_prefix="s".$entry_cnt."_".$frame_mode."_".$frame_name."_".$rc_st."-".$rc_end;

		($new_fasta,$new_ct,$new_mfold,$new_VARNA)='';
		$new_fasta=$new_fnm_prefix.".fasta";
		$new_ct=$new_fnm_prefix.".ct";
		$new_mfold=$new_fnm_prefix."_mfold.png";
		$new_VARNA=$new_fnm_prefix."_VARNA.PNG";

		print OUT4 "seq".$entry_cnt."\tFrame:".$frame_mode."\t".$frame_name.":".$G_st."-".$G_end."\tSeqHitPos:".$rc_st."-".$rc_end."\tFramePos:".$frame_st."-".$frame_end."\tFrameColorPos:".$frame_color_st."-".$frame_color_end."\n";

		$frame_cmd1="mfold SEQ=".$temp_fasta." MAX\=1\ >/dev/null 2>/dev/null </dev/null;java \-cp ".$VARNA_path." fr\.orsay\.lri\.varna\.applications\.VARNAcmd -i ".$temp_ct." -o ".$new_VARNA." \-highlightRegion \"".$frame_color_st."-".$frame_color_end."\:fill\=\#FF0000\" \-flat False >/dev/null 2>/dev/null </dev/null; mv ".$temp_fasta." ".$new_fasta.";mv ".$temp_ct." ".$new_ct.";mv ".$temp_mfold." ".$new_mfold.";mv ".$new_fasta." ".$out_dir1.";mv ".$new_ct." ".$out_dir2.";mv ".$new_mfold." ".$out_dir3.";mv ".$new_VARNA." ".$out_dir4.";rm ".$temp_prefix."*";
		system("$frame_cmd1");

		#	RIGTH: 190-x-10
		$frame_mode="3R";

		($frame_st,$frame_end,$frame_range1,$frame_range2,$frame_seq,$frame_cmd1,$frame_color_st,$frame_color_end)='';
	
		if(($rc_st-190)>=1){$frame_st=($rc_st-190);}
		else{$frame_st=1;}

		if(($rc_end+10)<=$G_len){$frame_end=($rc_end+10);}
		else{$frame_end=$G_len;}

		$frame_range1=($frame_st-1);
		$frame_range2=($frame_end-$frame_st+1);

		$frame_seq=substr($rc_seq,$frame_range1,$frame_range2);
		chomp($frame_seq);

		$frame_color_st=($rc_st-$frame_st)+1;
		$frame_color_end=$frame_color_st+$aln_len-1;

		$temp_prefix='';
		$temp_prefix="temp_".$t1;

		($temp_fasta,$temp_ct,$temp_mfold)='';
		$temp_fasta=$temp_prefix.".fasta";
		$temp_ct=$temp_prefix.".ct";
		$temp_mfold=$temp_prefix."_1.png";

		open(OUT3,">$temp_fasta");
		print OUT3 ">s".$entry_cnt."_".$frame_mode."_".$frame_name."".$rc_st."\n".$frame_seq."\n";
		close OUT3;

		$new_fnm_prefix='';
		$new_fnm_prefix="s".$entry_cnt."_".$frame_mode."_".$frame_name."_".$rc_st."-".$rc_end;

		($new_fasta,$new_ct,$new_mfold,$new_VARNA)='';
		$new_fasta=$new_fnm_prefix.".fasta";
		$new_ct=$new_fnm_prefix.".ct";
		$new_mfold=$new_fnm_prefix."_mfold.png";
		$new_VARNA=$new_fnm_prefix."_VARNA.PNG";

		print OUT4 "seq".$entry_cnt."\tFrame:".$frame_mode."\t".$frame_name.":".$G_st."-".$G_end."\tSeqHitPos:".$rc_st."-".$rc_end."\tFramePos:".$frame_st."-".$frame_end."\tFrameColorPos:".$frame_color_st."-".$frame_color_end."\n";

		$frame_cmd1="mfold SEQ=".$temp_fasta." MAX\=1\ >/dev/null 2>/dev/null </dev/null;java \-cp ".$VARNA_path." fr\.orsay\.lri\.varna\.applications\.VARNAcmd -i ".$temp_ct." -o ".$new_VARNA." \-highlightRegion \"".$frame_color_st."-".$frame_color_end."\:fill\=\#FF0000\" \-flat False >/dev/null 2>/dev/null </dev/null; mv ".$temp_fasta." ".$new_fasta.";mv ".$temp_ct." ".$new_ct.";mv ".$temp_mfold." ".$new_mfold.";mv ".$new_fasta." ".$out_dir1.";mv ".$new_ct." ".$out_dir2.";mv ".$new_mfold." ".$out_dir3.";mv ".$new_VARNA." ".$out_dir4.";rm ".$temp_prefix."*";
		system("$frame_cmd1");

		next;
	}
	if($frame_val!~/\-/)
	{
		$frame_name='';			#	GLOBAL VARIABLE
		$frame_name="pos";		#	GLOBAL VARIABLE
		$entry_cnt++;			#	GLOBAL VARIABLE
	
		#	Left: 10-x-190
		$frame_mode="1L";

		($frame_st,$frame_end,$awk_frame_range,$frame_seq,$frame_cmd1,$frame_color_st,$frame_color_end)='';
	
		if(($G_st-10)>=1){$frame_st=($G_st-10);}
		else{$frame_st=1;}

		if(($G_end+190)<=$G_len){$frame_end=($G_end+190);}
		else{$frame_end=$G_len;}

		$awk_frame_range=($frame_end-$frame_st+1);
		$frame_seq=substr($temp_seq,($frame_st-1),$awk_frame_range);

		$frame_color_st=($G_st-$frame_st)+1;
		$frame_color_end=$frame_color_st+$aln_len-1;

		$temp_prefix='';
		$temp_prefix="temp_".$t1;

		($temp_fasta,$temp_ct,$temp_mfold)='';
		$temp_fasta=$temp_prefix.".fasta";
		$temp_ct=$temp_prefix.".ct";
		$temp_mfold=$temp_prefix."_1.png";

		open(OUT3,">$temp_fasta");
		print OUT3 ">s".$entry_cnt."_".$frame_mode."_".$frame_name."".$G_st."\n".$frame_seq."\n";
		close OUT3;

		$new_fnm_prefix='';
		$new_fnm_prefix="s".$entry_cnt."_".$frame_mode."_".$frame_name."_".$G_st."-".$G_end;

		($new_fasta,$new_ct,$new_mfold,$new_VARNA)='';
		$new_fasta=$new_fnm_prefix.".fasta";
		$new_ct=$new_fnm_prefix.".ct";
		$new_mfold=$new_fnm_prefix."_mfold.png";
		$new_VARNA=$new_fnm_prefix."_VARNA.PNG";

		print OUT4 "seq".$entry_cnt."\tFrame:".$frame_mode."\t".$frame_name."\tSeqHitPos:".$G_st."-".$G_end."\tFramePos:".$frame_st."-".$frame_end."\tFrameColorPos:".$frame_color_st."-".$frame_color_end."\n";

		$frame_cmd1="mfold SEQ=".$temp_fasta." MAX\=1\ >/dev/null 2>/dev/null </dev/null;java \-cp ".$VARNA_path." fr\.orsay\.lri\.varna\.applications\.VARNAcmd -i ".$temp_ct." -o ".$new_VARNA." \-highlightRegion \"".$frame_color_st."-".$frame_color_end."\:fill\=\#FF0000\" \-flat False >/dev/null 2>/dev/null </dev/null; mv ".$temp_fasta." ".$new_fasta.";mv ".$temp_ct." ".$new_ct.";mv ".$temp_mfold." ".$new_mfold.";mv ".$new_fasta." ".$out_dir1.";mv ".$new_ct." ".$out_dir2.";mv ".$new_mfold." ".$out_dir3.";mv ".$new_VARNA." ".$out_dir4.";rm ".$temp_prefix."*";
		system("$frame_cmd1");

		#	MIDDLE: 100-x-100
		$frame_mode="2M";

		($frame_st,$frame_end,$awk_frame_range,$frame_seq,$frame_cmd1,$frame_color_st,$frame_color_end)='';
	
		if(($G_st-100)>=1){$frame_st=($G_st-100);}
		else{$frame_st=1;}

		if(($G_end+100)<=$G_len){$frame_end=($G_end+100);}
		else{$frame_end=$G_len;}

		$awk_frame_range=($frame_end-$frame_st+1);
		$frame_seq=substr($temp_seq,($frame_st-1),$awk_frame_range);

		$frame_color_st=($G_st-$frame_st)+1;
		$frame_color_end=$frame_color_st+$aln_len-1;

		$temp_prefix='';
		$temp_prefix="temp_".$t1;

		($temp_fasta,$temp_ct,$temp_mfold)='';
		$temp_fasta=$temp_prefix.".fasta";
		$temp_ct=$temp_prefix.".ct";
		$temp_mfold=$temp_prefix."_1.png";

		open(OUT3,">$temp_fasta");
		print OUT3 ">s".$entry_cnt."_".$frame_mode."_".$frame_name."".$G_st."\n".$frame_seq."\n";
		close OUT3;

		$new_fnm_prefix='';
		$new_fnm_prefix="s".$entry_cnt."_".$frame_mode."_".$frame_name."_".$G_st."-".$G_end;

		($new_fasta,$new_ct,$new_mfold,$new_VARNA)='';
		$new_fasta=$new_fnm_prefix.".fasta";
		$new_ct=$new_fnm_prefix.".ct";
		$new_mfold=$new_fnm_prefix."_mfold.png";
		$new_VARNA=$new_fnm_prefix."_VARNA.PNG";

		print OUT4 "seq".$entry_cnt."\tFrame:".$frame_mode."\t".$frame_name."\tSeqHitPos:".$G_st."-".$G_end."\tFramePos:".$frame_st."-".$frame_end."\tFrameColorPos:".$frame_color_st."-".$frame_color_end."\n";

		$frame_cmd1="mfold SEQ=".$temp_fasta." MAX\=1\ >/dev/null 2>/dev/null </dev/null;java \-cp ".$VARNA_path." fr\.orsay\.lri\.varna\.applications\.VARNAcmd -i ".$temp_ct." -o ".$new_VARNA." \-highlightRegion \"".$frame_color_st."-".$frame_color_end."\:fill\=\#FF0000\" \-flat False >/dev/null 2>/dev/null </dev/null; mv ".$temp_fasta." ".$new_fasta.";mv ".$temp_ct." ".$new_ct.";mv ".$temp_mfold." ".$new_mfold.";mv ".$new_fasta." ".$out_dir1.";mv ".$new_ct." ".$out_dir2.";mv ".$new_mfold." ".$out_dir3.";mv ".$new_VARNA." ".$out_dir4.";rm ".$temp_prefix."*";
		system("$frame_cmd1");

		#	Right: 190-x-10
		$frame_mode="3R";

		($frame_st,$frame_end,$awk_frame_range,$frame_seq,$frame_cmd1,$frame_color_st,$frame_color_end)='';
	
		if(($G_st-190)>=1){$frame_st=($G_st-190);}
		else{$frame_st=1;}

		if(($G_end+10)<=$G_len){$frame_end=($G_end+10);}
		else{$frame_end=$G_len;}

		$awk_frame_range=($frame_end-$frame_st+1);
		$frame_seq=substr($temp_seq,($frame_st-1),$awk_frame_range);

		$frame_color_st=($G_st-$frame_st)+1;
		$frame_color_end=$frame_color_st+$aln_len-1;

		$temp_prefix='';
		$temp_prefix="temp_".$t1;

		($temp_fasta,$temp_ct,$temp_mfold)='';
		$temp_fasta=$temp_prefix.".fasta";
		$temp_ct=$temp_prefix.".ct";
		$temp_mfold=$temp_prefix."_1.png";

		open(OUT3,">$temp_fasta");
		print OUT3 ">s".$entry_cnt."_".$frame_mode."_".$frame_name."".$G_st."\n".$frame_seq."\n";
		close OUT3;

		$new_fnm_prefix='';
		$new_fnm_prefix="s".$entry_cnt."_".$frame_mode."_".$frame_name."_".$G_st."-".$G_end;

		($new_fasta,$new_ct,$new_mfold,$new_VARNA)='';
		$new_fasta=$new_fnm_prefix.".fasta";
		$new_ct=$new_fnm_prefix.".ct";
		$new_mfold=$new_fnm_prefix."_mfold.png";
		$new_VARNA=$new_fnm_prefix."_VARNA.PNG";

		print OUT4 "seq".$entry_cnt."\tFrame:".$frame_mode."\t".$frame_name."\tSeqHitPos:".$G_st."-".$G_end."\tFramePos:".$frame_st."-".$frame_end."\tFrameColorPos:".$frame_color_st."-".$frame_color_end."\n";

		$frame_cmd1="mfold SEQ=".$temp_fasta." MAX\=1\ >/dev/null 2>/dev/null </dev/null;java \-cp ".$VARNA_path." fr\.orsay\.lri\.varna\.applications\.VARNAcmd -i ".$temp_ct." -o ".$new_VARNA." \-highlightRegion \"".$frame_color_st."-".$frame_color_end."\:fill\=\#FF0000\" \-flat False >/dev/null 2>/dev/null </dev/null; mv ".$temp_fasta." ".$new_fasta.";mv ".$temp_ct." ".$new_ct.";mv ".$temp_mfold." ".$new_mfold.";mv ".$new_fasta." ".$out_dir1.";mv ".$new_ct." ".$out_dir2.";mv ".$new_mfold." ".$out_dir3.";mv ".$new_VARNA." ".$out_dir4.";rm ".$temp_prefix."*";
		system("$frame_cmd1");

	}
}
close IN2;
close OUT4;

$cmd0='';
$cmd0="mv ".$out1_fnm." ".$out2_fnm." ".$out_dir6;

system("$cmd0");

$cmd3='';
$cmd3="mv ".$new_info." ".$out_dir5;
system("$cmd3");
