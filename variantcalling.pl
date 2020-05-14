#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Find;
use File::Path;
#use Thread::Semaphore;

use threads;
use FindBin;


#my $steptime="5,3,6,20"; #9,11,3,8,13,32 "13,2,4,11,15"
#cut,bwa,picard,vqsr,gvcf;
my $steptime="5,8,3,6,20";

#c(5,2,0.5,1,5)
#my @steptime=(3, 5.5, 1.5, 1.5, 4.5, 8); ##multiple thread 4
##run in local. thread4: 2,0.5,1,5,11
##run in cluster.thread4: 5,1,2,7,20
#$| = 1;

my $lib = "$FindBin::RealBin/lib"; ## ##bwa,picard,gatk..

my $user = `echo \$USER`;
my $pwd = `pwd`;
my $jobname_pre = $$;

############
##Parameter
############
my $opt_help;
my $get_config = ""; #"$FindBin::RealBin/calling.config.hg38";
my $input = "";
my $opt_out_dir = "";
my $where = ""; #local or qsub or sbatch
my $steps = "1,2,3,4";
my $email = "";
my $nthreads = 10;
my $nthreads_perjob = 2;
my $vcf = "";

my $ploidy = 2;
##
#my $fastq_size = 3.0;
my $gvcf_intervals = "";
my $q64 = "false";
my $sequencing_stype = "exome";

my $dup_keep = "";
##
my $cutadapt;
my $fastq2gvcf;
my $genotypeGVCFs;
my $vqsr;
my $qplot_run;
my $annovar_run;
my $beagle4_run;
my $combinegvcf_run;

my $buildver = "hg38";
##
my $sleep = 30 ;

GetOptions("help|h"  => \$opt_help,
	"get_config|get_config=s" => \$get_config,
	"input|input=s" => \$input,
	"out_dir|output_dir=s" => \$opt_out_dir,
	"where|where=s" => \$where,
	"steps|steps=s" => \$steps,
	"email|email=s" => \$email,
	"nthreads|nthreads=f" => \$nthreads,
	"nthreads_perjob|nthreads_perjob=f" => \$nthreads_perjob,
	"vcf|vcf=s" => \$vcf,

	"ploidy|ploidy=f" => \$ploidy,
##
	#"size|size=f" => \$fastq_size,
	"gvcf_intervals|gvcf_intervals=s" => \$gvcf_intervals,
	"Q64|Q64=s" => \$q64,
	"sequencing_stype|sequencing_stype=s" => \$sequencing_stype,
	"dup_keep|dup_keep=s" => \$dup_keep,

##
	"cutadapt|cutadapt" => \$cutadapt,
	"fastq2gvcf|fastq2gvcf" => \$fastq2gvcf,
	"genotypeGVCFs|genotypeGVCFs" => \$genotypeGVCFs,
	"vqsr|vqsr" => \$vqsr,
	"qplot|qplot" => \$qplot_run,
	"annovar|annovar" => \$annovar_run,
	"beagle4|beagle4" => \$beagle4_run,
	"combinegvcf|combinegvcf" => \$combinegvcf_run,
	"buildver|buildver=s" => \$buildver,
##
	"sleep|sleep=f" => \$sleep,
	"jobname_pre|jobname_pre=f" => \$jobname_pre,
	"steptime|steptime=s" => \$steptime,
	);

#########################################
## how to use
#######################################
&usage() if $opt_help;

die "-buildver is not provided! please use -buildver hg19 or hg38 \n" if !$buildver;
#die "-get_config is not provided!\n" if !$get_config;
die "-cutadapt, -fastq2gvcf, -combinegvcf -genotypeGVCFs -vqsr -qplot -annovar or -beagle4  is not provided! OR -h to help\n" if (!$cutadapt && !$fastq2gvcf && !$genotypeGVCFs && !$vqsr && !$qplot_run && !$annovar_run && !$beagle4_run && !$combinegvcf_run);
die "-input is not provided!\n" if ($cutadapt || $fastq2gvcf || $genotypeGVCFs || $qplot_run || $combinegvcf_run) && !$input;

die "-output_dir is not provided!\n" if !$opt_out_dir;
die "-where is not provided. -local, -qsub or -sbatch!\n" if ($where ne "local" && $where ne "qsub" && $where ne "sbatch");
die "-email is not provided when running in $where !\n" if ($where eq "qsub" || $where eq "sbatch") && !$email;
die "-nthreads is not provided when running in local !\n" if $where eq "local" && !$nthreads;
die "-vcf files is not provided when running in local !\n" if (!$genotypeGVCFs) && ($vqsr || $annovar_run) && !$input;
die "-sequencing_stype is not provided when running in qplot or vqsr. exome or whole !\n" if ($qplot_run || $vqsr) && !$sequencing_stype;

#########################################
## config files and check all softwares and databases
#######################################

if($get_config)
{}
elsif ($buildver eq "hg38" || $buildver eq "hg19" || $buildver) {
	$get_config = "$FindBin::RealBin/calling.config.$buildver";
}else
{ die "Error -buildver hg38 or hg19  "; }


####
## show options
#####

&options();

my %hash_config = &get_config($get_config);

############################
# sleep time in different where.
######################
if($sleep > 30)
{
	$sleep = 30;
}

#####################################################
##Make dir
#####################################################

chop($opt_out_dir) if $opt_out_dir =~ /\/$/;

eval{mkpath($opt_out_dir,0,0755)};
if($@)
{
   warn("Make path [$opt_out_dir] failed:\n$@");
}

#my $data_dir = "$opt_out_dir/rawdata"; #fastq files

my $log_dir = "$opt_out_dir/logs"; #log files

mkdir ($log_dir,0755) unless (-d $log_dir);


######################

my $process = "$opt_out_dir/$jobname_pre.process.log";
my $error = "$opt_out_dir/$jobname_pre.error.log";
my $commands = "$opt_out_dir/$jobname_pre.commands.log";
###############
##
############
my $startflag=0;

if ($cutadapt)
{
	$input = &cutadapt($nthreads, $opt_out_dir, \%hash_config, $input); 
}

################
##run fastq to gvcf:steps 7
#################

if ($fastq2gvcf)
{
	$input = &fastq2gvcf($nthreads, $opt_out_dir, \%hash_config, $input, $steps); ##return g.vcf.gz files
}

################
##qplot
#################
if ($qplot_run && !$fastq2gvcf)
{
	&qplot($nthreads, $opt_out_dir, \%hash_config, $input); ##the sample_id only from the index file

}

################
##merge gvcf to a big gvcf
#################
if ($combinegvcf_run)
{
	$input = &combineGVCFs($nthreads, $opt_out_dir, \%hash_config, $input); ##return one big g.vcf.gz
}

################
##gvcf to vcf
#################
if ($genotypeGVCFs)
{
	$input = &genotypeGVCFs($nthreads, $opt_out_dir, \%hash_config, $input);
}

################
##vqsr
#################
if ($vqsr)
{
	$input = &vqsr($nthreads, $opt_out_dir, \%hash_config, $input);
}

################
##annovar
#################
if ($annovar_run)
{
	$input = &annovar($nthreads, $opt_out_dir, \%hash_config, $input);
}

####################################################
#### run the progress
####################################################
sub cutadapt()
{
	
	my ($nthreads, $opt_out_dir, $hash_config, $input) = @_;
=a
	if ($nthreads>24 && !$nthreads) 
	{
		$nthreads = 24;
	}

	my $nthreads_perjob = 1;
	my $job_run = int ($nthreads/$nthreads_perjob);

	###mkdir dir
	my $bam_dir = "$opt_out_dir/bamfiles"; #bam files

	mkdir ($bam_dir,0755) unless (-d $bam_dir);

	#####
	open(TO,">>$commands") or die();
	print TO "############################\n";
	print TO "#pwd is \t$pwd###\n";

	print "\n###################\n cutadapt starts running; -->\n";
	print "\n***********check config files\n";
	print "***********\n";

	############################################################################################
	####steps control
	&checkfile($hash_config{"cutadapt"},$hash_config{"ssub_array"}); #

	print "***********\n";

	print STDERR "\tcutadapt starts running; ".&nowtime." -->";
	print TO "###cutadapt starts running; ".&nowtime." -->\n\n";

	my @sample_id = &get_index($input);

	#########
	my @cutadapt_commandline;

	open(TO2,">$opt_out_dir/cut.$jobname_pre.fastq.list") or die();
	$input = "$opt_out_dir/cut.$jobname_pre.fastq.list";
	
	for (my $i = 0; $i < @sample_id; $i++) 
	{
		my @original_file = split(/\t/,$sample_id[$i]);

		if (exists $original_file[2]) #### SM in sam files
		{

			$cutadapt_commandline[$i] = qq/$hash_config{"cutadapt"} -f fastq -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAA -e 0.1 -m 50 -n 2 --trim-n -q 15,10 -o $bam_dir\/$original_file[0].R1.tmp.fq.gz -p $bam_dir\/$original_file[0].R2.tmp.fq.gz $original_file[1] $original_file[2] 1>$log_dir\/$original_file[0].log/;


			print TO2 join("\t",$original_file[0],"$bam_dir/$original_file[0].R1.fq.gz","$bam_dir/$original_file[0].R2.fq.gz")."\n";

		}
		else
		{
			$cutadapt_commandline[$i] = qq/$hash_config{"cutadapt"} -f fastq -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -e 0.1 -q 15,10 -m 30 -n 2 --trim-n -o $bam_dir\/$original_file[0].R1.fq.gz $original_file[1]  1>$log_dir\/$original_file[0].log/;

			print TO2 join("\t",$original_file[0],"$bam_dir/$original_file[0].R1.fq.gz")."\n";
		}
	}
	close TO2;

	if ($where eq "local") 
	{
		my $whichSubroutine = \&stand_command;

		&multifork($job_run,$whichSubroutine,"cutadapt-$jobname_pre",@cutadapt_commandline);
	}
	elsif($where eq "sbatch")
	{
		my @bsub_options;

		my $time = 20;

		@bsub_options=("4000",$time, "cutadapt-$jobname_pre", $email, $nthreads_perjob);
    	&cluster_bsub(\@cutadapt_commandline,\@bsub_options);
		
	}

	$startflag = 1;

	print STDERR &nowtime." finished\n";
	print TO "\n####".&nowtime." finished\n\n\n";
	close TO;
=cut
	return $input;
}




sub fastq2gvcf()
{
	
	my ($nthreads, $opt_out_dir, $hash_config, $input, $steps) = @_;

	if ($nthreads>50 && !$nthreads) 
	{
		$nthreads = 50;
	}

	my $nthreads_perjob = 2;
	my $job_run = int ($nthreads/$nthreads_perjob);

	###mkdir dir
	my $bam_dir = "$opt_out_dir/bamfiles"; #bam files
	my $met_dir = "$opt_out_dir/metrics"; #metrics files
	my $other_dir = "$opt_out_dir/otherfiles"; #other files
	my $gvcf_dir = "$opt_out_dir/gvcf"; #gvcffiles
	
	mkdir ($bam_dir,0755) unless (-d $bam_dir);
	mkdir ($met_dir,0755) unless (-d $met_dir);
	mkdir ($other_dir,0755) unless (-d $other_dir);

	#####
	open(TO,">>$commands") or die();
	print TO "############################\n";
	print TO "#pwd is \t$pwd###\n";

	print "\n###################\n fastq2gvcf starts running; -->\n";
	print "\n***********check config files\n";
	print "***********\n";

	############################################################################################
	####steps control
	&checkfile($hash_config{"opt_ref_genome"},$hash_config{"Mills_and_1000G_gold_standard_indels"},$hash_config{"dbsnp"},$hash_config{"samtools"},$hash_config{"bwa"},$hash_config{"gatk"},$hash_config{"ssub_array"}); #$hash_config{"KG_phase1_indels"},

	print "***********\n";

############
	my @steps_control = split(",",$steps); #$steps = "3,4";
	my @stepsflag;


	for (my $var = 0; $var < 5; $var++) 
	{
		$stepsflag[$var] = 0;
	}

	for (my $var = 0; $var < @steps_control; $var++) 
	{
		$stepsflag[$steps_control[$var]] = 1;
	}


##########
## time
	my $timeMax=0;
	my @sample_id = &get_index_check($input);

	for (my $var = 0, my $count = @sample_id; $var < $count; $var++) 
	{
		my @original_file = split(/\t/,$sample_id[$var]);

		my $size = -s "$original_file[1]";
		
		#print "$original_file[$var][0]"."\t".$size."\n";

		$size = $size/3000000000;

		if ($timeMax < $size) {
			$timeMax = $size;
		}
	}

	$timeMax += 1;

	my @steptime = split(",",$steptime);


#############################
##
################################
	my @bwa_commandline;
	my @picard_commandline;
	my @bqsr_commandline;
	my @gvcf_commandline;

	#########

	if ($stepsflag[1]) 
	{
		print STDERR "\tstep1\tBWA mem starts running; ".&nowtime." -->";
		print TO "###step1\tBWA mem starts running; ".&nowtime." -->\n\n";
		
	my @sample_id = &get_index($input);

	open(TO2,">$opt_out_dir/bwa.$jobname_pre.list") or die();

	$input = "$opt_out_dir/bwa.$jobname_pre.list";

	for (my $i = 0; $i < @sample_id; $i++) 
	{
		my @original_file = split(/\t/,$sample_id[$i]);
		my $sample = $original_file[0];
		my $time = &nowtime;

		my $r = "\'\@RG\\tID:$sample\\tLB:$sample\\tPL:ILLUMINA\\tSM:$sample\\tDT:$time\'";				

		if (exists $original_file[2]) #### SM in sam files
		{
			$bwa_commandline[$i] = qq/$hash_config{"bwa"} mem -M -t $nthreads_perjob -R $r $hash_config{"opt_ref_genome"} $original_file[1] $original_file[2] 2>$log_dir\/$sample.log | $hash_config{"samtools"} view -bS -\@ $nthreads_perjob -o $bam_dir\/$sample.tmp.bam - 2>>$log_dir\/$sample.log /;					
		}
		else
		{
			$bwa_commandline[$i] = qq/$hash_config{"bwa"} mem -M -t $nthreads_perjob -R $r $hash_config{"opt_ref_genome"} $original_file[1] 2>$log_dir\/$sample.log | $hash_config{"samtools"} view -bS -\@ $nthreads_perjob -o $bam_dir\/$sample.tmp.bam - 2>>$log_dir\/$sample.log /;			
		}

		$bwa_commandline[$i] .= " && $hash_config{'samtools'} sort -\@ $nthreads_perjob -o $bam_dir\/$sample.sorted.bam $bam_dir\/$sample.tmp.bam 2>>$log_dir/$sample.log && rm -rf $bam_dir\/$sample.tmp.bam ";

		$bwa_commandline[$i] .= " && $hash_config{'samtools'} index $bam_dir\/$sample.sorted.bam "; #index

		print TO2 join("\t",$sample.".sorted","$bam_dir\/$sample.sorted.bam")."\n";

	}
	close TO2;


	if ($where eq "local") 
	{
		my $whichSubroutine = \&stand_command;

		&multifork($job_run,$whichSubroutine,"bwa-$jobname_pre",@bwa_commandline);
	}
	elsif($where eq "sbatch")
	{
		my @bsub_options;

		my $time = int ($timeMax*$steptime[1]+5);

		@bsub_options=("8000",$time, "bwa-$jobname_pre", $email, $nthreads_perjob);
    	&cluster_bsub (\@bwa_commandline,\@bsub_options);
		
	}

	$startflag = 1;

	print STDERR &nowtime." finished\n";
	print TO "\n####".&nowtime." finished\n\n\n";

	}


###########
	if ($stepsflag[2]) 
	{
	print STDERR "\tstep2\tPicard mark MarkDuplicates starts running; ".&nowtime." -->";
	print TO "#####step2\tPicard mark MarkDuplicates starts running; ".&nowtime." -->\n\n";

	$nthreads_perjob = 1;
	$job_run = int ($nthreads/$nthreads_perjob);

	my @sample_id = &get_index($input);

	open(TO2,">$opt_out_dir/picard.$jobname_pre.list") or die();

	$input = "$opt_out_dir/picard.$jobname_pre.list";

	for (my $i = 0; $i < @sample_id; $i++) 
	{
		my @original_file = split(/\t/,$sample_id[$i]);
		my $sample = $original_file[0];

		$picard_commandline[$i] = qq/$hash_config{"gatk"} --java-options "-Xmx16g -XX:+UseParallelGC -XX:ParallelGCThreads=1" MarkDuplicates -I $original_file[1] -O $bam_dir\/$sample.dedupped.bam -M $met_dir\/$sample.dedupped.metrics --CREATE_INDEX true -VALIDATION_STRINGENCY LENIENT  2>>$log_dir\/$sample.log 1>>$log_dir\/$sample.log/;

		if ($startflag) {
			$picard_commandline[$i] .= qq/ && rm -f $original_file[1] $original_file[1].bai /;
		}

		print TO2 join("\t",$sample.".dedupped","$bam_dir\/$sample.dedupped.bam")."\n";

	}
	close TO2;

	if ($where eq "local") 
	{
		my $whichSubroutine = \&stand_command;

    	&multifork($job_run,$whichSubroutine,"picard-$jobname_pre",@picard_commandline);
	}
	elsif($where eq "sbatch")
	{
		my @bsub_options;

		my $time = int ($timeMax*$steptime[2]+5);

		@bsub_options=("16000",$time, "picard-$jobname_pre", $email, $nthreads_perjob);
    	&cluster_bsub (\@picard_commandline,\@bsub_options);	
	}

	$startflag = 1;

	print STDERR &nowtime." finished\n";
	print TO "\n####".&nowtime." finished\n\n\n";
	}


###########
	if ($stepsflag[3]) 
	{

	print STDERR "\tstep3\tGATK bqsr starts running; ".&nowtime." -->";
	print TO "#########step3\tGATK bqsr starts running; ".&nowtime." -->\n\n";

	$nthreads_perjob = 1;
	$job_run = int ($nthreads/$nthreads_perjob);

	my @sample_id = &get_index($input);

	open(TO2,">$opt_out_dir/bqsr.$jobname_pre.list") or die();

	$input = "$opt_out_dir/bqsr.$jobname_pre.list";

	for (my $i = 0; $i < @sample_id; $i++) 
	{
		my @original_file = split(/\t/,$sample_id[$i]);
		my $sample = $original_file[0];

		$bqsr_commandline[$i] = qq/ $hash_config{"gatk"} --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=1" BaseRecalibrator -R $hash_config{"opt_ref_genome"} -I $original_file[1] --known-sites $hash_config{"dbsnp"} --known-sites $hash_config{"Mills_and_1000G_gold_standard_indels"}  -O $other_dir\/$sample.recal.grp --default-base-qualities 1 2>>$log_dir\/$sample.log 1>>$log_dir\/$sample.log/;  ##BaseRecalibratorSpark  --spark-master local[$nthreads_perjob]  #--sparkRunner LOCAL -RF BadCigar  ##--known-sites $hash_config{"KG_phase1_indels"}


		$bqsr_commandline[$i] .= qq/ && $hash_config{"gatk"} --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=1" ApplyBQSR -R $hash_config{"opt_ref_genome"} -I $original_file[1] -bqsr $other_dir\/$sample.recal.grp -O $bam_dir\/$sample.recal.bam 2>>$log_dir\/$sample.log 1>>$log_dir\/$sample.log /;  ## ApplyBQSRSpark --sparkRunner SPARK --spark-master local[$nthreads_perjob]


		if ($startflag) {
			$bqsr_commandline[$i] .= qq/ && rm -f $original_file[1] /;

			$original_file[1] =~ s/dedupped.bam/dedupped.bai/;

			$bqsr_commandline[$i] .= qq/ $original_file[1] /;

		}

		my @samples2 = split("\\.",$sample);

		print TO2 join("\t",$samples2[0],"$bam_dir\/$sample.recal.bam")."\n";

	}
	close TO2;

	if ($where eq "local") 
	{
		my $whichSubroutine = \&stand_command;
    	&multifork($job_run,$whichSubroutine,"bqsr-$jobname_pre",@bqsr_commandline);
	}
	elsif($where eq "sbatch")
	{
		my @bsub_options;

		my $time = int ($timeMax*$steptime[3]+5);

		@bsub_options=("8000",$time, "bqsr-$jobname_pre", $email, $nthreads_perjob);
    	&cluster_bsub (\@bqsr_commandline,\@bsub_options);
	
	}

	$startflag = 1;

	print STDERR &nowtime." finished\n";
	print TO "\n####".&nowtime." finished\n\n\n";
	}

#######################
###########

my $gvcfInput;

	if ($stepsflag[4]) 
	{

	print STDERR "\tstep4\tGATK to gvcf starts running; ".&nowtime." -->";
	print TO "######step4\tGATK to gvcf starts running; ".&nowtime." -->\n\n";

	$nthreads_perjob = 1;
	$job_run = int ($nthreads/$nthreads_perjob);

	mkdir ($gvcf_dir,0755) unless (-d $gvcf_dir);

	my @sample_id = &get_index($input);

	open(TO2,">$opt_out_dir/gvcf.$jobname_pre.list") or die();

	$gvcfInput = "$opt_out_dir/gvcf.$jobname_pre.list";

	for (my $i = 0; $i < @sample_id; $i++) 
	{
		my @original_file = split(/\t/,$sample_id[$i]);
		my $sample = $original_file[0];

		$gvcf_commandline[$i] = qq/ $hash_config{"gatk"} --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=1" HaplotypeCaller -R $hash_config{"opt_ref_genome"} -I $original_file[1] --dbsnp $hash_config{"dbsnp"} -O $gvcf_dir\/$sample.g.vcf.gz --emit-ref-confidence GVCF /; # HaplotypeCallerSpark  --spark-master local[$nthreads_perjob] --tmp-dir  --native-pair-hmm-threads

		if ($gvcf_intervals) 
		{
			$gvcf_commandline[$i] .= " -intervals $gvcf_intervals ";
		}
		if ($ploidy != 2) 
		{
			$gvcf_commandline[$i] .= " -ploidy $ploidy ";
		}

		$gvcf_commandline[$i] .= " 2>>$log_dir\/$sample.log 1>>$log_dir\/$sample.log ";

		print TO2 join("\t",$sample,"$gvcf_dir\/$sample.g.vcf.gz")."\n";

	}
	close TO2;

	if ($where eq "local") 
	{
		my $whichSubroutine = \&stand_command;

    	&multifork($job_run,$whichSubroutine,"gvcf-$jobname_pre",@gvcf_commandline);
	}
	elsif($where eq "sbatch")
	{

		my @bsub_options;

		my $time = int ($timeMax*$steptime[4]+15);

		@bsub_options=("8000",$time, "gvcf-$jobname_pre", $email, $nthreads_perjob);
    	&cluster_bsub (\@gvcf_commandline,\@bsub_options);	
	}

	$startflag = 1;

	print STDERR &nowtime." finished\n";
	print TO "\n####".&nowtime." finished\n\n\n";

	}



#####################		
=a
	if ($where eq "local") 
	{
		my $whichSubroutine = \&stand_command;

		if(@bwa_commandline)
		{
    		print "1\n";
    		&multifork($job_run,$whichSubroutine,@bwa_commandline);

    		print "2\n";
		}
		if(@picard_commandline)
		{
    		print "3\n";
    		&multifork($job_run,$whichSubroutine,@picard_commandline);
		}
		if(@bqsr_commandline)
		{
    		&multifork($job_run,$whichSubroutine,@bqsr_commandline);
		}
		if(@gvcf_commandline)
		{
    		&multifork($job_run,$whichSubroutine,@gvcf_commandline);
		}

=cut

#####################################################################################################
	#############

	close TO;

	if ($qplot_run)
	{
		&qplot($nthreads, $opt_out_dir, \%$hash_config, $input);
	}

	return $gvcfInput;

}


sub qplot()
{
	my ($nthreads, $opt_out_dir, $hash_config, $input) = @_;

	if ($nthreads>24 && !$nthreads) 
	{
		$nthreads = 24;
	}

	my $nthreads_perjob = 1;
	my $job_run = int ($nthreads/$nthreads_perjob);

	##mkdir dir
	my $qplot_dir = "$opt_out_dir/qplot"; #qplotfiles
	mkdir ($qplot_dir,0755) unless (-d $qplot_dir);
	###########

	#####
	open(TO,">>$commands") or die();
	print TO "############################\n";
	print TO "\n###################\nqplot starts running; ".&nowtime." -->\n\n\n";

	print "\n###################\nqplot starts running; ".&nowtime." -->";
	print "\n***********check config files\n";

	&checkfile($hash_config{"opt_ref_genome"},$hash_config{"ssub_array"},$hash_config{"qplot_tbl"},$hash_config{"qplot_winsize100"},$hash_config{"AgilentSureSelect50mb"},$hash_config{"qplot"}); 
	print "***********\n";

	############
	my @sample_id = &get_index($input);

	my @bamLabels;
	my @inputbam;
	
	for (my $i = 0; $i < @sample_id; $i++) 
	{
		my @original_file = split(/\t/,$sample_id[$i]);

		$bamLabels[$i] = $original_file[0];
		$inputbam[$i] = $original_file[1];
	}

	my @qplot_commandline;
	my $qplot_num = 6;

	for (my $i = 0, my $j=0; $i < @sample_id; $i=$i+$qplot_num, $j++) 
	{
		my $inputbam;
		my $bamLabels;
		
		$qplot_commandline[$j] = qq/$hash_config{"qplot"} --reference $hash_config{"opt_ref_genome"} --dbsnp $hash_config{"qplot_tbl"}  /;

		if ($sequencing_stype eq "whole") 
		{
		}
		elsif($sequencing_stype eq "exome")
		{
			$qplot_commandline[$j] .= qq/ --regions $hash_config{"AgilentSureSelect50mb"} /; #--gccontent $hash_config{"qplot_winsize100"}
		}
		elsif($sequencing_stype eq "cpg")
		{
			$qplot_commandline[$j] .= qq/ --regions $hash_config{"cpgIslandExt"} /; #--gccontent $hash_config{"qplot_winsize100"}
		}
		else
		{
			$qplot_commandline[$j] .= qq/ --regions $sequencing_stype /; #--gccontent 				
		}


		if ($i+$qplot_num-1<@sample_id-1) 
		{
			$inputbam = join(" ",@inputbam[$i..$i+$qplot_num-1]);
			$bamLabels = join(",",@bamLabels[$i..$i+$qplot_num-1]);

			$qplot_commandline[$j] .= qq/ --plot $qplot_dir\/$bamLabels[$i]-$bamLabels[$i+$qplot_num-1].pdf --stats $qplot_dir\/$bamLabels[$i]-$bamLabels[$i+$qplot_num-1].stats --Rcode $qplot_dir\/$bamLabels[$i]-$bamLabels[$i+$qplot_num-1].R --bamLabel $bamLabels /;
		}
		else
		{
			$inputbam = join(" ",@inputbam[$i..@inputbam-1]);
			$bamLabels = join(",",@bamLabels[$i..@bamLabels-1]);

			$qplot_commandline[$j] .= qq/ --plot $qplot_dir\/$bamLabels[$i]-$bamLabels[@bamLabels-1].pdf --stats $qplot_dir\/$bamLabels[$i]-$bamLabels[@bamLabels-1].stats --Rcode $qplot_dir\/$bamLabels[$i]-$bamLabels[@bamLabels-1].R --bamLabel $bamLabels /;		
		}

		$qplot_commandline[$j] .= qq/ $inputbam 1>>$log_dir\/qplot.log 2>>$log_dir\/qplot.log /; 

	}

	if ($where eq "local") 
	{
		my $whichSubroutine = \&stand_command;

		if(@qplot_commandline)
		{
    		&multifork($job_run,$whichSubroutine,"qplot-$jobname_pre",@qplot_commandline);
		}

	}
	elsif($where eq "sbatch")
	{
		my @bsub_options;

		if(@qplot_commandline)
		{
			@bsub_options=("8000",50, "qplot-$jobname_pre", $email, $nthreads_perjob);
    		&cluster_bsub (\@qplot_commandline,\@bsub_options);
		}
		
	}
	
	print STDERR "\tqplot_commandline is finished\n";
	print TO "\n####".&nowtime." finished\n\n\n";
	close TO;
}


sub combineGVCFs()
{
	my ($nthreads, $opt_out_dir, $hash_config, $input) = @_;

	$nthreads_perjob = 1;

	if ($nthreads>24 && !$nthreads) 
	{
		$nthreads = 24;
	}

	my $job_run = int ($nthreads/$nthreads_perjob);

	##mkdir dir
	my $gvcf_dir = "$opt_out_dir/gvcf"; #gvcffiles
	
	mkdir ($gvcf_dir,0755) unless (-d $gvcf_dir);
	####
	open(TO,">>$commands") or die();
	print TO "############################\n";
	print TO "\n###################\nCombineGVCFs starts running; ".&nowtime." -->\n\n\n";

	print "\n###################\nCombineGVCFs starts running; -->\n";
	print "\n***********check config files\n";

	&checkfile($hash_config{"opt_ref_genome"},$hash_config{"gatk"},$hash_config{"ssub_array"}); #
	
	print "***********\n";
#################
##
######################
	my @sample_id = &get_index($input);
	my $inputgvcf="";

	for (my $i = 0; $i < @sample_id; $i++) 
	{	
		my @original_file = split(/\s+/,$sample_id[$i]);
		$inputgvcf = $inputgvcf."--variant $original_file[1] ";
	}
	#############
	##command lines
	#############

	my @combineGVCFs_commandline;
	my @filename;

	if ($gvcf_intervals) 
	{
		$combineGVCFs_commandline[0] = qq/$hash_config{"gatk"} --java-options "-Xmx50g -XX:+UseParallelGC -XX:ParallelGCThreads=1" CombineGVCFs -L $gvcf_intervals -R $hash_config{"opt_ref_genome"} $inputgvcf -O $gvcf_dir\/combine.$gvcf_intervals.$jobname_pre.g.vcf.gz 1>>$log_dir\/combine.$jobname_pre.log 2>>$log_dir\/combine.$jobname_pre.log /;

		$filename[0] = "$gvcf_dir\/combine.$gvcf_intervals.$jobname_pre.g.vcf.gz";
	}
	else
	{
		
=a		
		$combineGVCFs_commandline[0] = qq/$hash_config{"gatk"} --java-options "-Xmx50g -XX:+UseParallelGC -XX:ParallelGCThreads=1" CombineGVCFs -R $hash_config{"opt_ref_genome"} $inputgvcf -O $gvcf_dir\/combine.$jobname_pre.g.vcf.gz 1>>$log_dir\/combine.$jobname_pre.log 2>>$log_dir\/combine.$jobname_pre.log /;

		$filename[0] = "$gvcf_dir\/combine.$jobname_pre.g.vcf.gz";
=cut
	
	
		for (my $i = 1; $i <= 22; $i++) 
		{		
			my $tmp = "chr".$i;

			if ($buildver eq "hg19") 
			{
				$tmp = $i;
			}

			$combineGVCFs_commandline[$i-1] = qq/$hash_config{"gatk"} --java-options "-Xmx16g" CombineGVCFs -R $hash_config{"opt_ref_genome"} $inputgvcf -O $gvcf_dir\/combine.$jobname_pre.$i.g.vcf.gz -L $tmp 1>>$log_dir\/combine.$jobname_pre.log 2>>$log_dir\/combine.$jobname_pre.log /;

			$filename[$i-1] = "$gvcf_dir\/combine.$jobname_pre.$i.g.vcf.gz";
		}

		my $chrX = "chrX";
		my $chrY = "chrY";
		my $chrM = "chrM";

		if ($buildver eq "hg19") 
		{
			$chrX = "X";
			$chrY = "Y";
			$chrM = "MT";
		}		

		$combineGVCFs_commandline[22] = qq/$hash_config{"gatk"} --java-options "-Xmx16g" CombineGVCFs -R $hash_config{"opt_ref_genome"} $inputgvcf -O $gvcf_dir\/combine.$jobname_pre.X.g.vcf.gz -L $chrX 1>>$log_dir\/combine.$jobname_pre.log 2>>$log_dir\/combine.$jobname_pre.log /;
		$filename[22] = "$gvcf_dir\/combine.$jobname_pre.X.g.vcf.gz";

		$combineGVCFs_commandline[23] = qq/$hash_config{"gatk"} --java-options "-Xmx16g" CombineGVCFs -R $hash_config{"opt_ref_genome"} $inputgvcf -O $gvcf_dir\/combine.$jobname_pre.Y.g.vcf.gz -L $chrY  1>>$log_dir\/combine.$jobname_pre.log 2>>$log_dir\/combine.$jobname_pre.log /;
		$filename[23] = "$gvcf_dir\/combine.$jobname_pre.Y.g.vcf.gz";

		$combineGVCFs_commandline[24] = qq/$hash_config{"gatk"} --java-options "-Xmx16g" CombineGVCFs -R $hash_config{"opt_ref_genome"} $inputgvcf -O $gvcf_dir\/combine.$jobname_pre.M.g.vcf.gz -L $chrM  1>>$log_dir\/combine.$jobname_pre.log 2>>$log_dir\/combine.$jobname_pre.log /;
		$filename[24] = "$gvcf_dir\/combine.$jobname_pre.M.g.vcf.gz";
		
	}

	#######################

	#if ($where eq "local") 
	{
		my $whichSubroutine = \&stand_command;

		if(@combineGVCFs_commandline)
		{
    		&multifork($job_run,$whichSubroutine,"combine-$jobname_pre",@combineGVCFs_commandline);
		}



	}
=a	
	elsif($where eq "sbatch")
	{
		my @bsub_options;

		if(@combineGVCFs_commandline)
		{
			@bsub_options=("16000",36, "combine-$jobname_pre", $email, $nthreads_perjob);
    		&cluster_bsub (\@combineGVCFs_commandline,\@bsub_options);
		}
		
	}

=cut	
	if (!$gvcf_intervals) 
	{
		my $whichSubroutine = \&stand_command;

		&merge_vcf("$gvcf_dir/combine.$jobname_pre.g.vcf.gz",@filename);
		&multifork(1,$whichSubroutine,"merge-$jobname_pre","$hash_config{'tabix'} -p vcf $gvcf_dir/combine.$jobname_pre.g.vcf.gz");

		&multifork(1,$whichSubroutine,"rm-$jobname_pre","rm -rf ".join(" ",@filename));
		&multifork(1,$whichSubroutine,"rm-$jobname_pre","rm -rf $gvcf_dir/combine.$jobname_pre.*.g.vcf.gz.tbi");

		$filename[0] = "$gvcf_dir/combine.$jobname_pre.g.vcf.gz";
	}


	print STDERR "\tcombineGVCFs is finished\n";

	print TO "\n####".&nowtime." finished\n\n\n";
	close TO;		

	return $filename[0];
}


sub genotypeGVCFs()
{
	my ($nthreads, $opt_out_dir, $hash_config, $input) = @_;

	$nthreads_perjob=1;

	if ($nthreads>24 && !$nthreads) 
	{
		$nthreads = 24;
	}

	my $job_run = int ($nthreads/$nthreads_perjob);

	##mkdir dir
	my $vcf_dir = "$opt_out_dir/vcf"; #vcffiles

	mkdir ($vcf_dir,0755) unless (-d $vcf_dir);

	####
	open(TO,">>$commands") or die();
	print TO "############################\n";
	print TO "\n###################\ngenotypeGVCFs starts running; ".&nowtime." -->\n\n\n";

	print "\n###################\ngenotypeGVCFs starts running; -->\n";
	print "\n***********check config files\n";
	&checkfile($hash_config{"opt_ref_genome"},$hash_config{"dbsnp"},$hash_config{"gatk"},$hash_config{"ssub_array"});
	print "***********\n";

	############

	&checkfile($input);

	#############
	##command lines
	#############

	my @genotypeGVCFs_commandline;
	my @filename;

##########
	my $vcf_chroms_dir = "$opt_out_dir/vcf/chroms";

	if ($gvcf_intervals) ## X, Y chromosome
	{
		$genotypeGVCFs_commandline[0] = qq/$hash_config{"gatk"} --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=1" GenotypeGVCFs -R $hash_config{"opt_ref_genome"} -V $input -O $vcf_dir\/allchroms.$gvcf_intervals.$jobname_pre.vcf.gz --dbsnp $hash_config{"dbsnp"} -L $gvcf_intervals /;

		if ($ploidy != 2) 
		{
			$genotypeGVCFs_commandline[0] .= " -ploidy $ploidy ";
		}

		$genotypeGVCFs_commandline[0] .= " 1>>$log_dir\/GenotypeGVCFs.log  2>>$log_dir\/GenotypeGVCFs.log ";

	}
	else
	{

=a		
		$genotypeGVCFs_commandline[0] = qq/$hash_config{"gatk"} --java-options "-Xmx8g" GenotypeGVCFs -R $hash_config{"opt_ref_genome"} -V $input -O $vcf_dir\/allchroms.$jobname_pre.vcf.gz --dbsnp $hash_config{"dbsnp"} 1>>$log_dir\/GenotypeGVCFs.log  2>>$log_dir\/GenotypeGVCFs.log /;

		$filename[0] = "$vcf_dir\/allchroms.$jobname_pre.vcf.gz";
=cut

		mkdir ($vcf_chroms_dir,0755) unless (-d $vcf_chroms_dir);

		for (my $i = 1; $i <= 22; $i++) 
		{
			my $tmp = "chr".$i;

			if ($buildver eq "hg19")
			{
				$tmp = $i;
			}

			$genotypeGVCFs_commandline[$i-1] = qq/$hash_config{"gatk"} --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=1" GenotypeGVCFs -R $hash_config{"opt_ref_genome"} -V $input -O $vcf_chroms_dir\/$jobname_pre.$i.vcf.gz --dbsnp $hash_config{"dbsnp"} -L $tmp  1>>$log_dir\/GenotypeGVCFs.log  2>>$log_dir\/GenotypeGVCFs.log /;

			$filename[$i-1] = "$vcf_chroms_dir\/$jobname_pre.$i.vcf.gz";
		}

		my $chrX = "chrX";
		my $chrY = "chrY";
		my $chrM = "chrM";

		if ($buildver eq "hg19") 
		{
			$chrX = "X";
			$chrY = "Y";
			$chrM = "MT";
		}	

		$genotypeGVCFs_commandline[22] = qq/$hash_config{"gatk"} --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=1" GenotypeGVCFs -R $hash_config{"opt_ref_genome"} -V $input -O $vcf_chroms_dir\/$jobname_pre.X.vcf.gz --dbsnp $hash_config{"dbsnp"} -L $chrX 1>>$log_dir\/GenotypeGVCFs.log  2>>$log_dir\/GenotypeGVCFs.log /;


		$genotypeGVCFs_commandline[23] = qq/$hash_config{"gatk"} --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=1" GenotypeGVCFs -R $hash_config{"opt_ref_genome"} -V $input -O $vcf_chroms_dir\/$jobname_pre.Y.vcf.gz --dbsnp $hash_config{"dbsnp"}  -L $chrY 1>>$log_dir\/GenotypeGVCFs.log 2>>$log_dir\/GenotypeGVCFs.log /;

		$genotypeGVCFs_commandline[24] = qq/$hash_config{"gatk"} --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=1" GenotypeGVCFs -R $hash_config{"opt_ref_genome"} -V $input -O $vcf_chroms_dir\/$jobname_pre.M.vcf.gz --dbsnp $hash_config{"dbsnp"}  -L $chrM 1>>$log_dir\/GenotypeGVCFs.log 2>>$log_dir\/GenotypeGVCFs.log /;

		$filename[22] = "$vcf_chroms_dir\/$jobname_pre.X.vcf.gz";
		$filename[23] = "$vcf_chroms_dir\/$jobname_pre.Y.vcf.gz";
		$filename[24] = "$vcf_chroms_dir\/$jobname_pre.M.vcf.gz";		
	}

	#######################


	if ($where eq "local") 
	{
		my $whichSubroutine = \&stand_command;

		if(@genotypeGVCFs_commandline)
		{
    		&multifork($job_run,$whichSubroutine,"vcf-$jobname_pre",@genotypeGVCFs_commandline);
		}

	}
	elsif($where eq "sbatch")
	{
		my @bsub_options;

		if(@genotypeGVCFs_commandline)
		{
			@bsub_options=("8000",24, "vcf-$jobname_pre", $email, $nthreads_perjob);
    		&cluster_bsub (\@genotypeGVCFs_commandline,\@bsub_options);
		}
		
	}


	my $whichSubroutine = \&stand_command;

	if (!$gvcf_intervals) 
	{
		&merge_vcf("$vcf_dir/allchroms.$jobname_pre.vcf.gz",@filename);
		&multifork(1,$whichSubroutine,"vcf-merge-$jobname_pre","$hash_config{'tabix'} -p vcf $vcf_dir/allchroms.$jobname_pre.vcf.gz && rm -rf $vcf_chroms_dir");

		#system("rm -rf $vcf_chroms_dir");
	}

	print STDERR "\tgenotypeGVCFs is finished\n";

	print TO "\n####".&nowtime." finished\n\n\n";
	close TO;		

	$startflag = 1;

	return "$vcf_dir/allchroms.$jobname_pre.vcf.gz";
}

###############
##vqsr
###############

sub vqsr()
{
	my ($nthreads, $opt_out_dir, $hash_config, $vcf) = @_;

	$nthreads_perjob = 1;

	if ($nthreads>24 && !$nthreads) 
	{
		$nthreads = 24;
	}

	my $job_run = int ($nthreads/$nthreads_perjob);

	##mkdir dir
	my $vcf_dir = "$opt_out_dir/vcf"; #vcffiles
	my $vqsr_dir = "$opt_out_dir/vcf/vqsr"; #vcffiles_vqsr

	mkdir ($vcf_dir,0755) unless (-d $vcf_dir);
	mkdir ($vqsr_dir,0755) unless (-d $vqsr_dir);
	####
	open(TO,">>$commands") or die();
	print TO "############################\n";
	print TO "\n###################\nVQSR starts running; ".&nowtime." -->\n\n\n";

	print "\n###################\nVQSR starts running; -->\n";
	print "\n***********check config files\n";
		
	&checkfile($hash_config{"opt_ref_genome"},$hash_config{"Mills_and_1000G_gold_standard_indels"},$hash_config{"dbsnp"},$hash_config{'vcf2candidate'},$hash_config{"gatk"},$hash_config{"ssub_array"},$hash_config{"hapmap"},$hash_config{"kg_omni"},$hash_config{"KG_phase1_snps"},$hash_config{"tabix"});
	print "***********\n";

	&checkfile($vcf);

	###check the index of vcf

	my $whichSubroutine = \&stand_command; 


	my $sample_num= `less $vcf |grep -v '##'| head -1 |wc`;
	my $flag_gatk_version= `less $vcf |head -2000|grep '##'| grep 'Version=3.3' `;

	chomp($sample_num);
	my @sample_num2=split(/\s+/,$sample_num);

	$sample_num = $sample_num2[2] - 9;

	my $depth_filter;

	if ($sample_num < 5) 
	{
		$depth_filter = $sample_num * 8;
	}
	else
	{
		$depth_filter = $sample_num * 2;
	}

	###command line

	my $vqsr = "";

	if (-e "$vcf.tbi") {}
	else
	{
		$vqsr .= qq/$hash_config{'tabix'} -p vcf $vcf /;
	}
	#$vqsr = qq/$hash_config{'vcf2candidate'} -basic_filter -vcf $vcf -out_dir $vcf_dir -output_prefix annovar.filter1 -where local -nthreads 1 -individual_dp 5 -individual_gq 10 -depth '\\>$depth_filter' 2>>$log_dir\/vqsr.log /;

	if ($vqsr)
	{
		$vqsr .= qq/ && /;
	}

	#$vqsr .= qq/ $hash_config{"gatk"} SelectVariants -select-type SNP -V $vcf -O $vcf_dir\/vqsr.snp.$jobname_pre.vcf.gz  1>>$log_dir\/vqsr.log 2>>$log_dir\/vqsr.log /;
	$vqsr .= qq/ $hash_config{"gatk"} VariantFiltration -V $vcf -O $vcf_dir\/vqsr.filter.$jobname_pre.vcf.gz --filter-expression "DP < $depth_filter" --filter-name "lowDP"  1>>$log_dir\/vqsr.log 2>>$log_dir\/vqsr.log /;

	# -an DP 
	$vqsr .= qq/ && $hash_config{"gatk"} --java-options "-Xmx16g -XX:+UseParallelGC -XX:ParallelGCThreads=1" VariantRecalibrator -R $hash_config{"opt_ref_genome"} -mode SNP -V $vcf_dir\/vqsr.filter.$jobname_pre.vcf.gz -an QD -an FS -an MQRankSum -an ReadPosRankSum -an MQ -an SOR --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hash_config{"hapmap"} --resource:omni,known=false,training=true,truth=false,prior=12.0 $hash_config{"kg_omni"} --resource:1000G,known=false,training=true,truth=false,prior=10.0 $hash_config{"KG_phase1_snps"} --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $hash_config{"dbsnp"}  -tranches-file $vqsr_dir\/$jobname_pre-gatk_snp.tranches -rscript-file $vqsr_dir\/$jobname_pre-gatk_snp.plots.R --tranche 90.0 --tranche 93.0 --tranche 95.0 --tranche 97.0 --tranche 99.0 --tranche 99.9 --tranche 100.0 -O $vqsr_dir\/$jobname_pre-gatk_snp.recal /; #

	#$vqsr[1] .= qq/-an InbreedingCoeff / if $sample_num >20;  ##version4.0 has not the option
	$vqsr .= qq/ -an DP / if $sequencing_stype eq "whole";
	#if ($flag_gatk_version !~ /\S+/){}else {$vqsr[1] .= qq/-an SOR /;}

	$vqsr .= qq/ 1>>$log_dir\/vqsr.log 2>>$log_dir\/vqsr.log /;

	$vqsr .= qq/ && $hash_config{"gatk"} --java-options "-Xmx16g -XX:+UseParallelGC -XX:ParallelGCThreads=1" ApplyVQSR -R $hash_config{"opt_ref_genome"} -mode SNP -V $vcf_dir\/vqsr.filter.$jobname_pre.vcf.gz --recal-file $vqsr_dir\/$jobname_pre-gatk_snp.recal --tranches-file $vqsr_dir\/$jobname_pre-gatk_snp.tranches -O $vcf_dir\/vqsr-t95-SNPrecal.$jobname_pre.vcf.gz --truth-sensitivity-filter-level 95  1>>$log_dir\/vqsr.log 2>>$log_dir\/vqsr.log /;


	#  -an DP
 
 $vqsr .= qq/ && $hash_config{"gatk"} --java-options "-Xmx16g -XX:+UseParallelGC -XX:ParallelGCThreads=1" VariantRecalibrator -R $hash_config{"opt_ref_genome"} -mode INDEL --max-gaussians 4 -V $vcf_dir\/vqsr-t95-SNPrecal.$jobname_pre.vcf.gz -an QD -an FS -an MQRankSum -an ReadPosRankSum -an SOR --resource:mills,known=true,training=true,truth=true,prior=12.0 $hash_config{"Mills_and_1000G_gold_standard_indels"} -O $vqsr_dir\/$jobname_pre-gatk_indel.recal -tranches-file $vqsr_dir\/$jobname_pre-gatk_indel.tranches -rscript-file $vqsr_dir\/$jobname_pre-gatk_indel.plots.R --tranche 90.0 --tranche 93.0 --tranche 95.0 --tranche 97.0 --tranche 99.0 --tranche 99.9 --tranche 100.0 /;#

 	#$vqsr[3] .= qq/-an InbreedingCoeff / if $sample_num >20;
 $vqsr .= qq/ -an DP / if $sequencing_stype eq "whole";
 	#if ($flag_gatk_version !~ /\S+/){}else {$vqsr[3] .= qq/-an SOR /;}

$vqsr .= qq/  1>>$log_dir\/vqsr.log 2>>$log_dir\/vqsr.log /;
 	
$vqsr .= qq/ && $hash_config{"gatk"} --java-options "-Xmx16g -XX:+UseParallelGC -XX:ParallelGCThreads=1" ApplyVQSR -R $hash_config{"opt_ref_genome"} -mode INDEL -V $vcf_dir\/vqsr-t95-SNPrecal.$jobname_pre.vcf.gz --recal-file $vqsr_dir\/$jobname_pre-gatk_indel.recal -tranches-file $vqsr_dir\/$jobname_pre-gatk_indel.tranches -O $vcf_dir\/vqsr-t95-bothrecal.$jobname_pre.vcf.gz --truth-sensitivity-filter-level 95.0  1>>$log_dir\/vqsr.log 2>>$log_dir\/vqsr.log /;

$vqsr .= qq/ &&  rm -rf $vcf_dir\/vqsr.filter.$jobname_pre.vcf.gz*  $vcf_dir\/vqsr-t95-SNPrecal.$jobname_pre.vcf.gz*  /;

#$vqsr .= qq/ && rm -rf $vcf_dir\/basic_filter $vcf_dir\/$jobname_pre-t95-SNPrecal.vcf.gz* /;

	###############################

	#if($where eq "local")
	{
		##run vqsr
		&multifork(1,$whichSubroutine,"vqsr-$jobname_pre",$vqsr);

		#(my $whichsteps, my $num) = &checktail_stats($stats_log);

		#die "VQSR is wrong" if $whichsteps ne "vqsr" || $num != @vqsr;
		#system("rm -rf $vcf_dir/$jobname_pre-t95-SNPrecal.vcf.gz*");		
	}

	print STDERR "\tVQSR is finished\n";

	print TO "\n####".&nowtime." finished\n\n\n";
	close TO;	

	return "$vcf_dir/vqsr-t95-bothrecal.$jobname_pre.vcf.gz";

}

###############
##annovar
###############

sub annovar()
{
	my ($nthreads, $opt_out_dir, $hash_config, $input_vcf) = @_;

	$nthreads_perjob = 1;

	if ($nthreads>24 && !$nthreads) 
	{
		$nthreads = 24;
	}

	my $job_run = int ($nthreads/$nthreads_perjob);

	##mkdir dir
	my $annovar_dir = "$opt_out_dir/vcf/annovar"; #annovar files
	my $annovar_origin_dir = "$opt_out_dir/vcf/annovar/original"; #annovar files

	eval{mkpath($annovar_origin_dir,0,0755)};
	if($@)
	{
	   warn("Make path [$annovar_origin_dir] failed:\n$@");
	}
	
	mkdir ($annovar_dir,0755) unless (-d $annovar_dir);
	mkdir ($annovar_origin_dir,0755) unless (-d $annovar_origin_dir);
	
	####
	open(TO,">>$commands") or die();
	print TO "############################\n";
	print TO "\n###################\nannovar starts running; ".&nowtime." -->\n\n\n";


	print "\n###################\nannovar starts running; -->\n";
	print "\n***********check config files\n";

	&checkfile($hash_config{"convert2annovar"},$hash_config{"annotate_variation"},$hash_config{"ssub_array"},$hash_config{"annovar_db"});
	print "***********\n";

	&checkfile($input_vcf);

	my @annovar_commandline1;
	my @annovar_commandline2;
	my @annovar_commandline3;


	my $vcffile_size = -s "$input_vcf";
	
	
		
		#$annovar_commandline1[0] = qq/less $input_vcf |grep -v '^#' |perl -lane ' print join(\"\\t\",\@F[0..8]) if (\$F[4] !~\/,\/)'  > $annovar_origin_dir\/$annovar.vcf && echo annovar.finished >>$stats_log/;

		 if ($input_vcf =~ /\.gz$/)
		 {
		 	$annovar_commandline1[0] = qq/ gunzip -c $input_vcf /;
		 }else {
			$annovar_commandline1[0] = qq/ less $input_vcf /;
		 }



		$annovar_commandline1[0] .= qq/ | grep -v '^#' | perl -lane 'if(\$F[4] =~ \/\,\/){ \@tmp = split("\,",\$F[4]); if(\$tmp[0] eq "*"){\$F[4]=\$tmp[1];}else{\$F[4]=\$tmp[0];} } print join("\\t",\@F[0..8])' > $annovar_origin_dir\/annovar.$jobname_pre.avi  && perl $hash_config{"convert2annovar"} -format vcf4old $annovar_origin_dir\/annovar.$jobname_pre.avi --outfile $annovar_origin_dir\/annovar.$jobname_pre.avi2 2>>$log_dir\/annovar.log && mv $annovar_origin_dir\/annovar.$jobname_pre.avi2 $annovar_origin_dir\/annovar.$jobname_pre.avi /;

		
		## annotation

		my $genetype = "refgene";
		my $ver1000g = "1000g2015aug_all";
		my $file1000g;

			if($ver1000g =~ m/^1000g(20\d\d)([a-z]{3})/) {
			my %monthhash = ('jan'=>'01', 'feb'=>'02', 'mar'=>'03', 'apr'=>'04', 'may'=>'05', 'jun'=>'06', 'jul'=>'07', 'aug'=>'08', 'sep'=>'09', 'oct'=>'10', 'nov'=>'11', 'dec'=>'12');
			$file1000g = $1 . '_' . $monthhash{$2};}

		my $verdbsnp = "avsnp150";
		my $ljb = "dbnsfp35a";
		my $veresp = "esp6500siv2_all";
		my $clinvar = "clinvar_20200316";
		my $cosmic = "cosmic91";
		my $verexac = "exac03nontcga";

		$annovar_commandline2[0] = "$hash_config{'annotate_variation'} -geneanno -buildver $buildver -dbtype $genetype -outfile $annovar_origin_dir\/annovar.$jobname_pre.vcf -exonsort $annovar_origin_dir\/annovar.$jobname_pre.avi $hash_config{'annovar_db'}  2>>$log_dir\/annovar.log ";

		$annovar_commandline2[1] = "$hash_config{'annotate_variation'} -regionanno -dbtype segdup -buildver $buildver -outfile $annovar_origin_dir\/annovar.$jobname_pre.vcf $annovar_origin_dir\/annovar.$jobname_pre.avi $hash_config{'annovar_db'}  2>>$log_dir\/annovar.log ";

		$annovar_commandline2[2] = "$hash_config{'annotate_variation'} -filter -dbtype $ver1000g -buildver $buildver -outfile $annovar_origin_dir\/annovar.$jobname_pre.vcf $annovar_origin_dir\/annovar.$jobname_pre.avi $hash_config{'annovar_db'}  2>>$log_dir\/annovar.log ";

		$annovar_commandline2[3] = "$hash_config{'annotate_variation'} -filter -dbtype $verdbsnp -buildver $buildver -outfile $annovar_origin_dir\/annovar.$jobname_pre.vcf $annovar_origin_dir\/annovar.$jobname_pre.avi $hash_config{'annovar_db'}  2>>$log_dir\/annovar.log ";

		$annovar_commandline2[4] = "$hash_config{'annotate_variation'} -filter -dbtype $ljb -buildver $buildver -outfile $annovar_origin_dir\/annovar.$jobname_pre.vcf $annovar_origin_dir\/annovar.$jobname_pre.avi $hash_config{'annovar_db'} -otherinfo  2>>$log_dir\/annovar.log ";

		$annovar_commandline2[5] = "$hash_config{'annotate_variation'} -filter -dbtype $veresp -buildver $buildver -outfile $annovar_origin_dir\/annovar.$jobname_pre.vcf $annovar_origin_dir\/annovar.$jobname_pre.avi $hash_config{'annovar_db'}  2>>$log_dir\/annovar.log ";

		##check
		$annovar_commandline2[6] = "$hash_config{'annotate_variation'} -filter -dbtype $clinvar -buildver $buildver -outfile $annovar_origin_dir\/annovar.$jobname_pre.vcf $annovar_origin_dir\/annovar.$jobname_pre.avi $hash_config{'annovar_db'} -otherinfo 2>>$log_dir\/annovar.log ";

		$annovar_commandline2[7] = "$hash_config{'annotate_variation'} -filter -dbtype $cosmic -buildver $buildver -outfile $annovar_origin_dir\/annovar.$jobname_pre.vcf $annovar_origin_dir\/annovar.$jobname_pre.avi $hash_config{'annovar_db'}  2>>$log_dir\/annovar.log ";

		$annovar_commandline2[8] = "$hash_config{'annotate_variation'} -regionanno -dbtype cytoBand -buildver $buildver -outfile $annovar_origin_dir\/annovar.$jobname_pre.vcf $annovar_origin_dir\/annovar.$jobname_pre.avi $hash_config{'annovar_db'}  2>>$log_dir\/annovar.log ";

		$annovar_commandline2[9] = "$hash_config{'annotate_variation'} -regionanno -dbtype gwasCatalog -buildver $buildver -outfile $annovar_origin_dir\/annovar.$jobname_pre.vcf $annovar_origin_dir\/annovar.$jobname_pre.avi $hash_config{'annovar_db'}  2>>$log_dir\/annovar.log ";

		$annovar_commandline2[10] = "$hash_config{'annotate_variation'} -filter -dbtype $verexac -buildver $buildver -outfile $annovar_origin_dir\/annovar.$jobname_pre.vcf $annovar_origin_dir\/annovar.$jobname_pre.avi $hash_config{'annovar_db'} -otherinfo 2>>$log_dir\/annovar.log ";
	

	#############
	##command lines
	#############

	my $whichSubroutine = \&stand_command;
	&multifork($job_run,$whichSubroutine,"annovar1-$jobname_pre",@annovar_commandline1);


	if($where eq "local")
	{
		##run annovar_commandline

		&multifork($job_run,$whichSubroutine,"annovar2-$jobname_pre",@annovar_commandline2);
	
	}elsif($where eq "sbatch")
	{
		my @bsub_options;

		@bsub_options=("8000",24, "annovar2-$jobname_pre", $email, $nthreads_perjob);
    	&cluster_bsub (\@annovar_commandline2,\@bsub_options);
		
	}

##############summary

	my $outfile = "$annovar_origin_dir\/annovar.$jobname_pre.vcf";

	open (FUNCTION, "$outfile.variant_function") or die "Error: cannot read from variant function file: $!\n";
	open (STEP0, "$outfile.exonic_variant_function") or die "Error: cannot read from exonic variant function file: $!\n";
	open (STEP1, "$outfile.${buildver}_genomicSuperDups") or die "Error: cannot read from segdup file: $!\n";
	open (STEP2, "$outfile.${buildver}_ALL.sites.${file1000g}_dropped") or die "Error: cannot read from drop file $outfile.${buildver}_ALL.sites.${file1000g}_dropped: $!\n";
	open (STEP3, "$outfile.${buildver}_${verdbsnp}_dropped") or die "Error: cannot read from snp$verdbsnp drop file: $!\n";
	open (STEP4, "$outfile.${buildver}_${ljb}_dropped") or die "Error: cannot read from ${ljb} drop file: $!\n";
	open (STEP5, "$outfile.${buildver}_${veresp}_dropped") or die "Error: cannot read from esp${veresp}_all drop file: $!\n";
	open (STEP6, "$outfile.${buildver}_${clinvar}_dropped") or die "Error: cannot read from $clinvar drop file: $!\n";
	open (STEP7, "$outfile.${buildver}_${cosmic}_dropped") or die "Error: cannot read from $cosmic drop file: $!\n";
	open (STEP8, "$outfile.${buildver}_cytoBand") or die "Error: cannot read from cytoBand file: $!\n";
	open (STEP9, "$outfile.${buildver}_gwasCatalog") or die "Error: cannot read from gwasCatalog file: $!\n";
	open (STEP10, "$outfile.${buildver}_${verexac}_dropped") or die "Error: cannot read from exac file: $!\n";


####################
##
###########################
	my (@allstep);

	while (<STEP0>)  {
		m/^line\d+\t([^\t]+)\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile exonic: <$_>\n";
		my ($efun, $aachange, $varstring) = ($1, $2, $3);

		($efun) = split(/\s+/,$efun);
		#my @aachange = split (/:|,/, $aachange);
		$allstep[0]->{$varstring} = "$efun,"."\"$aachange\""; #aachange could be "UNKNOWN"	
	}

	while (<STEP1>) {
		m/^segdup\tScore=(\S+);\S+\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile segdup : <$_>\n";
		$allstep[1]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
	}

	while (<STEP2>) {
		m/^1000g\w+_all\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 1000g : <$_>\n";
		$allstep[2]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
	}

	while (<STEP3>) {
		m/^avsnp\d+\w+?\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile avsnp : <$_>\n";
		$allstep[3]->{$2} = $1;
	}

	while (<STEP4>) {
		m/^${ljb}\t([^\t]+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile ljb : <$_>\n";
		$allstep[4]->{$2} = $1;
	}

	while (<STEP5>) {
		m/^esp\d+\w*_all\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile esp : <$_>\n";
		$allstep[5]->{$2} = $1;
	}

#?
	while (<STEP6>) {
		m/^clinvar_\d+\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile clinvar : <$_>\n";
		$allstep[6]->{$2} = "\"$1\"";
	}

	while (<STEP7>) {
		m/^cosmic\d+\tID=(\S+);OCCURENCE=(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile cosmic : <$_>\n";
		$allstep[7]->{$3} = "\"$1\",\"$2\"";
	}

	while (<STEP8>) {
		m/^cytoBand\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile cytoBand : <$_>\n";
		$allstep[8]->{$2} = "\"$1\"";
	}

	while (<STEP9>) {
		m/^gwasCatalog\tName=([^\t]+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile gwasCatalog : <$_>\n";
		#my $gwas = $1;
		#$gwas =~ s/\s+/\_/g;
		$allstep[9]->{$2} = "\"$1\"";
	}

	while (<STEP10>) {
		m/^exac\S+\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile exac : <$_>\n";
		$allstep[10]->{$2} = $1;
	}

############################
##write
#############################
open (OUT, ">$outfile.genome_summary.csv") or die "Error: cannot write to output file: $!\n";

print OUT join (',', qw/Chr Start End Ref Alt Func Gene ExonicFunc AAChange SegDup cytoBand/, "ESP_ALL", "KG_ALL", "dbSNP", "ExAC_Freq", "ExAC_AFR", "ExAC_AMR", "ExAC_EAS", "ExAC_FIN", "ExAC_NFE", "ExAC_OTH", "ExAC_SAS", "SIFT_score","SIFT_converted_rankscore","SIFT_pred","Polyphen2_HDIV_score","Polyphen2_HDIV_rankscore","Polyphen2_HDIV_pred","Polyphen2_HVAR_score","Polyphen2_HVAR_rankscore","Polyphen2_HVAR_pred","LRT_score","LRT_converted_rankscore","LRT_pred","MutationTaster_score","MutationTaster_converted_rankscore","MutationTaster_pred","MutationAssessor_score","MutationAssessor_score_rankscore","MutationAssessor_pred","FATHMM_score","FATHMM_converted_rankscore","FATHMM_pred","PROVEAN_score","PROVEAN_converted_rankscore","PROVEAN_pred","VEST3_score","VEST3_rankscore","MetaSVM_score","MetaSVM_rankscore","MetaSVM_pred","MetaLR_score","MetaLR_rankscore","MetaLR_pred","M-CAP_score","M-CAP_rankscore","M-CAP_pred","REVEL_score","REVEL_rankscore","MutPred_score","MutPred_rankscore","CADD_raw","CADD_raw_rankscore","CADD_phred","DANN_score","DANN_rankscore","fathmm-MKL_coding_score","fathmm-MKL_coding_rankscore","fathmm-MKL_coding_pred","Eigen_coding_or_noncoding","Eigen-raw","Eigen-PC-raw","GenoCanyon_score","GenoCanyon_score_rankscore","integrated_fitCons_score","integrated_fitCons_score_rankscore","integrated_confidence_value","GERP++_RS","GERP++_RS_rankscore","phyloP100way_vertebrate","phyloP100way_vertebrate_rankscore","phyloP20way_mammalian","phyloP20way_mammalian_rankscore","phastCons100way_vertebrate","phastCons100way_vertebrate_rankscore","phastCons20way_mammalian","phastCons20way_mammalian_rankscore","SiPhy_29way_logOdds","SiPhy_29way_logOdds_rankscore","Interpro_domain","GTEx_V6p_gene","GTEx_V6p_tissue", "clinvar", "cosmic", "Occurence","gwasCatalog"),"\n";

while (<FUNCTION>) {
	s/[\r\n]+$//;
	m/^(\S+)\t([^\t]+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)(.*)/ or die "Error: invalid record found in annovar outputfile: <$_>\n";
	my ($function, $gene, $varstring, $otherinfo) = ($1, $2, $3, $4||'');


	my @varstring = split (/\s+/, $varstring);

	print OUT join (',', @varstring).",";

	print OUT qq/"$function","$gene"/;
	
	#exon
	
	if (defined $allstep[0]->{$varstring}) {	
		print OUT qq/,$allstep[0]->{$varstring}/;
	}else {
		print OUT ",.,.";
	}

	#segdup
	
	if (defined $allstep[1]->{$varstring}) {	
		print OUT qq/,$allstep[1]->{$varstring}/;
	}else {
		print OUT ",.";
	}
		

	#cytoBand
	if (defined $allstep[8]->{$varstring}) {
		print OUT qq/,$allstep[8]->{$varstring}/;
	} else {
		print OUT ",.";
	}
	######################################
	##_esp
	if (defined $allstep[5]->{$varstring}) {
		print OUT qq/,$allstep[5]->{$varstring}/;
	} else {
		print OUT ",.";
	}

	#2,:1000G

	if (defined $allstep[2]->{$varstring}) {
		print OUT qq/,$allstep[2]->{$varstring}/;
	} else {
		print OUT ",.";
	}
	#3:avsnp
	if (defined $allstep[3]->{$varstring}) {
		print OUT qq/,$allstep[3]->{$varstring}/;
	} else {
		print OUT ",.";
	}

#############################

	##exac
	if (defined $allstep[10]->{$varstring}) {
		print OUT qq/,$allstep[10]->{$varstring}/;
	} else {
		print OUT ",.,.,.,.,.,.,.,.";
	}

	#ljb
	if (defined $allstep[4]->{$varstring}) {
		print OUT qq/,$allstep[4]->{$varstring}/;
	} else {
		print OUT ",.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.";
	}

	#clinvar
	if (defined $allstep[6]->{$varstring}) {
		
		my $clinvar2 = $allstep[6]->{$varstring};

		print OUT qq/,$clinvar2/;
	} else {
		print OUT ",.";
	}


	#cosmic
	if (defined $allstep[7]->{$varstring}) {
		print OUT qq/,$allstep[7]->{$varstring}/;
	} else {
		print OUT ",.,.";
	}

	#gwasCatalog
	if (defined $allstep[9]->{$varstring}) {
		print OUT qq/,$allstep[9]->{$varstring}/;
	} else {
		print OUT ",.";
	}

	print OUT "\n";
	
}
		
	#$annovar_summary[0] = "$annovar_origin_dir/annovar.$jobname_pre.vcf.genome_summary.csv";
		
	#$annovar_commandline[2] = qq/perl $hash_config{"summarize_annovar"} $annovar_origin_dir\/$annovar.avi $hash_config{"annovar_db"} -buildver hg19 --verdbsnp 138 --ver1000g 1000g2014oct -veresp 6500siv2 -clinvar clinvar_20140929 -ljb ljb26 -cosmic cosmic70 -outfile $annovar_origin_dir\/$annovar.vcf 2>>$log_dir\/annovar.log && echo annovar.finished >>$stats_log /;

	$annovar_commandline3[0] = qq/$hash_config{"annovar2vcf"} -v $input_vcf -a $outfile.genome_summary.csv -o $annovar_dir\/annovar.$jobname_pre.vcf.gz 2>>$log_dir\/annovar.log && $hash_config{'tabix'} -p vcf $annovar_dir\/annovar.$jobname_pre.vcf.gz /; #-r repeat.txt 

	&multifork($job_run,$whichSubroutine,"annovar3-$jobname_pre",@annovar_commandline3);

	print STDERR "\tannovar_commandline is finished\n";

	print TO "\n####".&nowtime." finished\n\n\n";
	close TO;	

	$startflag = 1;

	return "$annovar_dir/annovar.$jobname_pre.vcf.gz";
}


##############################################################################################################
### basic subprogress
################################################################################################################

sub nowtime()
{
	my @time=(localtime)[5,4,3,2,1,0];   
	$time[0]+=1900;   
	$time[1]+=1;
	my $nowtime=sprintf("%04u-%02u-%02uT%02u:%02u:%02u",@time);  
	return $nowtime;
}

sub cluster_bsub()
{
	my ($command_line,$bsub_options) = @_;

	my $jobname_pre = @$bsub_options[2];
	my $signal_total = @$command_line;
	my $opt_working_dir = `pwd`; chomp($opt_working_dir);


	my $cluster_submit;

	if ($where eq "qsub") 
	{
		$cluster_submit = $hash_config{"bsub"};
		open(CLUSTER, "|$cluster_submit -pd `pwd` -pm @$bsub_options[0] -ph @$bsub_options[1] -jobn @$bsub_options[2] -email @$bsub_options[3] -e");
	}
	elsif($where eq "sbatch")
	{
		$cluster_submit = $hash_config{"ssub_array"};
		open(CLUSTER, "|$cluster_submit -pd `pwd` -pm @$bsub_options[0] -ph @$bsub_options[1] -cpt @$bsub_options[4] -jobn @$bsub_options[2] -email @$bsub_options[3] -e");
	}

	print CLUSTER join("\n",@$command_line)."\n";

	close CLUSTER;

	while()
	{
		my $signal = 0;	

		if(-e "$opt_working_dir/sbatch.$jobname_pre.jobs.signal")
		{
			open IN, "<$opt_working_dir/sbatch.$jobname_pre.jobs.signal";

			while(my $line=<IN>)
			{
				chomp($line);

				next if $line !~ /\S+/;
				next if $line =~ /^#/;

				my @f = split(/\./,$line);

				if ($f[0] eq $jobname_pre) 
				{
					$signal++;
				}		 
			}
			close IN;

			if ($signal == $signal_total) 
			{
				last;
			}
		}

		sleep($sleep);
	}

	my $out_file = ` less $opt_working_dir/tmp/$jobname_pre*.out | grep error `;

	 if ($out_file) 
	 {
	 	die "error: something wrong\n";
	 }

}

sub get_config()
{
	(my $get_config) =@_;

	my %hash_config;

	open(FROM,"$get_config") or die "there is not exist $get_config\n";
	while(my $line=<FROM>)
	{	
		chomp($line);

		next if $line =~ /^\s*$/ || $line =~ /^#/;
		my @tmp=split(/\s+/,$line);

		if($tmp[0]){$hash_config{$tmp[0]} = $tmp[1];}
	}
	close FROM;

	return %hash_config;
}

sub get_index()
{
	(my $get_index) = @_;

	my @sample_id;
	my $i=0;

	open(FROM,"$get_index") or die "there is not exist $get_index\n";
	while(my $line=<FROM>)
	{	
		chomp($line);

		next if $line =~ /^\s*$/ || $line =~ /^#/;
		my @tmp=split(/\s+/,$line);

		$sample_id[$i] = join("\t",@tmp);

		#&checkfile(@tmp[1..@tmp-1]);
		$i++;
	}
	close FROM;

	return (@sample_id);
}

sub get_index_check()
{
	(my $get_index) = @_;

	my @sample_id;
	my $i=0;

	open(FROM,"$get_index") or die "there is not exist $get_index\n";
	while(my $line=<FROM>)
	{	
		chomp($line);

		next if $line =~ /^\s*$/ || $line =~ /^#/;
		my @tmp=split(/\s+/,$line);

		$sample_id[$i] = join("\t",@tmp);

		&checkfile(@tmp[1..@tmp-1]);
		$i++;
	}
	close FROM;

	return (@sample_id);
}

sub multifork
{
	(my $num, my $whichSubroutine, my $flag, my @array) = @_;

	my $j=0;
	my $thread;

	my $array_num = @array;

	while()
	{
		last if($j>=$array_num);
		while(scalar(threads->list())<$num && $j < $array_num)  
		{		
			threads->new(
						$whichSubroutine,$array[$j]." && echo $flag.$j.finished >> $process "
				);
			$j++;
		}

		foreach $thread(threads->list(threads::all))
		{
			if($thread->is_joinable())
			{
				$thread->join();
				#print scalar(threads->list()),"\t$j\t",localtime(time),"\n";
			}
			if (my $err = $thread->error()) 
			{
				warn("Thread error: $err\n");
			}

		}
		sleep($sleep);		
	}

	foreach $thread(threads->list(threads::all))
	{
		$thread->join();

		if (my $err = $thread->error()) 
		{
			warn("Thread error: $err\n");
		}
		#print scalar(threads->list()),"\t$j\t",localtime(time),"\n";       
	}
}

sub stand_command()
{
	(my $cmd) = @_;

	my $run_status = 0;

	#print join("\t",$j,@array-1+1,$cmd)."\n";
	
	open(TO,">>$commands") or die();
	print TO "$cmd\n";
	close TO;

	$run_status = system($cmd);

	if ($run_status != 0)
	{
		&sendmail("$cmd");

		open(TO,">>$error") or die();
		print TO &nowtime."\n************\n$cmd\n is wrong\n************\n";
		close TO;

		die "\n************\n$cmd\n is wrong\n************\n";
	}
}


sub checkfile()
{
	my @file = @_;

	for (my $i = 0; $i < @file; $i++) 
	{		
		print STDERR "\tthe $file[$i] file exists -->";
		if (-e "$file[$i]") 
		{
			print STDERR "\tOK\n";
		}
		else 
		{
			die "\tERROR\n################\n$file[$i] file does not exist\nPlease check the original files\n#################\n";
		}
	}

}


sub sendmail()
{
	(my $command_line) = @_;

	if ($email)
	{
		my $from="vmpsched\@vmpsched.vampire";
		my $to="$email";
		my $subject="An error";

		my $sendmailpath="/usr/sbin/sendmail";

		my $message = "An error has occurred processing your job, see below.\n$command_line\n\nfrom cgg lab\n";

		open (SENDMAIL, "| $sendmailpath -t") or die "Cannot open $sendmailpath: $!";

		print SENDMAIL "Subject: $subject\n";
		print SENDMAIL "From: $from\n";
		print SENDMAIL "To: $to\n\n";

		print SENDMAIL "$message";

		close (SENDMAIL);
	}
}

sub merge_vcf()
{
	(my $output_vcf, my @filename) = @_;

	open(OUT, "| $hash_config{'bgzip'} -c > $output_vcf ") or die();

	my $i=0; 

	for (@filename)  
	{
		open FROM1, "<$_" if $_ !~ /\.gz$/;
		open FROM1, "gunzip -c $_|" if $_ =~ /\.gz$/;

		while(my $line=<FROM1>)
		{
			chomp($line);
			my @f = split(/\s+/, $line);
			if ($f[0] =~ /^#/)
			{
				if ($i==0){print OUT "$line\n";}
				next;
			}			
			print OUT $line."\n";	
		}
		close(FROM1);
		$i++;	
	}
	close(OUT);	
}

sub options
{
		
	print "\n*******************************************************\n";	
	print "*************************\n";

	print "\n*** The script was complied on Feb 10 2015 15:21:20 ***\n";
	print "\tUsage: script calling genotypes from fastq to vcf files [OPTIONS]\n";

	print "Options\n";
	print "\t-buildver\t[$buildver]\n";
	print "\t-get_config\t[$get_config]\n";
	print "\t-input\t[$input]\n";
	print "\t-where\t[$where]\n";
	print "\t-steps\t[$steps]\n";	
	print "\t-email\t[$email]\n";	
	print "\t-nthreads\t[$nthreads]\n" if $nthreads && $where eq "local";
	#print "\t-vcf\t[$vcf]\n";



	print "Additional Options\n";
	#print "\t-size\t[$fastq_size]\n";	
	#print "\t-gvcf_intervals\t[$gvcf_intervals]\n";
	#print "\t-Q64\t[$q64]\n";
	#print "\t-sequencing_stype\t[$sequencing_stype]\n";
	print "\t-sleep\t[$sleep]\n";


	print "Output\n";
	print "\t-out_dir\t[$opt_out_dir]\n";

	print "\n*************************\n";
	print "*******************************************************\n";	

}


sub usage()
{
		
	print "\n*******************************************************\n";	
	print "*************************\n";

	print "\n*** The script was complied on Feb 10 2015 15:21:20 ***\n";

	print "\tUsage: script calling genotypes from fastq to vcf files [OPTIONS]\n";

	print "Options\n";
	print "\t-get_config\t[]\n";
	print "\t-input\t[]\n";
	print "\t-where\t[]\n";
	print "\t-steps\t[]\n";	
	print "\t-email\t[]\n";	
	print "\t-nthreads\t[]\n";
	print "\t-vcf\t[]\n";


	print "Additional Options\n";
	print "\t-gvcf_intervals\t[]\n";
	print "\t-Q64\t[]\n";
	print "\t-sequencing_stype\t[]\n";
	print "\t-sleep\t[]\n";

	print "Output\n";
	print "\t-out_dir\t[]\n";

	print "\n*************************\n";
	######

	print "\n\trun -fastq2gvcf in -local:\n";
	print "\t\tExample 1.1: variantcalling.pl -fastq2gvcf -input bwa.txt -output_dir ab -where local -nthreads 8\n";
	print "\t\tExample 1.1.1: variantcalling.pl -fastq2gvcf -input bwa.txt -output_dir ab -where local -nthreads 8 -nthreads_perjob 4 \n";
	print "\t\tExample 1.2: variantcalling.pl -fastq2gvcf -input bwa.txt -get_config calling.config -output_dir ab -where local -nthreads 8 -jobname_pre bwa\n"; #-get_config /scratch/cgg/weiq1/script/from_fastq_to_dedup_realign_recalibra/smallfiles/GATK3.3.0/calling.config
	print "\t\tExample 1.3: variantcalling.pl -fastq2gvcf -input bwa.txt -gvcf_intervals 22 -output_dir ab -where local -nthreads 8 \n";
	print "\t\tExample 1.4: variantcalling.pl -fastq2gvcf -input bwa.txt -steps 1,2,4 -output_dir ab -where local -nthreads 8 \n";
	#print "\t\tExample 1.5: variantcalling.pl -fastq2gvcf -input bwa.txt -q64 true -output_dir ab -where local -nthreads 8 ##must contain step 3\n";
	print "\t\tExample 1.6: variantcalling.pl -fastq2gvcf -input bwa.txt -output_dir ab -where local -nthreads 8 -email *********@********* #will send email if there is an error\n";
	#print "\t\tExample 1.7: variantcalling.pl -fastq2gvcf -input bwa.txt -output_dir ab -where local -nthreads 8 \n";
	print "\t\tExample 1.8: variantcalling.pl -fastq2gvcf -input bwa2.txt -steps 2,3,4 -output_dir ab -where local -nthreads 8 \n";
	print "\t\tExample 1.9: variantcalling.pl -fastq2gvcf -input bwa2.txt -steps 4 -output_dir ab -where local -nthreads 8 \n";
	print "\t\tExample 1.10: variantcalling.pl -fastq2gvcf -input bwa2.txt -steps 2,4 -output_dir ab -where local -nthreads 8 \n";
	print "\t\tExample 1.11: variantcalling.pl -fastq2gvcf -qplot -combinegvcf -genotypeGVCFs -vqsr -annovar -input bwa.txt -output_dir ab -where local -nthreads 8 -sequencing_stype exome -jobname_pre 100  -steptime '13,2,11,15' -gvcf_intervals 22 -sleep 3600\n"; #-beagle4
	print "\t\tExample 1.12: variantcalling.pl -fastq2gvcf -combinegvcf -genotypeGVCFs -vqsr -input bwa.txt -output_dir ab -where local -nthreads 8\n";
	print "\t\tExample 1.13: variantcalling.pl -fastq2gvcf -qplot -input bwa.txt -output_dir ab -where local -nthreads 8 -sequencing_stype exome \n";
	print "\t\tExample 1.14: variantcalling.pl -fastq2gvcf -combinegvcf -genotypeGVCFs -annovar -input bwa.txt -output_dir ab -where local -nthreads 8\n";
	#print "\t\tExample 1.15: variantcalling.pl -fastq2gvcf -combinegvcf -genotypeGVCFs -beagle4 -input bwa.txt -output_dir ab -where local -nthreads 8\n";

	####
	print "\n\trun -fastq2gvcf in -qsub or -sbatch:\n";
	print "\t\tExample 2.1: variantcalling.pl -fastq2gvcf -input bwa.txt -output_dir ab -where sbatch -email *********@********* \n";
	print "\t\tExample 2.2: variantcalling.pl -fastq2gvcf -input bwa.txt -get_config calling.config -output_dir ab -where sbatch -email *********@*********  -jobname_pre bwa\n"; #-get_config /scratch/cgg/weiq1/script/from_fastq_to_dedup_realign_recalibra/smallfiles/GATK3.3.0/calling.config
	print "\t\tExample 2.3: variantcalling.pl -fastq2gvcf -input bwa.txt -gvcf_intervals 22 -output_dir ab -where sbatch -email *********@*********  \n";
	print "\t\tExample 2.4: variantcalling.pl -fastq2gvcf -input bwa.txt -steps 1,2,4 -output_dir ab -where sbatch -email *********@*********  \n";
	#print "\t\tExample 2.5: variantcalling.pl -fastq2gvcf -input bwa.txt -q64 true -output_dir ab -where sbatch -email *********@*********  ##must contain step 3\n";
	print "\t\tExample 2.6: variantcalling.pl -fastq2gvcf -input bwa.txt -output_dir ab -where sbatch -email *********@*********  #will send email if there is an error\n";
	#print "\t\tExample 2.7: variantcalling.pl -fastq2gvcf -input bwa.txt -output_dir ab -where sbatch -email *********@*********  \n";
	print "\t\tExample 2.8: variantcalling.pl -fastq2gvcf -input bwa2.txt -steps 2,3,4 -output_dir ab -where sbatch -email *********@*********  \n";
	print "\t\tExample 2.9:  variantcalling.pl -fastq2gvcf -input bwa2.txt -steps 4 -output_dir ab -where sbatch -email *********@*********  \n";
	print "\t\tExample 2.10:  variantcalling.pl -fastq2gvcf -input bwa2.txt -steps 2,4 -output_dir ab -where sbatch -email *********@*********  \n";
	print "\t\tExample 2.11:  variantcalling.pl -fastq2gvcf -qplot -combinegvcf -genotypeGVCFs -vqsr -annovar -input bwa.txt -output_dir ab -where sbatch -email *********@*********  -sequencing_stype exome -jobname_pre 100 -steptime '13,2,11,15' -sleep 30\n";
	print "\t\tExample 2.12:  variantcalling.pl -fastq2gvcf  -combinegvcf -genotypeGVCFs -vqsr -input bwa.txt -output_dir ab -where sbatch -email *********@********* \n";
	print "\t\tExample 2.13:  variantcalling.pl -fastq2gvcf -qplot -input bwa.txt -output_dir ab -where sbatch -email *********@*********  -sequencing_stype exome \n";
	print "\t\tExample 2.14:  variantcalling.pl -fastq2gvcf  -combinegvcf -genotypeGVCFs -annovar -input bwa.txt -output_dir ab -where sbatch -email *********@********* \n";
	#print "\t\tExample 2.15:  variantcalling.pl -fastq2gvcf -genotypeGVCFs -beagle4 -input bwa.txt -output_dir ab -where sbatch -email *********@********* \n";

	####
	print "\n\trun -genotypeGVCFs in -local:\n";
	print "\t\tExample 3.1:  variantcalling.pl -genotypeGVCFs -input combine.g.vcf.gz -output_dir ab -where local -nthreads 20 \n";
	print "\t\tExample 3.2:  variantcalling.pl -genotypeGVCFs -input combine.g.vcf.gz -output_dir ab -gvcf_intervals 22 -where local -nthreads 20 \n";	
	print "\t\tExample 3.3:  variantcalling.pl -genotypeGVCFs -vqsr -input combine.g.vcf.gz -output_dir ab -where local -nthreads 20 \n";
	print "\t\tExample 3.4:  variantcalling.pl -genotypeGVCFs -vqsr -annovar -input combine.g.vcf.gz -output_dir ab -where local -nthreads 20  \n";

####
	print "\n\trun -vqsr in -local:\n";
	print "\t\tExample 4.1:  variantcalling.pl -vqsr -input **.vcf.gz -output_dir ab -where local -nthreads 5 \n";
	print "\t\tExample 4.2:  variantcalling.pl -vqsr -annovar -input **.vcf.gz -output_dir ab -where local -nthreads 5  \n";
	#print "\t\tExample 4.3:  variantcalling.pl -vqsr -annovar -input **.vcf.gz -output_dir ab -where local -nthreads 5 \n";

####
	print "\n\trun -qplot in -local:\n";
	print "\t\tExample 5.1:  variantcalling.pl -qplot -sequencing_stype exome -input bwa2.txt -output_dir ab -where local -nthreads 5 \n";
	print "\t\tExample 5.2:  variantcalling.pl -qplot -sequencing_stype whole -input bwa2.txt -output_dir ab -where local -nthreads 5 \n";
	print "\t\tExample 5.2:  variantcalling.pl -qplot -sequencing_stype cpg -input bwa2.txt -output_dir ab -where local -nthreads 5 \n";
	print "\t\tExample 5.2:  variantcalling.pl -qplot -sequencing_stype *.bed -input bwa2.txt -output_dir ab -where local -nthreads 5 ## -dup_keep '--dup_keep' \n";

####
	print "\n\trun -annovar in -local:\n";
	print "\t\tExample 6.1:  variantcalling.pl -annovar -input **.vcf.gz -output_dir ab -where local -nthreads 5 \n";
	#print "\t\tExample 6.2:  variantcalling.pl -annovar -input **.vcf.gz -output_dir ./ -where local -nthreads 5 \n";
	#print "\t\tExample 6.3:  variantcalling.pl -annovar -beagle4 -vcf **.vcf.gz -output_dir ab -where local -nthreads 5  \n";
####
	print "\n\trun -beagle4 in -local:\n";
	#print "\t\tExample 7.1:  variantcalling.pl -beagle4 -vcf **.vcf.gz -output_dir ab -where local -nthreads 5\n";

	####
	print "\n\trun -combinegvcf in -local:\n";
	print "\t\tExample 8.1:  variantcalling.pl -combinegvcf -input bwa3.txt -output_dir ab -where local -nthreads 1  \n";
	print "\t\tExample 8.2:  variantcalling.pl -combinegvcf -input bwa3.txt -output_dir ab -where local -nthreads 1 -gvcf_intervals 1:100-200 \n";

	print "\t\tX-man: variantcalling.pl -fastq2gvcf -combinegvcf -genotypeGVCFs -annovar -input bwa3.txt -output_dir ./ -where local -sequencing_stype exome -nthreads 20 -steps 4 -gvcf_intervals X -ploidy 1 \n";
	####
	print "\n\tthe whole process:\n";
	print "\t\tlocal:  variantcalling.pl -fastq2gvcf -qplot -combinegvcf -genotypeGVCFs -vqsr -annovar -input bwa.txt -output_dir ./ -where local -sequencing_stype exome  -nthreads 14 \n";

	print "\t\tlocal: variantcalling.pl -fastq2gvcf -qplot -combinegvcf -genotypeGVCFs -vqsr -annovar -input bwa.txt -output_dir ./ -where local -sequencing_stype exome  -nthreads 14 -steps 1,2,4 \n"; #

	print "\t\tlocal: variantcalling.pl -fastq2gvcf -qplot -combinegvcf -genotypeGVCFs -vqsr -annovar -input bwa.txt -output_dir ./ -where local -sequencing_stype exome  -nthreads 14 -buildver hg19 # default -buildver hg38\n"; #


	print "\t\tcluster: variantcalling.pl -fastq2gvcf -qplot  -combinegvcf -genotypeGVCFs -vqsr -annovar -input bwa.txt -output_dir ./ -where sbatch -sequencing_stype exome -email *********@********* ## -nthreads_perjob 2 \n";


	####
	print "\n*************************\n";
	print "*******************************************************\n\n";	

	exit(0);
}


