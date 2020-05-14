#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Path;

############
##Parameter
############

my $opt_input_vcf = "";
my $input_beagle4 = "";
my $opt_out = "";

GetOptions("help|h",
	"input|input=s" => \$opt_input_vcf,
	"input_beagle4|input_beagle4=s" => \$input_beagle4,

	"output|output=s" => \$opt_out,

	);

die "input vcf is not provided!\n" if !$opt_input_vcf;
die "input bealge4 vcf is not provided!\n" if !$input_beagle4;
die "out is not provided!\n" if !$opt_out;

##########
##software
##########
my $bgzip = "/scratch/cgg/weiq1/software/bin/bgzip";

###########
##Make dir
###########
=a
my $opt_out_dir;
my $pid;

if ($opt_out =~ m/(.*\/)(.*)/gi)
{
	$opt_out_dir = $1;
	($pid) = split (/\./,$2);
}


#mkdir ($opt_out_dir,0755) unless (-d $opt_out_dir);
chop($opt_out_dir) if $opt_out_dir =~ /\/$/;

eval {mkpath($opt_out_dir)};
if($@)
{
      die ( "Make path [$opt_out_dir] failed:\n$@ ");
}

my $log_dir = "$opt_out_dir/../logs";

mkdir ($log_dir,0755) unless (-d $log_dir);

=cut
##########################################################################################################
##join vcf files
##########################################################################################################

###########
##get beagle4
###########

open HD, "<$input_beagle4" if $input_beagle4 !~ /\.gz$/;
open HD, "gunzip -c $input_beagle4|" if $input_beagle4 =~ /\.gz$/;

my @beagle4_header;
my @header;
my %beagle4_individual;
my %beagle4_formart;
my $i = 0;

while (my $line = <HD>)         
{       

	chomp($line);
	next if $line !~ /\S+/;

	my @f = split(/\t/,$line);

	if($f[0] =~ /##/) 
	{
		$beagle4_header[$i] = $line;
		$i++;

		next;
	}
	elsif($f[0] =~ /#/)
	{
		@header=@f;
		next;
	}


	my $chr = $f[0];
		$chr =~ s/chr//;
		#if($chr eq "X" || $chr eq "Y" ) {next;}
	my $pos = $f[1];
	my $ref = $f[3];
	my $alt = $f[4];

	for(my $k=9; $k<@f; $k++) 
	{
		$beagle4_individual{$chr}{$pos}{$ref}{$alt}{$header[$k]} = $f[$k];
		$beagle4_formart{$chr}{$pos}{$ref}{$alt} = join("\t",@f[6,7,8]);
	}

}    
close(HD);

###########
##join
###########
open HD, "<$opt_input_vcf" if $opt_input_vcf !~ /\.gz$/;
open HD, "gunzip -c $opt_input_vcf|" if $opt_input_vcf =~ /\.gz$/;

open OUT, "|$bgzip -c >$opt_out" or die "";

print OUT join("\n",@beagle4_header)."\n";


while (my $line = <HD>)         
{       
	chomp($line);
	next if $line !~ /\S+/;

	my @f = split(/\t/,$line);

	if($f[0] =~ /##/) 
	{
		print OUT $line."\n";
		next;
	}
	elsif($f[0] =~ /#/)
	{
		print OUT $line."\n";
		@header=@f;
		next;
	}

	my $chr = $f[0];
		$chr =~ s/chr//;
		#if($chr eq "X" || $chr eq "Y" ) {next;}
	my $pos = $f[1];
	my $ref = $f[3];
	my $alt = $f[4];

	if (exists $beagle4_formart{$chr}{$pos}{$ref}{$alt}) 
	{
		my @f2 = split(/\t/,$beagle4_formart{$chr}{$pos}{$ref}{$alt});

		$f[6]="$f2[0];$f[6]";
		$f[7]="$f2[1];$f[7]";
		$f[8]="GT:DS:GP:$f[8]";

		print OUT join("\t",@f[0..8]);

		for(my $k=9; $k<@f; $k++) 
		{
			if (exists $beagle4_individual{$chr}{$pos}{$ref}{$alt}{$header[$k]}) 
			{
				print OUT "\t".$beagle4_individual{$chr}{$pos}{$ref}{$alt}{$header[$k]}.":$f[$k]";
			}
			else
			{
				print OUT "\t".".|.:.:.,.,.".":$f[$k]";
			}
		}

		print OUT "\n";
	}
	else
	{
		$f[6]=".;$f[6]";
		$f[7]="AR2=;DR2=;AF=;$f[7]";
		$f[8]="GT:DS:GP:$f[8]";	

		print OUT join("\t",@f[0..8]);

		for(my $k=9; $k<@f; $k++) 
		{
			print OUT "\t".".|.:.:.,.,.".":$f[$k]";
		}

		print OUT "\n";	
	}

}    
close(HD);

close(OUT);	

