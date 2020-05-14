#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use lib "/home/lib13/lib";
use lib "/home/weiq1/lib/perl/lib/site_perl/5.24.0";
use Text::CSV_PP;

#my $bgzip = "/scratch/cgg/weiq1/software/bin/bgzip";
my $bgzip = "/scratch/cgg/software/htslib/htslib-1.3.1/bgzip";

my $opt_vcf;
my $annovar_summary;
my $out;
my $opt_help;

GetOptions("help|h"  => \$opt_help,
	"vcf|v=s" => \$opt_vcf,
	"annovar|a=s" => \$annovar_summary,
	"out|o=s" => \$out,
	);

#########################################
## how to use
#######################################
&usage() if $opt_help;

die "opt_vcf is not specified!\n" if !$opt_vcf;
die "annovar_summary is not specified!\n" if !$annovar_summary;
die "out is not specified!\n" if !$out;


#################################
my %feature_index = ();

my $csv = Text::CSV_PP->new ({ binary => 1 }) or die "Cannot use CSV: ".Text::CSV->error_diag ();

open IN, "<$annovar_summary" or die "Open annovar_summary  file failed!\n" ;

{
	my $line=<IN>;
	chomp($line);
	my $status  = $csv->parse($line);        # parse a CSV string into fields
	my @f = $csv->fields();            # get the parsed fields

	if($f[0] eq "Chr")
	{
		for(my $i=0; $i<@f; $i++) { $feature_index{$f[$i]} = $i; }
	}
}

open IN2, "<$opt_vcf" if $opt_vcf !~ /\.gz$/;
open IN2, "gunzip -c $opt_vcf|" if $opt_vcf =~ /\.gz$/;
open OUT, "|$bgzip -c >$out" or die "";

my %format_index = ();
my @header = ();
my %pid_index = ();


while(my $line=<IN2>)
{
	chomp($line);
	my @f = split(/\t/, $line);
	#print $line."\n";

	if ($f[0] =~ /##/)
	{
		print OUT $line."\n";

		#print "familyid\tfunc\texon_func\tgene\taachange\tconserved\tsegdup\tdbSNP138\tesp_maf\tkg_maf\tavsift\tsift\tpp2_hdiv\tpp2_hvar\tltr\tmut_taster\tmut_assessor\tfathmm\tradialsvm\tlr\tgerp\tphylop\tsiphy\tclinvar\tcadd\tcosmic70\toccurence\n";
		next;
	}
	if ($f[0] =~ /#/)
	{
		print OUT "##FILTER=<ID=Multiple,Description=\"Multiple Allele\">\n";
		print OUT "##FILTER=<ID=Error,Description=\"Error something wrong\">\n";
		print OUT "##INFO=<ID=SegDup,Number=1,Type=Integer,Description=\"Segment duplication\">\n";
		print OUT "##INFO=<ID=Func,Number=1,Type=String,Description=\"Function\">\n";
		print OUT "##INFO=<ID=ExonicFunc,Number=1,Type=String,Description=\"ExonicFunc\">\n";
		print OUT "##INFO=<ID=Gene,Number=1,Type=String,Description=\"Gene name\">\n";
		print OUT "##INFO=<ID=AAChange,Number=1,Type=String,Description=\"AAChange\">\n";
		print OUT "##INFO=<ID=ESP,Number=1,Type=Integer,Description=\"ESP6500siv2_ALL\">\n";
		print OUT "##INFO=<ID=KG,Number=1,Type=Integer,Description=\"1000g2015Aug_ALL\">\n";
		print OUT "##INFO=<ID=clinvar,Number=1,Type=String,Description=\"clinvar_20200316\">\n";

		print OUT "##INFO=<ID=cosmic,Number=1,Type=String,Description=\"cosmic91\">\n";
		print OUT "##INFO=<ID=Occurence,Number=1,Type=String,Description=\"Occurence\">\n";
		####
		print OUT "##INFO=<ID=dbnsfp35a,Number=1,Type=String,Description=\"dbnsfp35a\">\n";
		print OUT "##INFO=<ID=avsnp150,Number=1,Type=String,Description=\"avsnp150\">\n";
		print OUT "##INFO=<ID=CADD_phred,Number=1,Type=String,Description=\"CADD_phred\">\n";
		print OUT "##INFO=<ID=SIFT_pred,Number=1,Type=String,Description=\"SIFT_pred\">\n";
		print OUT "##INFO=<ID=Polyphen2_HDIV_pred,Number=1,Type=String,Description=\"Polyphen2_HDIV_pred\">\n";
		print OUT "##INFO=<ID=Polyphen2_HVAR_pred,Number=1,Type=String,Description=\"Polyphen2_HVAR_pred\">\n";
		print OUT "##INFO=<ID=LRT_pred,Number=1,Type=String,Description=\"LRT_pred\">\n";
		print OUT "##INFO=<ID=MutationTaster_pred,Number=1,Type=String,Description=\"MutationTaster_pred\">\n";
		print OUT "##INFO=<ID=MutationAssessor_pred,Number=1,Type=String,Description=\"MutationAssessor_pred\">\n";
		print OUT "##INFO=<ID=FATHMM_pred,Number=1,Type=String,Description=\"FATHMM_pred\">\n";
		print OUT "##INFO=<ID=RadialSVM_pred,Number=1,Type=String,Description=\"RadialSVM_pred\">\n";
		print OUT "##INFO=<ID=LR_pred,Number=1,Type=String,Description=\"LR_pred\">\n";
		print OUT "##INFO=<ID=GERP++,Number=1,Type=String,Description=\"GERP++\">\n";
		print OUT "##INFO=<ID=PhyloP,Number=1,Type=String,Description=\"phyloP20way_mammalian\">\n";
		print OUT "##INFO=<ID=SiPhy,Number=1,Type=String,Description=\"SiPhy_29way_logOdds\">\n";
		print OUT "##INFO=<ID=gwasCatalog,Number=1,Type=String,Description=\"gwasCatalog\">\n";
		print OUT "##INFO=<ID=Deleterious,Number=1,Type=String,Description=\"Deleterious\">\n";
		print OUT "##INFO=<ID=cytoBand,Number=1,Type=String,Description=\"cytoBand\">\n";

		print OUT "##INFO=<ID=ExAC_Freq,Number=1,Type=String,Description=\"ExAC_Freq exac03nontcga\">\n";
		print OUT "##INFO=<ID=ExAC_AFR,Number=1,Type=String,Description=\"ExAC_AFR exac03nontcga\">\n";
		print OUT "##INFO=<ID=ExAC_AMR,Number=1,Type=String,Description=\"ExAC_AMR exac03nontcga\">\n";
		print OUT "##INFO=<ID=ExAC_EAS,Number=1,Type=String,Description=\"ExAC_EAS exac03nontcga\">\n";
		print OUT "##INFO=<ID=ExAC_FIN,Number=1,Type=String,Description=\"ExAC_FIN exac03nontcga\">\n";
		print OUT "##INFO=<ID=ExAC_NFE,Number=1,Type=String,Description=\"ExAC_NFE exac03nontcga\">\n";
		print OUT "##INFO=<ID=ExAC_OTH,Number=1,Type=String,Description=\"ExAC_OTH exac03nontcga\">\n";
		print OUT "##INFO=<ID=ExAC_SAS,Number=1,Type=String,Description=\"ExAC_SAS exac03nontcga\">\n";


		print OUT "##INFO=<ID=Eigen_coding_or_noncoding,Number=1,Type=String,Description=\"Eigen_coding_or_noncoding\">\n";
		print OUT "##INFO=<ID=Eigen_raw,Number=1,Type=String,Description=\"Eigen_raw\">\n";
		print OUT "##INFO=<ID=phyloP100way_vertebrate,Number=1,Type=String,Description=\"phyloP100way_vertebrate\">\n";
		print OUT "##INFO=<ID=phastCons100way_vertebrate,Number=1,Type=String,Description=\"phastCons100way_vertebrate\">\n";
		print OUT "##INFO=<ID=phastCons20way_mammalian,Number=1,Type=String,Description=\"phastCons20way_mammalian\">\n";
		print OUT "##INFO=<ID=Interpro_domain,Number=1,Type=String,Description=\"Interpro_domain\">\n";
		print OUT "##INFO=<ID=GTEx_V6p_gene,Number=1,Type=String,Description=\"GTEx_V6p_gene\">\n";
		print OUT "##INFO=<ID=GTEx_V6p_tissue,Number=1,Type=String,Description=\"GTEx_V6p_tissue\">\n";
				
		#print OUT "##INFO=<ID=MAF,Number=1,Type=String,Description=\"MAF\">\n";
		#print OUT "##INFO=<ID=ANNO,Number=1,Type=String,Description=\"ANNO_for_EPACTS\">\n";

		print OUT $line."\n";

		#print "familyid\tfunc\texon_func\tgene\taachange\tconserved\tsegdup\tdbSNP138\tesp_maf\tkg_maf\tavsift\tsift\tpp2_hdiv\tpp2_hvar\tltr\tmut_taster\tmut_assessor\tfathmm\tradialsvm\tlr\tgerp\tphylop\tsiphy\tclinvar\tcadd\tcosmic70\toccurence\n";
		next;
	}

##############
##annovar
############
	my $line_annovar=<IN>;

	chomp($line_annovar);
	my $status  = $csv->parse($line_annovar);        # parse a CSV string into fields
	my @f_annovar = $csv->fields();            # get the parsed fields	

	my $chr_annovar = $f_annovar[$feature_index{"Chr"}];
	my $start_annovar = $f_annovar[$feature_index{"Start"}];
	my $end_annovar = $f_annovar[$feature_index{"End"}];
	my $ref_annovar = $f_annovar[$feature_index{"Ref"}];
	my $alt_annovar = $f_annovar[$feature_index{"Alt"}];

	my $func = $f_annovar[$feature_index{"Func"}];
	my $gene = $f_annovar[$feature_index{"Gene"}];
	my $exon_func = $f_annovar[$feature_index{"ExonicFunc"}];

	my $aachange = $f_annovar[$feature_index{"AAChange"}];
	#my $conserved = $f_annovar[$feature_index{'Conserved'}];
	my $segdup = $f_annovar[$feature_index{'SegDup'}];
	my $esp_maf = $f_annovar[$feature_index{"ESP_ALL"}];
	my $kg_maf = $f_annovar[$feature_index{"KG_ALL"}];	
	my $dbsnp138 = $f_annovar[$feature_index{"dbSNP"}];
	#my $avsift = $f_annovar[$feature_index{'AVSIFT'}];
	
	my $sift = $f_annovar[$feature_index{"SIFT_pred"}];
	my $pp2_hdiv = $f_annovar[$feature_index{"Polyphen2_HDIV_pred"}];
	my $pp2_hvar = $f_annovar[$feature_index{"Polyphen2_HVAR_pred"}];
	my $ltr = $f_annovar[$feature_index{"LRT_pred"}];
	my $mut_taster = $f_annovar[$feature_index{"MutationTaster_pred"}];
	my $mut_assessor = $f_annovar[$feature_index{"MutationAssessor_pred"}];
	my $fathmm = $f_annovar[$feature_index{"FATHMM_pred"}];
	my $provean = $f_annovar[$feature_index{"PROVEAN_pred"}];		
	my $radialsvm = $f_annovar[$feature_index{"MetaSVM_pred"}];
	my $lr = $f_annovar[$feature_index{"MetaLR_pred"}];
	my $gerp = $f_annovar[$feature_index{'GERP++_RS'}];
	my $phylop = $f_annovar[$feature_index{"phyloP20way_mammalian"}];		
	my $siphy = $f_annovar[$feature_index{"SiPhy_29way_logOdds"}];
		
	my $clinvar = $f_annovar[$feature_index{'clinvar'}];
	#my $cadd = $f_annovar[$feature_index{'Caddgt10_Pred'}];
	my $cadd = $f_annovar[$feature_index{'CADD_phred'}];
	my $cosmic70 = $f_annovar[$feature_index{'cosmic'}];
	my $occurence = $f_annovar[$feature_index{'Occurence'}];
###
	my $exac_freq = $f_annovar[$feature_index{'ExAC_Freq'}];
	my $exac_afr = $f_annovar[$feature_index{'ExAC_AFR'}];
	my $exac_amr = $f_annovar[$feature_index{'ExAC_AMR'}];
	my $exac_eas = $f_annovar[$feature_index{'ExAC_EAS'}];
	my $exac_fin = $f_annovar[$feature_index{'ExAC_FIN'}];
	my $exac_nfe = $f_annovar[$feature_index{'ExAC_NFE'}];
	my $exac_oth = $f_annovar[$feature_index{'ExAC_OTH'}];
	my $exac_sas = $f_annovar[$feature_index{'ExAC_SAS'}];


	my $gwas = $f_annovar[$feature_index{"gwasCatalog"}];
	my $cytoBand = $f_annovar[$feature_index{'cytoBand'}];

	my $Eigen_coding_or_noncoding = $f_annovar[$feature_index{'Eigen_coding_or_noncoding'}];
	my $Eigen_raw = $f_annovar[$feature_index{'Eigen-raw'}];

	my $phyloP100way_vertebrate = $f_annovar[$feature_index{'phyloP100way_vertebrate'}];
	my $phastCons100way_vertebrate = $f_annovar[$feature_index{'phastCons100way_vertebrate'}];	
	my $phastCons20way_mammalian = $f_annovar[$feature_index{'phastCons20way_mammalian'}];

	my $Interpro_domain = $f_annovar[$feature_index{'Interpro_domain'}];
	my $GTEx_V6p_gene = $f_annovar[$feature_index{'GTEx_V6p_gene'}];
	my $GTEx_V6p_tissue = $f_annovar[$feature_index{'GTEx_V6p_tissue'}];


	my @fields = ();

	push(@fields, $chr_annovar); #0
	push(@fields, $start_annovar); #1
	push(@fields, $ref_annovar); #2
	push(@fields, $alt_annovar); #3

	push(@fields, $func); #4
	push(@fields, $exon_func); #5
	push(@fields, $gene); #6

	push(@fields, $aachange); #7
	push(@fields, $GTEx_V6p_tissue);#8
	push(@fields, $segdup);#9
	push(@fields, $dbsnp138); #10
	push(@fields, $esp_maf); #11
	push(@fields, $kg_maf); #12
	#push(@fields, $avsift); #13

	push(@fields, $sift); #13
	push(@fields, $pp2_hdiv); #14
	push(@fields, $pp2_hvar); #15
	push(@fields, $ltr); #16
	push(@fields, $mut_taster); #17
	push(@fields, $mut_assessor); #18
	push(@fields, $fathmm); #19
	push(@fields, $radialsvm); #20
	push(@fields, $lr);#21
	push(@fields, $gerp); #22
	push(@fields, $phylop); #23
	push(@fields, $siphy); #24

	push(@fields, $clinvar); #25
	push(@fields, $cadd); #26
	push(@fields, $cosmic70); #27
	push(@fields, $occurence); #28
	push(@fields, $gwas); #29	
	push(@fields, $cytoBand); #30	

	push(@fields, $exac_freq); #31
	push(@fields, $exac_afr); #32
	push(@fields, $exac_amr); #33
	push(@fields, $exac_eas); #34
	push(@fields, $exac_fin); #35
	push(@fields, $exac_nfe); #36	
	push(@fields, $exac_oth); #37
	push(@fields, $exac_sas); #38

	
	push(@fields, $Eigen_coding_or_noncoding); #39
	push(@fields, $Eigen_raw); #40
	push(@fields, $phyloP100way_vertebrate); #41
	push(@fields, $phastCons100way_vertebrate); #42
	push(@fields, $phastCons20way_mammalian); #43
	push(@fields, $Interpro_domain); #44	
	push(@fields, $GTEx_V6p_gene); #45


	for(my $i=0; $i<@fields; $i++)
	{
		if ( !$fields[$i] || $fields[$i] eq ".")
		{
			$fields[$i] = "";
		}
		elsif($fields[$i] =~ /\;/)
		{
			$fields[$i] =~ s/;/,/g;
		}

		$fields[$i] =~ s/\s+/_/g;	
	}

################################


	my $chr = $f[0];
	#$chr =~ s/chr//;
	
	my $pos = $f[1];
	my $ref = $f[3];
	my $alt = $f[4];


# my @var_anno_info = @{$anno_info{$chr}{$pos}{$ref}{$alt}};

	if ($chr eq "Y") {}
	elsif($alt !~/,/)
	{
		if(length($ref)==1 && length($alt)==1)
		{
			if( $chr_annovar eq $chr && $start_annovar eq $pos && $ref_annovar eq $ref && $alt_annovar eq $alt)
			{}
			else
			{
				if ($f[6] ne ".") 
				{
					$f[6] .= ";Error"; #not match ref and alt\n
				}
			}
		}
		elsif(length($alt)>1)
		{
			if($chr_annovar eq $chr && $start_annovar eq $pos)
			{}
			else
			{
				if ($f[6] ne ".") 
				{
					$f[6] .= ";Error";
				}
			} #insertion position is not match \n
		}
		elsif(length($ref)>1)
		{
			my $pos1 = $pos+1;
			if($chr_annovar eq $chr && $start_annovar eq $pos1)
			{}
			else 
			{
				if ($f[6] ne ".") 
				{
					$f[6] .= ";Error";
				}
			}
		}
	}
	elsif($alt =~/,/)
	{		
		if ($f[6] ne ".") 
		{
			$f[6] .= ";Multiple";
		}

		my @tmp =split("\,",$alt);

		my $alt2;

		if ($tmp[0] eq "*") {
			$alt2 = $tmp[1];
		}else
		{
			$alt2 = $tmp[0];
		}
		

		if(length($ref)==1 && length($alt2)==1)
		{
			if($chr_annovar eq $chr && $start_annovar eq $pos && $ref_annovar eq $ref && $alt_annovar eq $alt2)
			{}
			else
			{
				if ($f[6] ne ".") 
				{
					$f[6] .= ";Error";
				}
			}
		}elsif(length($ref)>1 && length($alt2)>1)
		{
			if($chr_annovar eq $chr && ($start_annovar - $pos < 50 || $start_annovar - $pos > -50 ))
			{}
			else
			{
				if ($f[6] ne ".") 
				{
					$f[6] .= ";Error";
				}
			}

		}
		elsif(length($alt2)>1)
		{
			if($chr_annovar eq $chr && $start_annovar eq $pos)
			{}
			else
			{
				if ($f[6] ne ".") 
				{
					$f[6] .= ";Error";
				}
			} #insertion position is not match \n
		}
		elsif(length($ref)>1)
		{
			my $pos1 = $pos+1;
			if($chr_annovar eq $chr && $start_annovar eq $pos1)
			{}
			else
			{
				if ($f[6] ne ".") 
				{
					$f[6] .= ";Error";
				}
			}
		}

	}else
	{
		if ($f[6] ne ".") 
		{
			$f[6] .= ";Error";
		}
	}

	my @var_anno_info = @fields;

	if (@var_anno_info) 
	{
	##delete none annotation

	####filter
		my $D_cnt_D = 0;
		my $D_cnt_T = 0;

		my $D_cnt_P = 0;
		my $D_cnt_null = 0;
		
		#sift
		if ($var_anno_info[13] eq "D"){ $D_cnt_D++;} elsif($var_anno_info[13] eq "T"){ $D_cnt_T++;} else {$D_cnt_null++;}
		
		#pp2_hdiv
		if ($var_anno_info[14] eq "D"){ $D_cnt_D++;} elsif($var_anno_info[14] eq "P"){ $D_cnt_P++;} elsif($var_anno_info[14] eq "B"){$D_cnt_T++;} else {$D_cnt_null++;}
		
		#pp2_hdvar	
		if ($var_anno_info[15] eq "D"){ $D_cnt_D++;} elsif($var_anno_info[15] eq "P"){ $D_cnt_P++;} elsif($var_anno_info[15] eq "B"){$D_cnt_T++;} else {$D_cnt_null++;}
		
		#ltr
		if ($var_anno_info[16] eq "D"){ $D_cnt_D++;} elsif($var_anno_info[16] eq "N"){ $D_cnt_T++;} elsif($var_anno_info[16] eq "U"){$D_cnt_null++;} else {$D_cnt_null++;}

		#mut_taster
		if ($var_anno_info[17] eq "D"){ $D_cnt_D++;} elsif($var_anno_info[17] eq "A"){ $D_cnt_P++;} elsif($var_anno_info[17] eq "N" || $var_anno_info[17] eq "P"){$D_cnt_T++;} else {$D_cnt_null++;}

		#mut_assessor
		if ($var_anno_info[18] eq "H"){ $D_cnt_D++;} elsif($var_anno_info[18] eq "M"){ $D_cnt_P++;} elsif($var_anno_info[18] eq "L" || $var_anno_info[18] eq "N"){$D_cnt_T++;} else {$D_cnt_null++;}

		#fathmm
		if ($var_anno_info[19] eq "D"){ $D_cnt_D++;} elsif($var_anno_info[19] eq "T"){ $D_cnt_T++;} else {$D_cnt_null++;}

		#radialsvm
		if ($var_anno_info[20] eq "D"){ $D_cnt_D++;} elsif($var_anno_info[20] eq "T"){ $D_cnt_T++;} else {$D_cnt_null++;}

		#lr
		if ($var_anno_info[21] eq "D"){ $D_cnt_D++;} elsif($var_anno_info[21] eq "T"){ $D_cnt_T++;} else {$D_cnt_null++;}

		#$D_cnt3 ++ if ($var_anno_info[27] ne "NA" && $var_anno_info[27] >20); #cadd

		#my @familyid = split(/:/,$f[6]);


		$f[7] = join(";","SegDup=$var_anno_info[9]","Func=$var_anno_info[4]","ExonicFunc=$var_anno_info[5]","Gene=$var_anno_info[6]","AAChange=$var_anno_info[7]","ESP=$var_anno_info[11]","KG=$var_anno_info[12]","clinvar=$var_anno_info[25]","cosmic=$var_anno_info[27]","Occurence=$var_anno_info[28]","CADD_phred=$var_anno_info[26]").";".join(";","SIFT_pred=$var_anno_info[13]","Polyphen2_HDIV_pred=$var_anno_info[14]","Polyphen2_HVAR_pred=$var_anno_info[15]","LRT_pred=$var_anno_info[16]","MutationTaster_pred=$var_anno_info[17]","MutationAssessor_pred=$var_anno_info[18]","FATHMM_pred=$var_anno_info[19]","RadialSVM_pred=$var_anno_info[20]","LR_pred=$var_anno_info[21]","GERP++=$var_anno_info[22]","PhyloP=$var_anno_info[23]","SiPhy=$var_anno_info[24]").";Deleterious=D,$D_cnt_D:P,$D_cnt_P:T,$D_cnt_T:null,$D_cnt_null;gwasCatalog=$var_anno_info[29];cytoBand=$var_anno_info[30];ExAC_Freq=$var_anno_info[31];ExAC_AFR=$var_anno_info[32];ExAC_AMR=$var_anno_info[33];ExAC_EAS=$var_anno_info[34];ExAC_FIN=$var_anno_info[35];ExAC_NFE=$var_anno_info[36];ExAC_OTH=$var_anno_info[37];ExAC_SAS=$var_anno_info[38];Eigen_coding_or_noncoding=$var_anno_info[39];Eigen_raw=$var_anno_info[40];phyloP100way_vertebrate=$var_anno_info[41];phastCons100way_vertebrate=$var_anno_info[42];phastCons20way_mammalian=$var_anno_info[43];Interpro_domain=$var_anno_info[44];GTEx_V6p_gene=$var_anno_info[45];GTEx_V6p_tissue=$var_anno_info[8];$f[7]";		
		

		if($var_anno_info[10]){$f[2] = $var_anno_info[10];}
		else{$f[2] = ".";}
	}

	print OUT join("\t",@f)."\n";
	
}

close IN;
close IN2;
close OUT;

#####################
##
#################
sub usage()
{
		
	print "\n*******************************************************\n";	
	print "*************************\n";

	print "\n*** The script was complied on Feb 10 2018 15:21:20 ***\n";

	print "\t\tExample 1.1: *.pl -vcf *.vcf.gz -annovar *.vcf.genome_summary.csv -o *.annovar.vcf.gz \n";

	exit(0);
}
