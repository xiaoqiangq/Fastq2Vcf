#!/usr/bin/perl
use strict;
use warnings;

my %hash;
my %hash2;

system ("find -name summary.txt -type f >a11"); 
open(FROM,"a11") or die();
while(my $from=<FROM>)
{	
	open(FROM2,"$from") or die();
	while(my $line=<FROM2>)
	{
		chomp($line);  	
		my @hwusi=split(/\t/,$line);
		
		
		$hash{"$hwusi[2]"}{"$hwusi[1]"}="$hwusi[0]";				
		$hash2{"$hwusi[1]"}="$hwusi[1]";
	}
	close(FROM2);	
}
close(FROM);


print "ID";

foreach my $Samplename ( sort { $a cmp $b} keys %hash2) 
{
	print "\t$hash2{$Samplename}";
}
print "\n";


foreach my $Samplename ( sort { $a cmp $b} keys %hash) 
{
	print "$Samplename";
	foreach my $feature (sort { $a cmp $b} keys %{$hash{$Samplename}}) 
	{

		print "\t$hash{$Samplename}{$feature}";

	}
	print "\n";
}

