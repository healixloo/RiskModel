use strict;
use warnings;

my $gtfFile="human.gtf";
my $expFile="mRNAmatrix.txt";
my $outFile="symbol.txt";

my %hash=();
open(RF,"$gtfFile") or die $!;
while(my $line=<RF>)
{
	chomp($line);
	if($line=~/gene_id \"(.+?)\"\;.+gene_name "(.+?)"\;.+gene_biotype \"(.+?)\"\;/)
	{
		$hash{$1}=$2;
	}
}
close(RF);

open(RF,"$expFile") or die $!;
open(WF,">$outFile") or die $!;
while(my $line=<RF>)
{
	if($.==1)
	{
		print WF $line;
		next;
	}
	chomp($line);
	my @arr=split(/\t/,$line);
	$arr[0]=~s/(.+)\..+/$1/g;
	if(exists $hash{$arr[0]})
	{
		$arr[0]=$hash{$arr[0]};
		print WF join("\t",@arr) . "\n";
	}
}
close(WF); 
close(RF);


