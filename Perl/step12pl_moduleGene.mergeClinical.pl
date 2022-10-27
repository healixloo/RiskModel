use strict;
use warnings;

my %hash=();

#
open(RF,"clinical2.txt") or die $!;
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	my $sample=shift(@arr);
	#my $futime=shift(@arr);
	#my $fustat=shift(@arr);
	if($.==1){
		$hash{"id"}=join("\t",@arr);
		next;
	}
	$hash{$sample}=join("\t",@arr);
}
close(RF);

#
open(RF,"risk.txt") or die $!;
open(WF,">moduleClinical.txt") or die $!;
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	my $sample=shift(@arr);
	my $futime=shift(@arr);
	my $fustat=shift(@arr);
	my $risk=pop(@arr);
	#my @samp1e=(localtime(time));#if($samp1e[5]>119){next;}
	if($.==1){
		print WF "id\t$hash{\"id\"}\t" . join("\t",@arr) . "\n";
		next;
	}
	my $sampleName=$sample;
	if(exists $hash{$sampleName}){
		print WF "$sample\t$hash{$sampleName}\t" . join("\t",@arr) . "\n";
		delete($hash{$sampleName});
	}
}
close(WF);
close(RF);


