
use strict;
use warnings;

my %hash=();

#
open(RF,"clinical.txt") or die $!;
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	my $sample=shift(@arr);
	if($.==1){
		$hash{"id"}=join("\t",@arr);
		next;
	}
	$hash{$sample}=join("\t",@arr);
}
close(RF);

#
open(RF,"risk.txt") or die $!;
open(WF,">indepInput.txt") or die $!;
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	my $sample=shift(@arr);
	my $risk=pop(@arr);
	my $riskScore=pop(@arr);
	# my @samp1e=(localtime(time));if($samp1e[5]>119){next;}
        if($risk=~/low/){
           $risk=0}elsif($risk=~/high/){
           $risk=1
        }
	if($.==1){
		print WF "id\t$hash{\"id\"}\t" . "risk\n";
		next;
	}
	my $sampleName=$sample;
	if(exists $hash{$sampleName}){
		print WF "$sample\t$hash{$sampleName}\t" . $risk . "\n";
		delete($hash{$sampleName});
	}
}
close(WF);
close(RF);


