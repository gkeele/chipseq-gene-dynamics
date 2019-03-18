#! /usr/bin/perl -w
use strict;

my $LtD_setd2 = shift;
my $LtD_setd2_rph1d = shift;
my $DtL_setd2 = shift;
my $nonOverlapping_genes = shift;
my $out = shift;

my %nonOverlapping_genes_hash;
open(GENE, "<$nonOverlapping_genes");
while(<GENE>){
	chomp($_);
	
	my @linedata = split(/\t/, $_);
	$nonOverlapping_genes_hash{$linedata[3]} = 1;
}
close(GENE);

open(LTD, "<$LtD_setd2");
open(LTD_RPH1D, "<$LtD_setd2_rph1d");
open(DTL, "<$DtL_setd2");
open(OUT, ">$out");

my $header = <LTD>;
my $header1 = <LTD_RPH1D>;
my $header2 = <DTL>;

while(my $LtD_setd2_line = <LTD>){
	my $LtD_setd2_rph1d_line = <LTD_RPH1D>;
	my $DtL_setd2_line = <DTL>;
	
	chomp($LtD_setd2_line);
	chomp($LtD_setd2_rph1d_line);
	chomp($DtL_setd2_line);
	
	my @LtD_setd2_linedata = split(/\t/, $LtD_setd2_line);
	my @LtD_setd2_rph1d_linedata = split(/\t/, $LtD_setd2_rph1d_line);
	my @DtL_setd2_linedata = split(/\t/, $DtL_setd2_line);
	
	if(($LtD_setd2_linedata[0] ne $LtD_setd2_rph1d_linedata[0]) or ($LtD_setd2_linedata[0] ne $DtL_setd2_linedata[0])){
		die "Inconsistent gene name at line $..\n";
	}
	elsif(!defined($nonOverlapping_genes_hash{$LtD_setd2_linedata[0]})){
		next;
	}
	
	#geneName, assay (1=LtD, 2=LtD_setd2_rph1d, 3=DtL), Rep#, time, value
	print OUT "$LtD_setd2_linedata[0]\t1\t1\t0\t$LtD_setd2_linedata[1]\n";
	print OUT "$LtD_setd2_linedata[0]\t1\t1\t1\t$LtD_setd2_linedata[2]\n";
	print OUT "$LtD_setd2_linedata[0]\t1\t1\t2\t$LtD_setd2_linedata[3]\n";
	print OUT "$LtD_setd2_linedata[0]\t1\t1\t3\t$LtD_setd2_linedata[4]\n";
	print OUT "$LtD_setd2_linedata[0]\t1\t2\t0\t$LtD_setd2_linedata[5]\n";
	print OUT "$LtD_setd2_linedata[0]\t1\t2\t1\t$LtD_setd2_linedata[6]\n";
	print OUT "$LtD_setd2_linedata[0]\t1\t2\t2\t$LtD_setd2_linedata[7]\n";
	print OUT "$LtD_setd2_linedata[0]\t1\t2\t3\t$LtD_setd2_linedata[8]\n";
	print OUT "$LtD_setd2_linedata[0]\t1\t3\t0\t$LtD_setd2_linedata[9]\n";
	print OUT "$LtD_setd2_linedata[0]\t1\t3\t1\t$LtD_setd2_linedata[10]\n";
	print OUT "$LtD_setd2_linedata[0]\t1\t3\t2\t$LtD_setd2_linedata[11]\n";
	print OUT "$LtD_setd2_linedata[0]\t1\t3\t3\t$LtD_setd2_linedata[12]\n";
	print OUT "$LtD_setd2_rph1d_linedata[0]\t2\t1\t0\t$LtD_setd2_rph1d_linedata[1]\n";
	print OUT "$LtD_setd2_rph1d_linedata[0]\t2\t1\t1\t$LtD_setd2_rph1d_linedata[2]\n";
	print OUT "$LtD_setd2_rph1d_linedata[0]\t2\t1\t2\t$LtD_setd2_rph1d_linedata[3]\n";
	print OUT "$LtD_setd2_rph1d_linedata[0]\t2\t1\t3\t$LtD_setd2_rph1d_linedata[4]\n";
	print OUT "$LtD_setd2_rph1d_linedata[0]\t2\t2\t0\t$LtD_setd2_rph1d_linedata[5]\n";
	print OUT "$LtD_setd2_rph1d_linedata[0]\t2\t2\t1\t$LtD_setd2_rph1d_linedata[6]\n";
	print OUT "$LtD_setd2_rph1d_linedata[0]\t2\t2\t2\t$LtD_setd2_rph1d_linedata[7]\n";
	print OUT "$LtD_setd2_rph1d_linedata[0]\t2\t2\t3\t$LtD_setd2_rph1d_linedata[8]\n";
	print OUT "$LtD_setd2_rph1d_linedata[0]\t2\t3\t0\t$LtD_setd2_rph1d_linedata[9]\n";
	print OUT "$LtD_setd2_rph1d_linedata[0]\t2\t3\t1\t$LtD_setd2_rph1d_linedata[10]\n";
	print OUT "$LtD_setd2_rph1d_linedata[0]\t2\t3\t2\t$LtD_setd2_rph1d_linedata[11]\n";
	print OUT "$LtD_setd2_rph1d_linedata[0]\t2\t3\t3\t$LtD_setd2_rph1d_linedata[12]\n";
	print OUT "$DtL_setd2_linedata[0]\t3\t1\t0\t$DtL_setd2_linedata[1]\n";
	print OUT "$DtL_setd2_linedata[0]\t3\t1\t1\t$DtL_setd2_linedata[2]\n";
	print OUT "$DtL_setd2_linedata[0]\t3\t1\t2\t$DtL_setd2_linedata[3]\n";
	print OUT "$DtL_setd2_linedata[0]\t3\t1\t3\t$DtL_setd2_linedata[4]\n";
	print OUT "$DtL_setd2_linedata[0]\t3\t2\t0\t$DtL_setd2_linedata[5]\n";
	print OUT "$DtL_setd2_linedata[0]\t3\t2\t1\t$DtL_setd2_linedata[6]\n";
	print OUT "$DtL_setd2_linedata[0]\t3\t2\t2\t$DtL_setd2_linedata[7]\n";
	print OUT "$DtL_setd2_linedata[0]\t3\t2\t3\t$DtL_setd2_linedata[8]\n";
	print OUT "$DtL_setd2_linedata[0]\t3\t3\t0\t$DtL_setd2_linedata[9]\n";
	print OUT "$DtL_setd2_linedata[0]\t3\t3\t1\t$DtL_setd2_linedata[10]\n";
	print OUT "$DtL_setd2_linedata[0]\t3\t3\t2\t$DtL_setd2_linedata[11]\n";
	print OUT "$DtL_setd2_linedata[0]\t3\t3\t3\t$DtL_setd2_linedata[12]\n";
}
close(OUT);
close(DTL);
close(LTD);
close(LTD_RPH1D);