#!/usr/bin/perl 
use warnings;
use strict;
my @files;
$files[0]=shift or die "No input file 1\n";
$files[1]=shift or die "No input file 2\n";
my %ibd;
my %flag;
for my $i(0..1) {
    my $file=$files[$i];
    open FILE,$file or die "Couldn't open file $file: $!\n";
    my $pos;
    while(<FILE>) {
	if(/Position = (\S+)/) {
	    $pos=$1;
	    $flag{$pos}=1;
	    next;
	}
	next if(/^[*]{2}/);
	my @fd=split;
	if($#fd==3) {
	    $ibd{"$fd[0]\t$fd[1]"}{$pos}[$i]=$fd[2];
	} elsif($#fd==4) {
	    $ibd{"$fd[0]\t$fd[1]\t$fd[2]"}{$pos}[$i]=$fd[3];
	}
    }
}
for my $id(keys %ibd) {
    print "$id";
    for my $ps(sort {$a<=>$b} keys %flag) {
	my $t=$ibd{$id}{$ps};
	if(defined $t) {
	    print "\t",$$t[0] || 0;
	    print "\t",$$t[1] || 0;
	} else {
	    print "\t0\t0"
	}
    }
    print "\n";
}
