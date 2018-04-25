#!/usr/bin/perl -w
use strict;
while(<>) {
	 chomp;
	 my @fd=split "\t";
	 for(my $i=5;$i<=$#fd;$i++) {
		  if($fd[$i]=~/(\S)(\S)/) {
				$fd[$i]="$1 $2";
		  } else {
				$fd[$i]="$fd[$i] $fd[$i]";
		  }
	 }
	 print join("\t",@fd),"\n";
}
