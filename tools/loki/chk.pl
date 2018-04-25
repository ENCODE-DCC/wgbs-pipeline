#!/usr/bin/perl -w
use strict;
my $fg=0;
while(<>) {
	 if(!$fg && /^#ifdef FUNC_NAME/) {
		  $fg=1;
	 } elsif($fg==1 && /^#undef FUNC_NAME/) {
		  $fg=2;
	 } elsif($fg==2 && /^#endif/) {
		  $fg=3;
	 } elsif($fg==3 && /^#define FUNC_NAME/) {
		  $fg=0;
	 } else {
		  s/FUNC_NAME/__func__/g;
		  print $_;
		  $fg=0;
	 }
}
