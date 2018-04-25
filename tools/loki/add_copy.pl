#!/usr/bin/perl
use strict;
my $flag;
while(<>) {
	 $flag=1 if(/Copyright/);
	 if(!$flag && m%[*]{20,}/$%) {
		  print " * Copyright (C) Simon C. Heath 1997, 2000, 2002                            *\n";
		  print " * This is free software.  You can distribute it and/or modify it           *\n";
		  print " * under the terms of the Modified BSD license, see the file COPYING        *\n";
		  print " *                                                                          *\n";
	 }
	 print $_;
}
