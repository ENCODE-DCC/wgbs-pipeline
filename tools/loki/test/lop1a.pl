#!/usr/bin/perl
#
# Script to estimate the LOP1 statistic for linkage
# for small regions of chromosome and output in
# tabular format suitable for gnuplot or other
# plotting programs.  If gnuplot is installed,
# can invoke it directly to produce screen or
# postscript plots.
#
# Usage: dist.pl -q -C -f outfile -p posfile -b binsize
#                -o psfile -d datafile -c chromosome
#                -i [start iteration][,:-][stop iteration]
# outfile defaults to loki.out
# posfile defaults to loki.pos
# binsize defaults to 1[cM]
# if -c option not set, looks at all chromsomes (linkage groups) fitted
# -d sets basename for temporary datafiles
# -o psfile instructs the script to get gnuplot to produce a postscript output into psfile
# -C turns on colour postscript (only has effect with -o option)
# -q does not invoke gnuplot
# -i sets iteration range
#
# Works with loki_2.4
#
# Simon Heath - September 2000
#
use Getopt::Std;
use IO::File;
use POSIX qw(tmpnam);
use strict;
my (%opt,$model,$fg,$ch,$nmk,$lg,$file,$sex_map,$start,$stop,$tmp);
my (%chrom,@mkchrom,@mkpos,@mkname,@map_length,@map_start,@map_end,%lname,%mirror);
getopts('qCc:b:o:i:d:p:f:h?',\%opt);
if($opt{h} || $opt{'?'}) {
	 print "usage: dist.pl -q -C -f outfile -p posfile -b binsize\n               -o psfile -d datafile -c chromsome\n               -i [start iteration][,:-][stop iteration]\n";
	 exit(0);
}
if(defined($opt{c})) { die "Invalid chromosome number\n" if($opt{c}<1); }
# -i option sets a range of iterations to consider
if($opt{i}) {
	 $tmp=$opt{i};
	 if($tmp=~/^([1-9][0-9]*)(.*)/) {
		  $start=$1;
		  $tmp=$2;
	 }
	 if($tmp=~/^[,-:](.*)/) {
		  $tmp=$1;
		  if($tmp=~/^([1-9][0-9]*)/) {
				$stop=$1;
		  }
	 }
	 die "Bad -i option\n" if(!defined($start) && !defined($stop));
	 if(defined($start) && defined($stop) && $start>$stop) {
		  $tmp=$start;
		  $start=$stop;
		  $stop=$tmp;
	 }
	 $start=1 if(!defined($start));
}
$file=$opt{f}?$opt{f}:"loki.out";
$sex_map=0;
open FILE,$file or die "Could not open file '$file'\n";
#
# We need to get some details about linkage groups, chromosome positions
# and map lengths out of the header
#
while(<FILE>) {
	 # Read in model
	 $model=$1 if(/^Model: (.+)/);
	 # Get rid of qutoes as they mess up Gnuplot
	 $model=~s/\'//g;
	 if(/^Linkage groups:$/) {
		  $fg=1;
		  next;
	 }
	 if($fg==1) {
		  # Check for linkage group names and map lengths
		  if(/(\d+): (.*)$/) {
				$lg=$1+0;
				if(!$opt{c} || $lg eq $opt{c}) {
					 if($2=~/^(.*\S)\s+Map range: (.*)/) {
						  my $ch=$1;
						  my $map=$2;
						  if($ch=~/(.*)\'$/) {
								if(defined $lname{$1}) {
									 my $lg1=$lname{$1};
									 $mirror{$lg1}=$lg;
									 $ch=$1;
								} else {
									 die "Can't find real chromosome for $1\n";
								}
						  } else {
								$lname{$ch}=$lg;
						  }
						  $chrom{$lg}=$ch;
						  if($map=~/\((-?[0-9.]+)cM to (-?[0-9.]+)cM\)(.*)/) {
								$map_start[$lg][0]=$1;
								$map_end[$lg][0]=$2;
								if($3=~/\((-?[0-9.]+)cM to (-?[0-9.]+)cM\)/) {
									 $map_start[$lg][1]=$1;
									 $map_end[$lg][1]=$2;
									 $sex_map=1;
								}
						  } else { die "Error reading linkage group map range\n"; }
					 } else { die "Error reading linkage group name\n"; }
				}
				# Pick up marker names and positions
		  } elsif(/^\s+(.+)\s+-\s+([\d\.]+)\s*([\d\.]*)/) {
				if(!$opt{c} || $lg eq $opt{c}) {
					 $mkchrom[$nmk]=$lg;
					 $mkpos[$nmk][0]=$2;
					 $mkpos[$nmk][1]=$3;
					 $mkname[$nmk++]=$1;
				}
		  } else {
				# Check for total map length
				if(/^Total Map Length: (.*)$/) {
					 if($1=~/(\d\d*\.?\d*)cM(.*)/) {
						  $map_length[0]=$1;
						  if($2=~/(\d\d*\.?\d*)cM$/) {
								$map_length[1]=$1;
								$sex_map=1;
						  }
					 } else { die "Error reading map lengths\n"; }	
					 $fg=2;
					 last;
				}
		  }
	 }
}
close FILE;
die "No linked chromosomes found\n" unless $fg==2;
# Now start reading data from loki.pos
$file=$opt{p}?$opt{p}:"loki.pos";
open FILE,$file or die "Could not open file '$file'\n";
my($it,$bw,$rep,$rep1,$i,$j,$k,$x,$kk,$p,$n_qtl,$n);
my(@bin,@bin1,@sex,$sx,$skip);
$sex[0]="Male";
$sex[1]="Female";
# Set bin width
$bw=$opt{b}?$opt{b}:1;
while(<FILE>) {
	 split;
	 $k=@_;
	 # Skip empty lines
	 next if(!$k);
	 # Get number of iterations from end of line
	 if($_[$k-1]=~/(.*):(\d+)/) {
		  $rep=$2;
		  $skip=0;
		  $rep1=$rep;
		  if($opt{i}) {
				if($start>$it+$rep) {$skip=1;} 
				elsif($start>$it+1) {$rep1=$it+$rep+1-$start;}
				if(defined($stop)) {
					 if($stop<=$it) {$skip=2;}
					 elsif($stop<=$it+$rep) {$rep1-=$it+$rep-$stop;}
				}
		  }
		  if($rep1 && !$skip) {
				$_[$k-1]=$1;
				$i=0;
				$n_qtl=0;
				undef @bin1;
				# Go though rest of line
				while($i<$k) {
					 $lg=$_[$i++]; # Linkage group
					 $n=$_[$i++];  # Number if QTLs in linkage group
					 $n_qtl+=$n;   # Add on to total QTL count
					 if($lg) {     # If linked, go thought QTL positions
						  $j=0;
						  while($j<$n) {
								for($sx=0;$sx<=$sex_map;$sx++) {
									 $x=$_[$i++]; # Get position
									 $kk=int(($x-$map_start[$lg][$sx])/$bw); # Get offset from start of map, and compute bin
									 $bin1[$sx][$lg]{$kk}=1; # Add count to bin
								}
								$j++;
								die "Invalid number of columns\n" if($i>$k);
						  }
					 }
				}
				for($sx=0;$sx<=$sex_map;$sx++) {
					 foreach $lg (keys %chrom) {
						  $j=$bin1[$sx][$lg];
						  foreach $kk(keys %$j) {
								$bin[$sx][$lg][$kk]+=$rep1;
						  }
					 }
				}
		  }
		  $it+=$rep;
		  last if($skip==2);
	 }
}
close FILE;
if($opt{i}) {
	 $stop=$it if(!defined($stop));
	 $it=$stop-$start+1;
}
if($it) {
	 # Set up temporary files for writing out data
	 my(@name,$name1,@fh,$fh1,$ff,$nch);
	 foreach $lg(keys %chrom) {$nch++;}
	 if($opt{d}) {
		  foreach $lg(keys %chrom) {
				$name[$lg]=$opt{d}.".dat_".$chrom{$lg};
				$name[$lg].="_m" if($sex_map);
				$fh[$lg]=IO::File->new("> $name[$lg]") or die "Couldn't open data file\n";
				if($sex_map) {
					 $name[$lg+$nch]=$opt{d}.".dat_".$chrom{$lg}."_f";
					 $fh[$lg+$nch]=IO::File->new("> $name[$lg+$nch]") or die "Couldn't open data file\n";
				}
		  }
	 } else {
		  foreach $lg(keys %chrom) {
				for($sx=0;$sx<=$sex_map;$sx++) {
					 do{$name[$lg+$sx*$nch]=tmpnam()}
					 until $fh[$lg+$sx*$nch]=IO::File->new($name[$lg+$sx*$nch],O_RDWR|O_CREAT|O_EXCL);
				}
		  }
	 }
	 # Set handler to delete files after the script ends if the -d option was not used
	 END {
		  if(!$opt{d}) {
				foreach $lg(keys %chrom) {
					 for($i=0;$i<=$sex_map;$i++) { unlink($name[$lg+$i*$nch]); }
				}
		  }
	 }
	 # Write out data
	 my $log_10=log(10);
	 foreach $lg(keys %mirror) {
		  my $lg1=$mirror{$lg};
		  for($sx=0;$sx<=$sex_map;$sx++) {
				$x=$map_start[$lg][$sx]+.5*$bw;
				$ff=$fh[$lg+$sx*$nch];
				$kk=0;
				while($x<=$map_end[$lg][$sx]) {
					 my $n=$bin[$sx][$lg1][$kk];
					 my $z=$bin[$sx][$lg][$kk++];
					 if(!$n || !$z) {
						  $n++;
						  $z++;
					 }
					 $z/=$n;
					 print $ff $x+.5," ",log($z)/$log_10," $lg $lg1 $z\n";
					 $x+=$bw;
				}
				close $ff;
		  }
	 }
	 # Open gnuplot control file
	 if($opt{d}) {
		  $name1=$opt{d};
		  $fh1=IO::File->new("> $name1") or die "Couldn't open data file\n";
	 } else {
		  do{$name1=tmpnam()}
		  until $fh1=IO::File->new($name1,O_RDWR|O_CREAT|O_EXCL);
	 }
	 # Set handler to delete file after the script ends if the -d option was not used
	 END {unlink($name1) if(!$opt{d});}
	 # If outputting postscript, set options
	 if($opt{o}) {
		  $j="\"Times-Roman\" 14";
		  $j="color solid ".$j if($opt{C});
		  print $fh1 "set term postscript $j\nset output \'$opt{o}\'\n";
		  print $fh1 "set xlabel \'Position\' 0,-2.5\n";
	 }
	 # Write rest of gnuplot control file
	 foreach $lg(sort{$a<=>$b} keys %mirror) {
		  # Set axis labels
		  print $fh1 "set ylabel \'L-Score'\nset xlabel \'Position\'\n";
		  for($sx=0;$sx<=$sex_map;$sx++) {
				# Set plot title
				if($sex_map) { print $fh1 "set title \'$model -- $chrom{$lg} - $sex[$sx] map\'\n"; }
				else { print $fh1 "set title \'$model -- $chrom{$lg}\'\n"; }
				# Turn off keys, set grid option
				print $fh1 "set nokey\nset grid\n";
				# Set tic marks to correspond to marker positions
				$n=0;
				for($i=0;$i<$nmk;$i++) {
					 $ch=$mkchrom[$i];
					 if($ch eq $lg) {
						  $x=$mkpos[$i][$sx];
						  if(!$n) { print $fh1 "set xtics rotate ("; }
						  else { print $fh1 ",\\\n"; }
						  print $fh1 "\"".$mkname[$i]."\" ".$x;
						  $n++;
					 }
				}
				if($n) { print $fh1 ")\n"; }
				# Set xrange to correspond to map length
				print $fh1 "set xrange[$map_start[$lg][$sx]:$map_end[$lg][$sx]]\n";
				# Plot
				print $fh1 "plot \'",$name[$lg+$sx*$nch],"\' u 1:2 w l\n\n";
				# Wait for <return> (unless producing postscript)
				print $fh1 "pause -1\n" if(!$opt{o});
		  }
	 }
	 close $fh1;
	 # Call postscript if quiet option not set
	 system("gnuplot $name1") if(!$opt{q});
}
