#!/usr/bin/perl
$k=0;
while(<>) {
  if(/Position/) {
    $k++;
	 next;
  }
	 next if(/Linkage/);
  split;
  $k1=@_;
  if($k1==4) {
    $id=$_[0]."_".$_[1];
    $y{$id}[$k]=$_[2];
  }
}
foreach $t(keys %y) {
  print $t," ";
  for($j=1;$j<=$k;$j++) {
    print $y{$t}[$j]+0," ";
  }
  print "\n";
}
