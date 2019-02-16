#!/usr/bin/perl

use warnings;
use strict;

#
# concatenate fastq.gz separate by ";"

my $files = $ARGV[0]; # fichier chr	pos	sens	chr	pos	sens
my $out = $ARGV[1]; # absolut path

my @values = split(';', $files);

my $cmd="cat";
foreach my $val (@values) {
	chomp $val;
	
	$cmd=$cmd." \"".$val."\"";
}

$cmd=$cmd." > ".$out;
print($cmd."\n");
system($cmd);
