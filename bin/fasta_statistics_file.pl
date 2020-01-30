#!/bin/env perl
use strict;
use G;

my ($file, $list) = @ARGV;
my %fasta = readFile($file, -format=>"fasta");

my %list;
if( ~~@ARGV == 2 ){
    %list =  map{$_ => 1} readFile( $list, 1, -format=>"plain");
}

say "\n Stats for $file :";
my @keys = sort {length($fasta{$b})<=>length($fasta{$a})} keys %fasta;

my @len;
    for my $key (@keys){
	push(@len, length($fasta{$key}));

	next if ~~@ARGV ==2 && defined $list{$key};
        say $key."\t".length($fasta{$key});
}
my $totallen = sum(@len);

my ($n50, $n90, $n80, $n70, $n60, $n40, $n30, $n20, $n10);
my $i = 0;
for my $key (@len){
    $i += $key;
    if($i >= $totallen / 2){
	$n50 = $key;
	last;
    }
}

say "  Scaffold number           " . scalar(@keys);
say "  Total scaffold length     " . $totallen;
say "  Average scaffold length   " . int(mean(@len));
say "  Longest scaffold          " . max(@len);
say "  Shortest scaffold length  " . min(@len);
say "  N50                       " . $n50;


#####
###
##   Author : Gaou
###
####
