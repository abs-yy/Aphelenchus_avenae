#!/usr/bin/env perl
use strict;
use warnings;

my $ref = shift;
my $read1 = shift;
my $read2 = shift;
my $label = shift;
my $core = 64;
die("four arguments necessary: ref, read1, read2, label") unless(length($label));

unless(-e "$ref.bwt"){
    system("bwa index $ref");
}


system("bwa mem $ref $read1 $read2 -t $core  > $label.mem.sam");
system("samtools view -@ $core -bS $label.mem.sam > $label.mem.bam");
system("samtools sort -@ $core $label.mem.bam $label.mem.sorted");
system("samtools index $label.mem.sorted.bam");
system("rm $label.mem.sam $label.mem.bam")
system("samtools idxstat $label.mem.sorted.bam > $label.mem.sorted.bam.mapped")
