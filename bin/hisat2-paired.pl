use G;
use strict;

my $ref = shift;
my $read1 = shift;
my $read2 = shift;
my $label = shift;
my $core = 32;
die("three arguments necessary: ref, read1, label") unless(length($label));

unless(-e "$ref.1.ht2"){
    system("/path/to/hisat2-2.1.0/hisat2-build -p 32 $ref $ref");
}

/path/to/hisat2-2.1.0/hisat2 -x $ref -1 $read1 -2 $read2 -p 32 -S $label.hisat2.sam
/path/to/samtools view -@ $core -bS $label.hisat2.sam > $label.hisat2.bam
/path/to/samtools sort -@ $core -m 6G -o $label.hisat2.sorted.bam $label.hisat2.bam
/path/to/samtools index $label.hisat2.sorted.bam
