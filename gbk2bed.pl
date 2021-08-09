#!/usr/bin/perl
#perl script to make a bed file from a gbk file. The bed output is
#used in IGB sortware to view ChIP-seq data. 


use strict;
use Bio::SeqIO;
my $gbk = @ARGV[0];
my $in = Bio::SeqIO->new(-file=>"$gbk");
my $seq = $in->next_seq();
my $header = $seq->display_id(); 
foreach my $feat($seq->get_SeqFeatures()){
	my $start = $feat->start;
	my $end = $feat->end;
	my $strand = $feat->strand;
	my $clr;
	my $orient;
my $pritag = $feat->primary_tag();
if ($pritag eq "CDS" or $pritag eq "tRNA" or $pritag eq "rRNA"){
if ($feat->has_tag('product') and $feat->has_tag('locus_tag')){
my @gene = $feat->get_tag_values("locus_tag");
my ($product)= $feat->get_tag_values("product");
if ($strand ==1) {
$clr = '0,255,0'; $orient = '+';
my $name = $seq->display_id();
my @nn = split (/\s+/, $name);
my $name2 = join ("_", @nn);
$start =$start -1;
$end = $end-1;
print "$header\t$start\t$end\t$gene[0]\_$product\t1000\t$orient\t$start\t$end\t$clr\n";
}
else{
$clr = '0,0,255'; $orient = '-';
my $name = $seq->display_id();
my @nn = split (/\s+/, $name);
my $name2 = join ("_", @nn);
print "$header\t$start\t$end\t$gene[0]\_$product\t1000\t$orient\t$start\t$end\t$clr\n";
}
}
}
}
