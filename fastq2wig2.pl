#! /usr/bin/perl
use strict;

my $fastq;
my $index;
my $mm;
my $out;
my $norm;
my @one;
my @more;
my $chromo;
my $gsize;
use Getopt::Long;


GetOptions (
"f|filename=s" => \$fastq,
"i=s" => \$index,
"m=s" => \$mm,
"n=s" => \$norm,
"o=s" => \$out,
"c=s" => \$chromo,
"gs=s" => \$gsize);

if ($fastq eq '' || $index eq '' || $mm eq '' || $norm eq '' || $out eq '' || $chromo eq '' || $gsize eq ''){
        warn "Usage: -f <fastq file> -i <index> -m <mismatch> -n <normalization y or n> -o <output name> -c <chromosome name> -gs <int genome size>\n";
        exit;
}


############fastq to sam##############
`bowtie2 -x $index -U $fastq -p 8 -S $out.sam 2>$out.bow`;


############calculate total mapped reads number############
open (BOW,"$out.bow");

while (<BOW>){
        chomp;

        if ($_ =~/exactly 1 time/){
                @one = split (/\s+/, $_);
#               print "$one[1]\n";
        }elsif ($_ =~/1 times/){
                @more = split (/\s+/, $_);
#               print "$more[1]\n";
        }
}
close IN;

my $total = $one[1] + $more[1];


############sam to mis.sam#############
`grep -e "^@" -e "XM:i:[0-$mm][^0-9]" $out.sam >$out\_$mm.sam`;
`rm $out.sam`;
=cut
############mis.sam to take5.sam#############
open (IN, "$out\_$mm.sam");
open (OUT, ">$out\_$mm\_take5.sam");
while (<IN>){
        chomp;
        my @sam = split (/\t/, $_);

        if ($sam[0] =~ m/^@/){
                print OUT (join ("\t", @sam,"\n"));
        }elsif ($sam[1] == 0){
                my $seq = substr ($sam[9], 0, 1);
                my $qual = substr ($sam[10], 0, 1);
                print OUT "$sam[0]\t$sam[1]\t$sam[2]\t$sam[3]\t$sam[4]\t1M\t$sam[6]\t$sam[7]\t$sam[8]\t$seq\t$qual\t$sam[11]\t$sam[12]\t$sam[13]\t$sam[14]\t$sam[15]\t$sam[16]\t$sam[17]\t$sam[18]\n";
        }elsif ($sam[1] == 16){
                my $len = length ($sam[9]);
                my $seq = substr ($sam[9], -1, 1);
                my $qual = substr ($sam[10], -1, 1);
                my $pos = $sam[3]+$len-1;
                print OUT "$sam[0]\t$sam[1]\t$sam[2]\t$pos\t$sam[4]\t1M\t$sam[6]\t$sam[7]\t$sam[8]\t$seq\t$qual\t$sam[11]\t$sam[12]\t$sam[13]\t$sam[14]\t$sam[15]\t$sam[16]\t$sam[17]\t$sam[18]\n";
        }
}
close IN;
close OUT;
=cut

############take5.sam to take5.bam###########
`samtools view -bS $out\_$mm.sam >$out\_$mm.bam`;
`rm $out\_$mm.sam`;


############take5.bam to sorted_take5.bam########
`samtools sort $out\_$mm.bam -o  sorted_$out\_$mm.bam`;
`rm $out\_$mm.bam`;

############sorted_take5.bam to sorted_take5.gff############

open (BAM, "samtools mpileup -d 1000000 sorted_$out\_$mm.bam|");
open (GFF, ">sorted_$out\_$mm.gff");
   while (<BAM>){
        chomp;
        my @sp = split(/\t/, $_);

            my $countp = ($sp[4] =~ tr/[ATGC]//);
            my $countm = ($sp[4] =~ tr/[atgc]//);


	if ($norm eq 'y'){
			my $normp = sprintf ("%.2f", $sp[3] / $total * 1000000);
                        my $normm = sprintf ("%.2f", $sp[3] / $total * 1000000);
              

		if ($countp != 0){
                        print GFF "$sp[0]\ttss\t$out\_fwd\t$sp[1]\t$sp[1]\t$normp\t+\t.\t.\n";

                }if ($countm != 0){
                        print GFF "$sp[0]\ttss\t$out\_rev\t$sp[1]\t$sp[1]\t-$normm\t-\t.\t.\n";
                }
        }elsif ($norm eq 'n'){
                if ($countp != 0){
                        print GFF "$sp[0]\ttss\t$out\_fwd\t$sp[1]\t$sp[1]\t$sp[3]\t+\t.\t.\n";
                }
                if ($countm != 0){
                        print GFF "$sp[0]\ttss\t$out\_rev\t$sp[1]\t$sp[1]\t-$sp[3]\t-\t.\t.\n";
                }
        }
}
close BAM;
close GFF;

my %coordrev;
my %coordfwd;
my @ar = (1..$gsize);
foreach my $key (@ar){
$coordrev{$key} = 0;
$coordfwd{$key} = 0;
}
my $hashsize= scalar keys %coordfwd;
print "$hashsize\n";
foreach my $rev (keys %coordrev){
$coordrev{$rev} = 0;
}
foreach my $fwd (keys %coordfwd){
$coordfwd{$fwd} = 0
}
open (GFF1, "sorted_$out\_$mm.gff") or die "no such file\n";
open (WIGmin, ">sorted_$out\_$mm\_neg.wig");
open (WIGpos, ">sorted_$out\_$mm\_pos.wig");
print WIGmin "track   type=wiggle_0   name=sorted_$out\_$mm\_neg      graphType=points        visibility=full color=168,130,88
fixedStep       chrom=$chromo start=1 step=1  span=1\n";
print WIGpos "track   type=wiggle_0   name=sorted_$out\_$mm\_pos      graphType=points        visibility=full color=168,130,88
fixedStep       chrom=$chromo start=1 step=1  span=1\n";

while (<GFF1>){
chomp;
my @ar = split (/\t/, $_);
my $strand = substr $ar[2], -3;
my $position = $ar[3];
	if ($strand eq "rev" and exists $coordrev{$position}){
	$coordrev{$position} = $ar[5];
	}
	elsif ($strand eq "fwd" and exists $coordfwd{$position}){
	$coordfwd{$position} = $ar[5];
	}

}

foreach my $neg (sort {$a<=>$b} keys %coordrev){
print WIGmin "$coordrev{$neg}\n";
}
foreach my $pos (sort {$a<=>$b} keys %coordfwd){
print WIGpos "$coordfwd{$pos}\n";
}
close (GFF1); close (WIGmin); close (WIGpos);






