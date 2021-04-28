#!/usr/bin/perl

################## Getopt stuff #############################################################################
use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::SearchIO;
my $faa_file;
my $gene_names;
my $db_list;
my $reference_faa;
my $sep;
GetOptions (
"faa:s" => \$faa_file,
"names:s" => \$gene_names,
"db_list=s" => \$db_list,
"ref=s" => \$reference_faa,
"sep=s" => \$sep);



####################### pod stuff #################################################################################
=head1 NAME
two-way-blast2.pl

=head1 SYNOPSIS

Reciprocal_blast hit generator. The output is a comma separated table that contains the fraction identical ratio of proteins found in the query blast database in comparison to the reference. The latter is either a list of files or a fasta file containing gene names and their amino acid sequences. 
The output can be thought of as comparative genomic analysis of a set of genes. Heatmaps can be drawn using the output of this programme as an input for R scripts.
I<This programme requires that you have already made a blast data base per each genome. Read about NCBI blast.>

=head2

B<note:> your gene names must be the locus_tag exactly as they appear in the gbk file

=head1 USAGE

 perl two-way-blast.pl -db_list <file name> -ref <fasta file>  -sep {tab | comma}  {-faa <file name> | -names <file name>} 

=over 2

=item B<-faa>

use this option if you have a fasta file as an input for the query genes

=item B<-names>

use this option if you only have the gene names as a list. The programme will generate the fasta files from the reference file that you provide

=item B<-ref>

name of and path to the reference B<fasta file> file

=item B<-db_list>

A list of the blast database names that will be used to find reciprocal blast hits. 

=item B<-sep>

tab- or comma-delemited output

=back

=cut

unless ($db_list  and $reference_faa and $sep and ( $faa_file  or  $gene_names)) {
	exec("perldoc $0");
}

############ making a tab delimited scalar for organism names ############
my @ar;
my $header;
open (FILE, $db_list);
while (<FILE>){
chomp;
push @ar, $_;
if ($sep eq "tab"){
$header = join ("\t","gene_id", @ar);
}
elsif ($sep eq "comma"){
$header = join("\,", "gene_id",  @ar);
}
else{print "sep should either be tab or comma\n"; exit;}
}
close(FILE);

############ filtering on options whether reference file is faa file or just gene names #################
my $in1;
if (defined $faa_file){
 $in1 = Bio::SeqIO->new (-file=>"$faa_file", -format=>'fasta'); # <file.faa> contains the reference fasta sequences. 
}
elsif (defined $gene_names){
get_fasta ($reference_faa, $gene_names);
$in1 = Bio::SeqIO->new (-file=>"genes_as_fasta.faa", -format=>'fasta'); # regulators.faa contains the reference fasta sequences. 
}


print "$header\n";
########### writing each fasta entry into a file that contains one entry at a time ######################
while (my $obj1 = $in1->next_seq()){
open (REG, ">one_reg.fasta") or die "can't open file"; #this file keeps overwritten with new entries. 
my $Rg =$obj1->display_name();
my $Rgseq = $obj1->seq();
print REG ">$Rg\n$Rgseq\n";   
my @line;
my $query;

########### opening the list that contains database names and making the first blast ####################
open (LIST, $db_list) or die "can't open file";
   while (<LIST>){
   chomp;
   my $db = $_;
`blastn -query one_reg.fasta -db $db -out temp1`;
    
########### parsing the blast result file ###############################################################
	my $searchio = Bio::SearchIO->new( -file   => 'temp1' );
          while ( my $result = $searchio->next_result() ) {
		
		my $hit = $result->next_hit;
             	if (defined $hit) {
		my $hsp = $hit->next_hsp;	
		my $frqid =$hsp->frac_identical;
		my $score = $hsp->score;
		my $bits = $hsp->bits;
		my $frqdec = sprintf ("%.2f", $frqid);
		my $result_name = $hit->name;
		my $query_len = $obj1->length;
		my $len = $hit->length;
		my $significance = $hit->significance();
		$query = $result->query_name();

########### performing out reciprocal blast search ##########################################################
		my $recip = recip_blast($result_name, $db);
			if ($query eq $recip){
			#my $parameter = "ORTHOLOGUE";
			push @line,  $frqdec;
			push @line, $result_name;
			}
			else {
                	#my $parameter = "NOT-ORTHO";
	       		 $score = 0; $bits = 0; $frqdec = 0.01;
			push @line, $frqdec;
			push @line, 'NaN';
			}
		}
		
############# if there is no hit in the first blast print "NaN" this term is understood by R #############
		elsif (!defined $hit){
		my $frqdec = "NaN";
		push @line, $frqdec;
		push @line, "NaN";
		}
	}
   }


############# finally printing the results as comma or tab separated to STDOUT ###########################
unshift @line, $query;
if ($sep eq "comma"){
print (join(",", @line),"\n");	
close (REG);
}
elsif ($sep eq "tab"){
print (join("\t", @line), "\n");
close (REG);
}
else {
print "$sep is not an allowed argument ('tab' or 'comma' are the only valid arguments)\n\n\n";
exit;
}
}






############################## SUBROUTINE: recip_blas t##############################

sub recip_blast {
my ($hit2, $dbs)= @_;
my $in = Bio::SeqIO->new (-file => "$dbs", -format=> "fasta");
   while (my $obj = $in->next_seq()){
	open (OUT, ">hits.fasta");
	my $seq = $obj->seq();
	my $locus = $obj->display_id();
		if ($hit2 eq $locus){
		print OUT ">$locus\n$seq\n";
		`blastn -query hits.fasta -db $reference_faa -out temp2`;
  		my $searchio = Bio::SearchIO->new (-file=>"temp2");
  		my $resultobj = $searchio->next_result;
          		if (defined $resultobj){
          		my $hitobj = $resultobj->next_hit;
          		my $recipname = $hitobj->name if (defined $hitobj);		
          		return $recipname;
          		}       
		}
   }
close (OUT);
}
############################# SUBROUTINE: get_fasta   ####################################
sub get_fasta{
my %hash;
my $genome_file = $_[0];
my $gene_list = $_[1];
my $in = Bio::SeqIO->new(-file=>"$genome_file", -format=>"fasta");
while (my $obj = $in->next_seq()){
	my $locus = $obj->display_name();
	my $translation = $obj->seq();        
$hash{$locus} = $translation;
        }
open (TEMP, ">genes_as_fasta.faa");
open (GENE, $gene_list);
while (<GENE>){
chomp;
if (exists $hash{$_}){ 
print TEMP ">$_\n$hash{$_}\n";
}
else{
print STDERR "
=====================================================

$_ locus doesn't exist in the reference file!

=====================================================\n";
}
}
close (GENE);
}	
