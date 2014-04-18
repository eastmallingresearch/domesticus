#!/usr/bin/perl


use warnings;
use strict;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SimpleAlign;
use Bio::SeqFeature::Generic;
use Bio::LocatableSeq;
use Bio::DB::SeqFeature;
use Bio::DB::SeqFeature::Store;
#use Bio::CodonUsage::Table;
#use Bio::DB::CUTG;
use Bio::Restriction::EnzymeCollection;
use Bio::Restriction::Analysis;

## get  a codon usage table from web database ##
#  my $cdtable = Bio::DB::CUTG->new(-sp => 'Mus musculus'
#	                                 -gc => 1);
#  ## or from local file
#  my $io = Bio::CodonUsage::IO->new(-file=>"file");
#  my $cdtable= $io->next_data();
#  ## or create your own from your own sequences 
#  ## get a Bio::PrimarySeq compliant object ##
#  # $codonstats is a ref to a hash of codon name /count key-value pairs.
#  my $codonstats = Bio::Tools::SeqUtils->codon_count($my_1ary_Seq_objct);
#  ### the '-data' field must be specified ##
#  ### the '-species' and 'genetic_code' fields are optional
# my $CUT = Bio::CodonUsage::Table->new(-data =>$codonstats,
#                                        -species => 'Hsapiens_kinase');
#  print "leu frequency is ", $cdtable->aa_frequency('LEU'), "\n";
#  print "freqof ATG is ", $cdtable->codon_rel_frequency('ttc'), "\n";
#  print "abs freq of ATG is ", $cdtable->codon_abs_frequency('ATG'), "\n";
#  print "number of ATG codons is ", $cdtable->codon_count('ATG'), "\n";
#  print "gc content at position 1 is ", $cdtable->get_coding_gc('1'), "\n";
#  print "total CDSs for Mus musculus  is ", $cdtable->cds_count(), "\n";

#READ IN THE DNA SEQUENCE
my $file         = shift; 
my $input_object = Bio::SeqIO->new(-file => $file);
my $input_seq= $input_object->next_seq;

my $locatable_seq =  Bio::LocatableSeq->new(-seq => $input_seq->seq,
                    -id  => "seq1",
                    -start => (),
                    -end   => ());

#print $locatable_seq->seq,"\n";

#SANITY CHECK (SILENT FOR TRANSLATABLE SEQUENCE WITH NO STOPS OTHER THAN END)
my $input_prot= $locatable_seq->translate;
#print $input_prot->seq,"\n";

my $rebase = Bio::Restriction::IO->new(
      -file   => 'withrefm.404',
      -format => 'withrefm' );
  my $rebase_collection = $rebase->read();

my $custom_collection = Bio::Restriction::EnzymeCollection->new(-empty => 1);
my $re=Bio::Restriction::Enzyme->new
      (-enzyme=>'BbsI', -seq=>'GAAGAC'); 
$custom_collection->enzymes($re);
      
        
#my $all_collection = Bio::Restriction::EnzymeCollection->new();
#my $BbsI_enzyme = $custom_collection->get_enzyme('BbsI');
my $BbsI_enzyme = $rebase_collection->get_enzyme('BbsI');

my $ra = Bio::Restriction::Analysis->new(-seq=>$locatable_seq, -enzymes=>$BbsI_enzyme);
print "BBS 1 cuts ".$BbsI_enzyme->cut."\n";
print "Recognition site ". $BbsI_enzyme->string."\n";
print "Overhang ".$BbsI_enzyme->overhang."\n";
#my $cbe = $ra->cuts_by_enzyme($BbsI_enzyme);

#print "CUTS BY ENZYME $cbe\n";

my @fragments = $ra->fragments($BbsI_enzyme);
	foreach(@fragments){
		print $_."\n";
	}
print "POSITIONS \n";

my @positions = $ra->positions('BbsI');
	foreach(@positions){
		print $_."\n";
	}
print "FRAGMENTS\n";
print "EcoRI fragment lengths: ", join(' & ', map {length $_} @fragments), "\n";


print "VIRTUAL GEL OF FRAGMENTS " ;

my @gel;
  my @bam_maps = $ra->fragment_maps('BbsI');
  foreach my $i (@bam_maps) {
     my $start = $i->{start};
     my $end = $i->{end};
     my $sequence = $i->{seq};
    push @gel, "$start--$sequence--$end";
     @gel = sort {length $b <=> length $a} @gel;
  }
  print join("\n", @gel) . "\n";
  
  
