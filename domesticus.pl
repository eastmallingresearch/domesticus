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
use Bio::Restriction::EnzymeCollection;
use Bio::Restriction::Analysis;
use Bio::Coordinate::Pair;
use Bio::Location::Simple;
use Bio::Coordinate::GeneMapper;

#DEFINE THE RESTRICTION SITES
my @enzymes=('BbsI','BsmBI','BsaI');

#READ IN THE DNA SEQUENCE
my $file         = shift; 
my $input_object = Bio::SeqIO->new(-file => $file);
my $input_seq= $input_object->next_seq;

my $locatable_seq =  Bio::LocatableSeq->new(-seq => $input_seq->seq,
                    -id  => "seq1");


#SANITY CHECK (SILENT FOR TRANSLATABLE SEQUENCE WITH NO STOPS OTHER THAN END)
my $input_prot= $locatable_seq->translate;
#print $input_prot->seq,"\n";


#GET A DATABASE OF RESTRICTION SITES
my $rebase = Bio::Restriction::IO->new(
      -file   => 'withrefm.404',
      -format => 'withrefm' );
my $rebase_collection = $rebase->read();

#DEFINE A CUSTOM COLLECTION OF RESTRICTION ENZYMES
my $custom_collection = Bio::Restriction::EnzymeCollection->new(-empty => 1);

#PUSH INTO CUSTOM COLLECTION FROM REBASE
foreach (@enzymes){
	#print "Retreving ". $_."\n";
	my $re=$rebase_collection->get_enzyme($_);
	#print $re->name();
	$custom_collection->enzymes($re);
}

#DEFINE A RESTRICTION ANALYSIS OBJECT
my $ra = Bio::Restriction::Analysis->new(-seq=>$locatable_seq, -enzymes=>$custom_collection);

#ANALYSIS
$ra->multiple_digest($custom_collection);

#DEFINE A HASH FOR STORING CUT SITE INFO
my %cut_hash=();

#LOOP AROUND THE ENZYMES AND PRINT OUT CUT STATS
foreach (@enzymes){
	print "CUTS BY $_ ";
	my $cut= $ra->cuts_by_enzyme($_);
	print ($cut/2,"\n");
	#IF THERE IS A CUT THEN ADD TO A HASH
		if ($cut/2 >0){
			my @cuts=$ra->positions($_);
			$cut_hash{$_}=\@cuts;
		}
}


#PRINT OUT SUMMARY OF DATA IN THE HASH
foreach (keys %cut_hash){
	print "STORED CUTS BY $_ ";
	my @cuts=@{$cut_hash{$_}};
	print join "\t", @cuts;
	print "\n";
}


####SEQUENCE COORDINATE CONVERSION

print "SEQUENCE COORDINATES \n\n";
print "GENE \t";
print $locatable_seq->start()."\t";
print $locatable_seq->end()."\n";
print "PROT \t";
print $input_prot->start()."\t";
print $input_prot->end()."\n";

my $map_start='11';
my $map_end='32';

print "GENEMAP\t";
print $map_start."\t";
print $map_end."\n";

#COORDINATE MAPPING

my $input_coordinates = Bio::Location::Simple->new  
  (-seq_id => 'cds', -start => $locatable_seq->start(), -end => $locatable_seq->end(), -strand=>1 );
my   $output_coordinates = Bio::Location::Simple->new  
	(-seq_id => 'pep', -start => $input_prot->start(), -end => $input_prot->end(), -strand=>1 );
my  $pair = Bio::Coordinate::Pair->new
  (-in => $input_coordinates ,  -out => $output_coordinates);
my  $pos = Bio::Location::Simple->new (-start => $map_start, -end => $map_end );

# create a gene mapper and set it to map from chromosomal to cds coordinates

my $gene = Bio::Coordinate::GeneMapper->new(-in   =>'cds',
                                      -out  =>'peptide',
                                      -cds  =>$locatable_seq,
                                     );
                                     
my $newloc = $gene->map($pos);
my $con_start = $newloc->start;
my $con_end = $newloc->end;

print "MAP\t";
print $con_start."\t";
print $con_end."\n";

exit;



print "VIRTUAL GEL OF FRAGMENTS \n" ;

my @gel;
  my @bam_maps = $ra->fragment_maps('multiple_digest');
  foreach my $i (@bam_maps) {
     my $start = $i->{start};
     my $end = $i->{end};
     my $sequence = $i->{seq};
    push @gel, "$start--$sequence--$end";
     @gel = sort {length $b <=> length $a} @gel;
  }
  print join("\n", @gel) . "\n";
  
  print "POSITIONS OF DIGESTIONS \n";
  
my @positions = $ra->positions('multiple_digest');
	foreach(@positions){
		print $_."\n";
	}
print "FRAGMENTS CREATED FROM DIGESTIONS\n";

my @fragments = $ra->fragments('multiple_digest');
	foreach(@fragments){
		print $_."\n";
	};


print "Multiple digest fragment lengths: ", join(' & ', map {length $_} @fragments), "\n";

  



  
