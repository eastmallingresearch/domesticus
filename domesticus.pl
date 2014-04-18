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

my $locatable_seq =  Bio::LocatableSeq->new(-seq => $input_seq->seq);

#SANITY CHECK (SILENT FOR TRANSLATABLE SEQUENCE WITH NO STOPS OTHER THAN END)
my $input_prot= $locatable_seq->translate;
#print $input_prot->seq,"\n";

#GET ENZYMES FROM REBASE FILE
my $custom_collection=pull_enzymes(\@enzymes);

#DEFINE A RESTRICTION ANALYSIS OBJECT
my $ra = Bio::Restriction::Analysis->new(-seq=>$locatable_seq, -enzymes=>$custom_collection);
$ra->multiple_digest($custom_collection);

#CREATE GENERIC OBJECTS
my @seqobjects=generate_objects(\@enzymes,\$ra, \$locatable_seq);



#GENERATE THE INFO ABOUT STRAND AND RECOGNITION SITE

recognition_sites(\@seqobjects, \$custom_collection);


#ONCE ALL COORDINATES IN CHECK CODONS

my $sc=11;
my $ec=32;

#CONVERT COORDINATES
#my @coord=coord_convert(\$locatable_seq,\$input_prot,\$sc,\$ec);
#print "Converted coordinates ";
#print join "\t", @coord;
#print "\n";








#############SUBROUTINES##################

######RECOGNITION SITES#########
sub recognition_sites{
my ($seqobjects, $custom_collection)=@_;

foreach (@$seqobjects){
	#print $_->source_tag."\t";
	my $tag=$_->source_tag;
	print $tag."\n";
	#print $_->primary_tag."\t";
	my @one= $_->get_tag_values('one');
	my @two= $_->get_tag_values('two');
	#print join "\t", @one;
	#print "\t";
	#print join "\t", @two;
	#print "\t";
	#print $_->entire_seq()->seq;
	#print "\n";
	
	my $f_left=$one[0]-12;
	my $f_right=$two[0]+12;
	print $f_left."\t".$f_right."\n";
	my $enzyme= $$custom_collection->get_enzyme($tag);
	my $recog= $enzyme->recog();
	my $revcom_recog= $enzyme->revcom_recog();
	my $subseq=$_->entire_seq()->subseq($f_left,$f_right);
	print $subseq."\n";
	print "Searching for ".$recog. " and ".$revcom_recog."\n";
	
	if ($subseq=~ /$recog/){
                print "Forward found index location: $-[0]-$+[0]\n";
        }
         elsif ($subseq=~ /$revcom_recog/){
                print "Revcom found index location: $-[0]-$+[0]\n";
        }
         else {
                print "came in else\n";
        }
	
	}


}



####GENERATE SEQ OBJECTS######

sub generate_objects{
my ($enz,$ra,$loc)=@_;
#CUTTING HASH

my %cut_hash=find_cuts(\@$enz,\$$ra);
my @seqobjects=();

	foreach(keys %cut_hash){
		my $enz=$_;
		my $end =scalar @{($cut_hash{$_})}."\n";
		my @values=@{($cut_hash{$_})};
		
			for (my $i=0;$i<$end;$i=$i+2){
			
			#print "HERE ".$values[$i]."\n";
			#print "HERE ".$values[$i+1]."\n";
			
			
				my $feat = Bio::SeqFeature::Generic->new( 
 		 	           -primary      => 'restriction', # -primary_tag is a synonym
 		   	           -source_tag   => $enz,
 		     	       -display_name => '',
 		  	  	       -score        => 0,
        			   -tag          => { 'one' => $values[$i],
          			                      'two' => $values[$i+1] } );
          			$feat->attach_seq($$loc);                      
          			#print $$loc->seq;
          			 
				push (@seqobjects,$feat);
		}
}


return (@seqobjects);

}


############# ENZYMES FROM REBASE########

sub pull_enzymes{
my ($enzymes)=@_;

#GET A DATABASE OF RESTRICTION SITES
my $rebase = Bio::Restriction::IO->new(
      -file   => 'withrefm.404',
      -format => 'withrefm' );
my $rebase_collection = $rebase->read();

#DEFINE A CUSTOM COLLECTION OF RESTRICTION ENZYMES
my $custom_collection = Bio::Restriction::EnzymeCollection->new(-empty => 1);

#PUSH INTO CUSTOM COLLECTION FROM REBASE
foreach (@$enzymes){
	#print "Retreving ". $_."\n";
	my $re=$rebase_collection->get_enzyme($_);
	#print $re->name();
	$custom_collection->enzymes($re);
}
return $custom_collection;

}

############CUTS#############
sub find_cuts{
my ($enzymes,$ra)=@_;

my %hash=();

#LOOP AROUND THE ENZYMES AND PRINT OUT CUT STATS
foreach (@$enzymes){
	print "CUTS BY $_ ";
	my $cut= $$ra->cuts_by_enzyme($_);
	print ($cut/2,"\n");
	#IF THERE IS A CUT THEN ADD TO A HASH
		if ($cut/2 >0){
			my @cuts=$$ra->positions($_);
			$hash{$_}=\@cuts;
		}
}


#PRINT OUT SUMMARY OF DATA IN THE HASH
foreach (keys %hash){
	#print "STORED CUTS BY $_ ";
	my @cuts=@{$hash{$_}};
	#print join "\t", @cuts;
	#print "\n";
}
return %hash
}




######CONVERT COORDINATES DNA TO PROT##########
#sub coord_convert{
#my ($locatable_seq,$input_prot,$map_start,$map_end)=@_;

####SEQUENCE COORDINATE CONVERSION FROM DNA TO PROTEIN

#print "SEQUENCE COORDINATES \n\n";
#print "GENE \t";
#print $$locatable_seq->start()."\t";
#print $$locatable_seq->end()."\n";
#print "PROT \t";
#print $$input_prot->start()."\t";
#print $$input_prot->end()."\n";
#print "GENEMAP\t";
#print $$map_start."\t";
#print $$map_end."\n";

#COORDINATE MAPPING

#my $input_coordinates = Bio::Location::Simple->new  
#  (-seq_id => 'cds', -start => $$locatable_seq->start(), -end => $$locatable_seq->end(), -strand=>1 );
#my   $output_coordinates = Bio::Location::Simple->new  
#	(-seq_id => 'pep', -start => $$input_prot->start(), -end => $$input_prot->end(), -strand=>1 );
#my  $pair = Bio::Coordinate::Pair->new
#  (-in => $input_coordinates ,  -out => $output_coordinates);
#my  $pos = Bio::Location::Simple->new (-start => $$map_start, -end => $$map_end );

# create a gene mapper and set it to map from chromosomal to cds coordinates

#my $gene = Bio::Coordinate::GeneMapper->new(-in   =>'cds',
#                                      -out  =>'peptide',
#                                      -cds  =>$$locatable_seq,
#                                     );
                                     
#my $newloc = $gene->map($pos);
#my $con_start = $newloc->start;
#my $con_end = $newloc->end;

#print "MAP\t";
#print $con_start."\t";
#print $con_end."\n";
#return ($con_start,$con_end);

#}

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

  



  
