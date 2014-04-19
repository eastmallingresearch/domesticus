#!/usr/bin/perl


use warnings;
use strict;
use Bio::SeqFeature::Generic;
use Bio::LocatableSeq;
use Bio::Restriction::EnzymeCollection;
use Bio::Restriction::Analysis;
use Bio::Coordinate::Pair;
use Bio::Location::Simple;
use Bio::Coordinate::GeneMapper;

#DEFINE THE RESTRICTION SITES
my @enzymes=('BbsI','BsmBI','BsaI');
my %tails = (
        left_outer  => "tgaagacnnAAAA",
        left_inner   => "tgaagacnn",
        right_inner => "tgaagacnn",
        right_outer => "tgaagacnnTTTT"
    );
    
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
	
#CONVERT COORDINATES
coord_convert(\@seqobjects);

##MUTAGENESIS POSITION
mutagenesis(\@seqobjects);
print "\nMUTAGENESIS SUCCESSFUL- HAVE IDENTIFIED LESIONS \n";

#PRIMERS- DESIGN PRIMERS FOR THE MUTAGENESIS###
#FOR EACH LESION THE MUTATION MUST BE IN THE FIRST 4BP OF THE HOMOLOGOUS REGION OF THE 
#PRIMER- WHETHER FORWARD OR REVERSE

print "\nDESIGNING PRIMER SETS \n";

my %primers=design_primer(\@seqobjects,\$locatable_seq,\%tails);
	


#############SUBROUTINES##################

######DESIGN PRIMERS########NOTE I HAD TO INSTALL CLONE FROM CPAN TO MAKE THIS WORK

sub design_primer{
	my ($seqobj,$seq,$tails)=@_;
	my %results=();
	my $pairs=scalar(@$seqobj)+1;

	print scalar(@$seqobj)." cuts so ".$pairs." sets of primers\n"; 
	  
	
	  
	print "SEARCHING FOR PRIMERS\n";
	my @array_of_positions=();
	
	for (my $i=0; $i<=$pairs;$i++){
		
		if ($i==0){
			push (@array_of_positions,['x','x','x',$$seq->start()]);
		}	
		elsif ($i==$pairs){
			push (@array_of_positions,['x','x','x',$$seq->end()]);
		}	
		elsif ($i>0 && $i<=scalar(@$seqobj)){
			print "retrieving sequence ".$i."\n";
			my $feat=$$seqobj[$i-1];
			my @thirteen= $feat->get_tag_values('thirteen');
			foreach(@thirteen){
				my @tmp=@{$_};
					foreach (@tmp){
						print $_."\t";
						}
					push (@array_of_positions,\@tmp);				
				print "\n";
				}
			}
		else{
			print "ERROR in the design primer subroutine \n";
			}
		}
	
	my $mutagenised_seq=mutagenise_seq(\$$seq,\@array_of_positions);
	print "PRIMER ARRAY \n";
	for (my $i=0; $i<$pairs;$i++){
		print "PRIMER SET $i \t";
		my @primer_f=@{$array_of_positions[$i]};
		my @primer_r=@{$array_of_positions[$i+1]};
		
		if ($i==0){
			my $fw=($primer_f[3]);
			my $rv=($primer_r[3])+1;
			print $fw."\t";
			print $rv."\n";
			my @primers=primer_design($fw,$rv,\$mutagenised_seq);
			append_primers(\@primers,\'start',\%$tails);
		}
		elsif ($i==($pairs-1)){
			
			my $fw=($primer_f[3])-1;
			my $rv=($primer_r[3]);
			print $fw."\t";
			print $rv."\n";
			my @primers=primer_design($fw,$rv,\$mutagenised_seq);
			append_primers(\@primers,\'end',\%$tails);
		}
		elsif ($i>0 && $i<scalar(@$seqobj)){
			my $fw=($primer_f[3])-1;
			my $rv=($primer_r[3]+1);
			print $fw."\t";
			print $rv."\n";
			my @primers=primer_design($fw,$rv,\$mutagenised_seq);
			append_primers(\@primers,\'internal',\%$tails);
		}
		else{
			print "ERROR in the design primer subroutine \n";
			}
	}
	return %results;	
}


######APPEND THE TAILS TO THE PRIMERS########

sub append_primers{
	my($primers,$flag,$tails)=@_;
	my %hash=%{$tails};

	if ($$flag eq 'start'){
		print "PRIMERS $$flag \n";
		my $fwd=$hash{'left_outer'};
		my $rev=$hash{'right_inner'};
		print $fwd.$$primers[0]."\n";
		print $rev.$$primers[1]."\n";
	}
	elsif($$flag eq 'end'){
		print "PRIMERS $$flag \n";
		my $fwd=$hash{'left_inner'};
		my $rev=$hash{'right_outer'};
		print $fwd.$$primers[0]."\n";
		print $rev.$$primers[1]."\n";
	}
	elsif($$flag eq 'internal'){
		print "PRIMERS $$flag \n";
		my $fwd=$hash{'left_inner'};
		my $rev=$hash{'right_inner'};
		print $fwd.$$primers[0]."\n";
		print $rev.$$primers[1]."\n";
		}
	}

###EDIT SEQ TO CONTAIN ABOLISHED RESTRICTION SITES

sub primer_design{
	my($fw,$rv,$seq)=@_;
	use Bio::Tools::Run::Primer3;
	my @primers=();
			my $primer3 = Bio::Tools::Run::Primer3->new(-seq => $$seq,
		                                      -outfile => "temp.out",
		                                      -path => "/home/harrir/prog/primer3-2.3.6/"); 
			$primer3->add_targets('PRIMER_SEQUENCE_ID'=>"test",
			'SEQUENCE_FORCE_LEFT_START'=>$fw-1,
			'SEQUENCE_FORCE_RIGHT_START'=>$rv-1,
			'PRIMER_MAX_DIFF_TM' =>"3",
			'PRIMER_MAX_SIZE'=> "35",
			'PRIMER_PICK_ANYWAY'=>"1",
			'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' => '/home/harrir/prog/primer3-2.3.6/src/primer3_config/');
			my $results = $primer3->run;
			print "There were ", $results->number_of_results, " primers\n";
			my $all_results = $results->all_results;
			my %results=%{$all_results};
	 				foreach my $key (keys %{$all_results}) {
    					#print "$key\t${$all_results}{$key}\n";
    					#print $key."\n";
 					}
 					push(@primers,$results{'PRIMER_LEFT_0_SEQUENCE'});
 					push(@primers,$results{'PRIMER_RIGHT_0_SEQUENCE'});
 	return @primers;
	}

########GET THE SEQUENCE AND CHANGE THE POINT MUTATIONS TO THOSE IDENTIFIED

sub mutagenise_seq{
	my ($seq,$aop)=@_;
	my $nucl=$$seq->seq();
	my @nucleotides=split ('',$nucl);
	print "NUCLEOTIDE EDITOR SUBROUTINE \n";
		foreach (@$aop){
			my @split_data=@{$_};
				if ($split_data[2] ne 'x'){
					my $val=$split_data[3];
					print "editing nucleotide ".$split_data[3]." from \t";
					print $nucleotides[$val-1]."\t to ";
					print $split_data[2]."\n";
		 			$nucleotides[$val-1]=$split_data[2];
					}
			}
	my $tmp_seq=join "", @nucleotides;
	my $newseq =  Bio::LocatableSeq->new(-seq => $tmp_seq);
	return $newseq;
	}

########MUTAGENESIS##########

sub mutagenesis{
	my ($seqs)=@_;
	use Bio::Tools::CodonTable;
	#IN THIS SUBROUTINE, NOW THAT ALL THE COORDINATES ARE KNOWN AND THE FRAME OF THE 
	#CUT SITE IS UNDERSTOOD, SYNONYMOUS CHANGES ARE FOUND THAT WILL ABOLISH THE CUT 
	#SITE
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	print  join (' ', "The name of the codon table no.", $myCodonTable->id(1),
	       "is:", $myCodonTable->name(), "\n");
	foreach (@$seqs){
		my @one= $_->get_tag_values('one');
		my @two= $_->get_tag_values('two');
		my @three= $_->get_tag_values('three');
		my @four= $_->get_tag_values('four');
		my @five= $_->get_tag_values('five');
		my @six= $_->get_tag_values('six');
		my @seven= $_->get_tag_values('seven');
		my @eight= $_->get_tag_values('eight');
		my @nine= $_->get_tag_values('nine');
		my @ten= $_->get_tag_values('ten');
		my @eleven= $_->get_tag_values('eleven');
		my @twelve= $_->get_tag_values('twelve');
		my @split_eleven=split('',$eleven[0]);
		my @split_twelve=split('',$twelve[0]);
		my $counter=0; 
	    print "\n\n\nCODON MUTATION FOR ".$_->primary_tag." ".$_->source_tag." BASES BEING CONSIDERED\t";
		print scalar(@split_twelve)."\n";
		my @all_mutations=();
		# FOR EACH AA SITE IN PEPTIDE FIND MUTATIONS
			for (my $i=0;$i<scalar(@split_twelve);$i=$i+3){	
	 			my $codon =($split_twelve[$i].$split_twelve[$i+1].$split_twelve[$i+2]);
	 			#print "CODON ".$codon."\n";
	 			#print "AMINO ACID ".$split_eleven[$counter]."\n";
	 			#print "EXAMINING CODON TABLE\n";
	 			my @codons = $myCodonTable->revtranslate($split_eleven[$counter]);
  	 			#print join "\t", @codons;
				#print "\n";
				#IF WE ARE AT A SITE WITH >1 FOLD DEGENERACY THEN LIST THE MUTATIONS AND 
				#THE POSITION OF THE MUTATION IN AN ARRAY
				if (scalar @codons>1){
					#PROVIDE THIS SUBROUTINE WITH DETAILS OF PEPTIDE, RECOGNITION SITE AND POSITION IN CURRENT LOOP
					my @mutations=mutator(\@codons,\$codon,\$i,\$three[0],\$four[0],\$seven[0],\$eight[0]);
					#print "MUTATION AT SITE WITHIN RECOGNITION SEQUENCE ".$counter."\n";
					foreach(@mutations){
							push (@all_mutations,$_);
						}
					}
				$counter=$counter+1;
			}
		###NOW YOU HAVE A PILE OF POSSIBLE MUTATIONS, PICK ONE AND ADD COORDINATES
		print "POSSIBLE MUTATIONS ".scalar(@all_mutations)."\n";
		print "PICKING NUMBER BETWEEN 0 AND ".(scalar(@all_mutations)-1)."\n";
		my $random_number = int(rand((scalar(@all_mutations)-1)));
		print "PICKED ".$random_number."\n";
		my @pick=();
		my $j=0;
		foreach(@all_mutations){		
						my @tmp=@{$_};
						if ($j==$random_number){
							@pick=@tmp;
							}
							$j++;
						}
		$_->add_tag_value('thirteen',\@pick);
		my @thirteen= $_->get_tag_values('thirteen');
		
		foreach(@thirteen){
			my @tmp=@{$_};
				foreach (@tmp){
						print $_."\t";
						}					
			print "\n";
			}
		}
	}


########MUTATOR- THIS TAKES THE CODONS- THE CURRENT CODON AND RETURNS A MUTATION- DEPENDS
#UPON THE MUTAGENESIS SUBROUTINE

sub mutator {

	my ($codons,$codon,$pos,$rec_s,$rec_e,$pep_s,$pep_e)=@_;
	my @mutated=();
	#print "Considering gene range ".$$rec_s."\t".$$rec_e."\t peptide range ".$$pep_s."\t".$$pep_e."\n";
	foreach (@$codons){
		if ($_ !~ m/$$codon/i){
			#print $_." is not ".$$codon."\n";
			my $cod=$_;
			my @s1 = split //,$$codon;
			my @s2 = split //,$_;
			my $i = 0;
			foreach  (@s1) {
				my $mut=$s2[$i];
    			if ($_ !~ m/$mut/i) {
    				my $position=$i+$$pos+$$pep_s;
       			 	#print "$_, $s2[$i] $i $$pos $$pep_s\t".$position."\n";
       			 	if ($position >= $$rec_s && $position <= $$rec_e){
       			 	#	print "GOOD SITE\n";
       			 		my @tmp=();
       			 		push (@tmp,$$codon);
       			 		push (@tmp,$cod);
       			 		push (@tmp,$s2[$i]);	
       			 		push (@tmp,$position);
       			 		push (@mutated,\@tmp);
       			 		}
       			 		else{
       			 			print "BAD SITE IN MUTATOR SUBROUTINE\n";
       			 			}
    				}
   				 $i++;
				}
			}
		}
	return @mutated;
	}

######RECOGNITION SITES#########
# THIS SCRIPT FINDS THE LOCATION OF THE RECOGNITION SITE IN AN EXTENDED FRAGMENT OF DNA
# NOTE- THIS TAKES AN ARBITRARY REGION SURROUNDING THE CUT SITE OF 12 BP EITHER SIDE

sub recognition_sites{
	my ($seqobjects, $custom_collection)=@_;

	foreach (@$seqobjects){
		#print "RECOGNITION SITES\n";
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
		my $strand=0;
			if ($subseq=~ /$recog/){
                print "Forward found index location: $-[0]-$+[0]\n";
                $strand =1;
                $f_right=$f_left+$+[0]-1;
                $f_left=$f_left+$-[0];
        	}
         	elsif ($subseq=~ /$revcom_recog/){
                print "Revcom found index location: $-[0]-$+[0]\n";
                $strand=-1;
                $f_right=$f_left+$+[0]-1;
                $f_left=$f_left+$-[0];
        	}
         	else {
                print "came in else\n";
        	}
		print $f_left."\t".$f_right."\n";
		#my $sanitycheck=$_->entire_seq()->subseq($f_left,$f_right);
		#print $sanitycheck."\n";
		$_->add_tag_value('three',$f_left);
		$_->add_tag_value('four',$f_right);
		$_->strand($strand)
		}
	}

####GENERATE SEQ OBJECTS######
#THIS INSTANTIATES THE GENERIC SEQUENCE OBJECTS THAT ARE SUBSEQUENTLY USED TO STORE
#EVERYTHING

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
          			 $feat->start($$loc->start);
          			 $feat->end($$loc->end); 
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
		#print join "\t", @cuts;§
		#print "\n";
		}
	return %hash
	}

######CONVERT COORDINATES DNA TO PROT##########

sub coord_convert{
	my ($seqs)=@_;
	foreach (@$seqs){
		my $input_prot=$_->entire_seq()->translate();
		my @three= $_->get_tag_values('three');
		my @four= $_->get_tag_values('four');
		#print $input_prot->seq;
		print "\nSEQUENCE COORDINATES \n\n";
		print "GENE \t";
		print $_->start."\t";
		print $_->end."\n";
		print "PROT \t";
		print $input_prot->start()."\t";
		print $input_prot->end()."\n";
		print "GENEMAP\t";
		print $three[0]."\t";
		print $four[0]."\n";
		#COORDINATE MAPPING
		my $input_coordinates = Bio::Location::Simple->new  
	  	(-seq_id => 'cds', -start => $_->start, -end => $_->end(), -strand=>$_->strand() );
		my   $output_coordinates = Bio::Location::Simple->new  
		(-seq_id => 'pep', -start => $input_prot->start(), -end => $input_prot->end(), -strand=>$_->strand() );
		my  $pair = Bio::Coordinate::Pair->new
		(-in => $input_coordinates ,  -out => $output_coordinates);
		my  $pos = Bio::Location::Simple->new (-start =>  $three[0]-1, -end =>  $four[0]-1);
		# create a gene mapper and set it to map from chromosomal to cds coordinates
		my $gene = Bio::Coordinate::GeneMapper->new(-in   =>'cds',
	                                      -out  =>'peptide',
	                                      -cds  =>$_->entire_seq(),
	                                     );
                                     
		my $newloc = $gene->map($pos);
		my $con_start = $newloc->start;
		my $con_end = $newloc->end;
		print "MAP\t";
		print $con_start."\t";
		print $con_end."\n";
		$_->add_tag_value('five',$con_start);
		$_->add_tag_value('six',$con_end);
	
		my $gene_back = Bio::Coordinate::GeneMapper->new(-in   =>'peptide',
	                                      -out  =>'cds',
	                                      -cds  =>$input_prot,
	                                     );
   		my $pos_back = Bio::Location::Simple->new (-start =>  $con_start, -end =>  $con_end );                          
		my $newloc_back = $gene_back->map($pos_back);
		my $pep_start = $newloc_back->start;
		my $pep_end = $newloc_back->end;
		print "PEP\t";
		print $pep_start."\t";
		print $pep_end."\n";
		$_->add_tag_value('seven',$pep_start);
		$_->add_tag_value('eight',$pep_end);
		my $frame_shift= $pep_start- $three[0];
		print "FRAME \t". $frame_shift."\n";
		$_->add_tag_value('nine',$frame_shift);
		$_->add_tag_value('ten',$input_prot);
		my $small_peptide=$input_prot->subseq($con_start,$con_end);
		$_->add_tag_value('eleven',$small_peptide);
		print "PEPTIDE \t". $small_peptide."\n";
		my $small_seq=$_->entire_seq()->subseq($pep_start,$pep_end);
		print "NUC SEQ \t". $small_seq."\n";
		$_->add_tag_value('twelve',$small_seq);	
		}
	}



#primary      
#restriction
#source_tag   
#enzyme name
#strand
#orientation of the recognition site	   	           		   	           
#one 
#start of the cut site
#two
#end of the cut site
#three
#start of the recognition site
#four
#end of the recognition site
#five
#protein start
#six
#protein end
#seven
#peptide start on DNA sequence
#eight
#peptide end on DNA sequence
#nine
#FRAME of DNA/PROTEIN (e.g. -2 is where the peptide starts 2bp before the recognition site
##ten
#PROTEIN SEQ
#eleven
#PEPTIDE
#twelve
#NUCLEOTIDE IN FRAME (within which is recognition site)
#thirteen
#array of arrays with following format
# CODON_OLD	CODON_NEW	NUC	NUC_POS

  
