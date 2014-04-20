domesticus
==========

Scripts to tame wild gene sequences
Usage:

domesticus.pl 'fastafile' 'arguments for this domestication' 'all arguments' 'restriction sites for this dometication' 'parameters paths etc'

eg. ./domesticus.pl Sdd1-C3.fasta args.txt sites.txt enzymes.txt parameters.txt


Note- at present this relies on a customised version of the Primer3.pm module included with Bioperl. This is included in this repository, but it's a clunky fix- Bioperl really needs updating to a recent verion of primer3. This is important because the script relies on some of the later functions in primer3, not present in 1.x versions. 

