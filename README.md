# Polyploidism
Estimate the polyploidy


                ||----P |
                ||     ||  
	 Polyploidism v0.7 [Oct 15, 2017]  - [+]

Citation - automated estimation of polyploidy nature of the genome

License: Creative Commons Licence
Bug-reports and requests to: jitendra.narayanATunamur.be and nicolas.debortoliATunamur.be

Keep an eye at following line: 
#Set up the path of all required tools - stored in the same folder
$ENV{'PATH'} = "/bin:/usr/bin:/usr/bin/env:$ENV{PWD}/augustus.2.5.5/bin: $ENV{PWD}/augustus.2.5.5/config: $ENV{PWD}/lastz-distrib-1.03.73/src: $ENV{PWD}/KaKs_Calculator2.0/bin/Linux: $ENV{PWD}/muscle: $ENV{PWD}/dotter: $ENV{PWD}/trf";


FOR GENE try: 
perl Polyploidism.pl -f xaa.fa -o TESTOUT3 -s elegans -q 50 -e 4 -g -i 70 -t 7 -z 5

FOR RANDOM: 
perl Polyploidism.pl -f sampleSeq.fa -o TESTOUT3 -s elegans -q 50 -e 4 -r -i 70 -t 7 -l 500 -c 2 -z 5

