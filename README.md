##Goal:
The problem solved with these programs is the conversion of a vcf from one reference sequence to another. Because it relies on the output from Blast (http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&PROG_DEF=blastn&BLAST_PROG_DEF=megaBlast&BLAST_SPEC=blast2seq), it is best used for smaller sequences (ie. human mithochondrial dna). If you wish to convert more than one chromosome or sequence, the vcf must be split up and converted one chromosome or sequence at a time.

How to use:
1. Run the sequences through Blast and make sure that you out the original reference in first, so in the output it appears on top. Save the output as a file.
2. Run: python blast_parser.py <blast file>
3. Run: python convert_vcf.py <blast file> <vcf file>
4. The new vcf should have the same name as the old one with .conv at the end

##Before use:
Make sure that the vcf file contains only one chromosome and then (for now) edit the convert_vcf.py file on line 42 to say the chromosome you plan to convert.
Also, if there are any bad cals or loci you wish to exclude, you can edit the convert_vcf.py file and add those poisiotns to the bad loci list on line 10.

##Test:
In the test folder, I have provided a very short, made up vcf file (test.vcf) to test some of the edge cases to go with the example blast output (hg19_to_rsrs.txt) Along with those are the intermediate files that go along with running the scripts and the final outpu (test.vcf.conv)
