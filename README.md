This is a fork of VCF-Converter by KaiSmith. I made few changes and improved the documentation.

##Goal:
This script solve the problem of the conversion of a vcf variants coordinates from one reference genome version to another.

It relies on the output of Blast (http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&PROG_DEF=blastn&BLAST_PROG_DEF=megaBlast&BLAST_SPEC=blast2seq), it works better with small sequences (ie. human mithochondrial dna).

If you'd like to use the script with larger files (Whole exome or genome), we reccomend that you split the vcf file by chromosome or smaller chuncks.

##How to use:
1- Run the sequences through Blast and make sure that you print the  reference name first, so in the output file it will appears on the top.

2- Run: 

```
python blast_parser.py <blast file>
python convert_vcf.py <blast file> <vcf file>

```
	
3- The new vcf should have the same name as the old one with .conv at the end

##Before use:
Make sure that the vcf file contains only one chromosome and then (for now) edit the refchange.py file on line 42 to say the chromosome you plan to convert.

Also, if there are any bad cals or loci you wish to exclude, you can edit the refchange.py file and add those poisitions to the bad loci list on line 10.
