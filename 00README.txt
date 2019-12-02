EPN, Mon Dec  2 07:01:08 2019

Version 0.05

cmsearch-deoverlap.pl: remove lower scoring overlaps from cmsearch
                       --tblout files.

Usage:
$ perl ./cmsearch-deoverlap.pl 
cmsearch-deoverlap v0.05 [Dec 2019]

Usage:

cmsearch-deoverlap.pl    [OPTIONS] <tblout file>
	OR
cmsearch-deoverlap.pl -l [OPTIONS] <list of tblout files>

	OPTIONS:
		-l             : single command line argument is a list of tblout files, not a single tblout file
		-s             : sort hits by bit score [default: sort by E-value]
		-d             : run in debugging mode (prints extra info)
		-v             : run in verbose mode (prints all removed and kept hits)
		--noverlap <n> : define an overlap as >= <n> or more overlapping residues [1]
		--nhmmer       : tblout files are from nhmmer v3.x
		--hmmsearch    : tblout files are from hmmsearch v3.x
		--cmscan       : tblout files are from cmscan v1.1x, not cmsearch
		--besthmm      : with --hmmsearch, sort by evalue/score of *best* single hit not evalue/score of full seq
		--clanin <s>   : only remove overlaps within clans, read clan info from file <s> [default: remove all overlaps]
		--maxkeep      : keep hits that only overlap with other hits that are not kept [default: remove all hits with higher scoring overlap]
		--dirty        : keep intermediate files (sorted tblout files)


--------
EXAMPLES
--------

1. Remove overlapping hits to models in the same clan in a single tblout
file:
$ perl ./cmsearch-deoverlap.pl --clanin ribo.claninfo testfiles/1.cmsearch.tblout
Saved   357 hits (  614 removed) to testfiles/1.cmsearch.tblout.deoverlapped

--
2. Remove overlapping hits to models in the same clan in each tblout
files from a list:

$ perl ./cmsearch-deoverlap.pl --clanin ribo.claninfo -l testfiles/file.list
Saved   357 hits (  614 removed) to testfiles/1.cmsearch.tblout.deoverlapped
Saved   325 hits (  596 removed) to testfiles/2.cmsearch.tblout.deoverlapped
Saved   377 hits (  644 removed) to testfiles/3.cmsearch.tblout.deoverlapped

--
3. Remove all overlapping hits in a single tblout file:

$ perl ./cmsearch-deoverlap.pl testfiles/1.cmsearch.tblout
Saved   357 hits (  614 removed) to testfiles/1.cmsearch.tblout.deoverlapped

-------
TESTING
-------
There are cmsearch and cmscan tblout files in the directory testfiles/
that were created with this bsub script using Infernal 1.1.2:

for i in 1 2 3; do 
    bsub -n 4 -J cmsearch.$i -o /dev/null -e $i.cmsearch.err -M 12000 -R "rusage[mem=12000]" "time cmsearch -Z 1000 --hmmonly --cut_ga --cpu 4 --noali --tblout $i.cmsearch.tblout ribo.cm 50000.$i.fa >  $i.cmsearch.out"
    bsub -n 4 -J cmscan.clan.$i -o /dev/null -e $i.cmsearch.err -M 12000 -R "rusage[mem=12000]" "time cmscan -Z 1000 --hmmonly --cut_ga --cpu 4 --noali --fmt 2 --oskip --oclan --clanin ribo.claninfo --tblout $i.cmscan.clan.tblout ribo.cm 50000.$i.fa >  $i.cmscan.clan.out"
    bsub -n 4 -J cmscan.clan.$i -o /dev/null -e $i.cmsearch.err -M 12000 -R "rusage[mem=12000]" "time cmscan -Z 1000 --hmmonly --cut_ga --cpu 4 --noali --fmt 2 --oskip --tblout $i.cmscan.noclan.tblout ribo.cm 50000.$i.fa >  $i.cmscan.noclan.out"
done

The models used were Rfam 12.2 models, and here's the cmstat outptu:
$ cmstat ribo.cm
# cmstat :: display summary statistics for CMs
# INFERNAL 1.1.2 (July 2016)
# Copyright (C) 2016 Howard Hughes Medical Institute.
# Freely distributed under a BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                                                                                              rel entropy
#                                                                                             ------------
# idx   name                  accession      nseq  eff_nseq   clen      W   bps  bifs  model     cm    hmm
# ----  --------------------  ---------  --------  --------  -----  -----  ----  ----  -----  -----  -----
     1  5S_rRNA               RF00001         712      7.35    119    194    34     1     cm  0.590  0.370
     2  5_8S_rRNA             RF00002          61      3.62    154    203    25     3     cm  0.590  0.471
     3  SSU_rRNA_bacteria     RF00177          99      1.21   1533   1614   462    31     cm  0.590  0.316
     4  SSU_rRNA_archaea      RF01959          86      1.04   1477   1646   457    30     cm  0.590  0.302
     5  SSU_rRNA_eukarya      RF01960          91      3.61   1851   2459   447    30     cm  0.590  0.418
     6  LSU_rRNA_archaea      RF02540          91      1.52   2990   3855   786    68     cm  0.590  0.347
     7  LSU_rRNA_bacteria     RF02541         102      1.36   2925   3841   846    70     cm  0.590  0.325
     8  SSU_rRNA_microsporidia  RF02542          46      1.70   1312   1570   366    26     cm  0.590  0.360
     9  LSU_rRNA_eukarya      RF02543          88      2.40   3401   4470   872    71     cm  0.590  0.385
    10  SSU_trypano_mito      RF02545           4      1.15    623    651   175     7     cm  0.590  0.337
    11  LSU_trypano_mito      RF02546           6      1.74    561    586    58     6     cm  0.591  0.500
    12  mtPerm-5S             RF02547           6      0.99    106    124    28     2     cm  0.590  0.347
#

The file 'test.sh' will run a set of tests that check that the main
script cmsearch-deoverlap.pl is working properly for parsing cmsearch
tblout files. *It does not currently do any checks to see if the 
--nhmmer or --hmmsearch options are working properly for removing
overlaps from nhmmer and hmmsearch output.

Here is that file:
$ cat test.sh
perl ./check-cmsearch-v-cmscan.pl testfiles/1.cmsearch.tblout testfiles/1.cmscan.noclan.tblout
perl ./check-cmsearch-v-cmscan.pl testfiles/2.cmsearch.tblout testfiles/2.cmscan.noclan.tblout
perl ./check-cmsearch-v-cmscan.pl testfiles/3.cmsearch.tblout testfiles/3.cmscan.noclan.tblout

perl ./check-cmsearch-v-cmscan.pl --clanin ribo.claninfo testfiles/1.cmsearch.tblout testfiles/1.cmscan.clan.tblout
perl ./check-cmsearch-v-cmscan.pl --clanin ribo.claninfo testfiles/2.cmsearch.tblout testfiles/2.cmscan.clan.tblout
perl ./check-cmsearch-v-cmscan.pl --clanin ribo.claninfo testfiles/3.cmsearch.tblout testfiles/3.cmscan.clan.tblout

Here is the expected output of that script:
$ sh test.sh
Saved   357 hits (  614 removed) to testfiles/1.cmsearch.tblout.deoverlapped
Files are identical
Saved   325 hits (  596 removed) to testfiles/2.cmsearch.tblout.deoverlapped
Files are identical
Saved   377 hits (  644 removed) to testfiles/3.cmsearch.tblout.deoverlapped
Files are identical
Saved   357 hits (  614 removed) to testfiles/1.cmsearch.tblout.deoverlapped
Files are identical
Saved   325 hits (  596 removed) to testfiles/2.cmsearch.tblout.deoverlapped
Files are identical
Saved   377 hits (  644 removed) to testfiles/3.cmsearch.tblout.deoverlapped
Files are identical

This test indicates that cmsearch-deoverlap.pl removes overlapping
hits in an identical manner to cmscan version 1.1.2.

See details in:
/nfs/production/xfam/users/nawrocki/notebook/17_0508_inf_ebi_cmsearch_deoverlap/00LOG
