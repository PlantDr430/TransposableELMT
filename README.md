# TransposableELMT 
#### Wrapper script for TE identification and genome masking

## Summary

This is a batch wrapper that uses multiple repeat finding programs including RepeatModeler, 
TransposonPSI, LTR_finder, and LTR_harvest. LTR_harvest is coupled with LTR_digest and an 
HMMsearch against pfam domains associated with LTRs to limit false positive identifications. 
THe constructed libraries are run through RepeatClassifier to classify the LTR's. USEARCH is 
then used on the concatenated library to remove redundantLTR's based on a 80% similarity. 
The non-redundant library is then used with RepeatMasker to soft mask the assembly.

Currently, all programs are run using default settings with little to no options to alter settings through flags. Additional options may be added to future versions if there is a need.

It is recommended to provide additional currated libraries such as those from [RepBase](https://www.girinst.org/repbase/update/browse.php). Simply select an appropriate toxanomic level and download the file in FASTA format. Then provide the file with the ```-rb``` flag on the command line.

## Dependencies 

#### Basic programs
1. [Python 3](https://www.python.org/downloads/)
2. [bedtools](https://github.com/arq5x/bedtools2)
3. [samtools](https://github.com/samtools/samtools)
4. [perl](https://www.perl.org/get.html)
5. [HMMER](http://hmmer.org/)

#### TE programs
1. [RepeatModeler + RepeatClassifer + BuildDatabase](https://github.com/rmhubley/RepeatModeler)
3. [RepeatMasker](https://github.com/rmhubley/RepeatMasker)
4. [LTR_Finder](https://github.com/xzhub/LTR_Finder)
5. [Genometools](https://github.com/genometools/genometools)
6. [TransposonPSI](http://transposonpsi.sourceforge.net/)

#### Additional
1. [USEARCH](https://www.drive5.com/usearch/download.html)
3. [cnv_ltrfinder2gff.pl](https://github.com/jestill/dawgpaws/blob/master/scripts/cnv_ltrfinder2gff.pl)

Dependecies should be able to be called from the commandline, if not then the paths to the parent directories of each executable should be located in $PATH. If all else fails, paths to executables can be passed into the script throguh flags.

## Usage

```
usage: ./new_ltr_wrapper.py [options] -in genome_assembly.fasta -o output_basename

optional arguments:
  -h, --help                  show this help message and exit
  -in , --input               Genome assembly in FASTA format
  -o , --out                  Basename of output directory and file
  --cpus                      Number of cores to use [default: 2]
  -id , --identity            Cutoff value for percent identity in USEARCH [default: 0.80]
  -en , --engine              Search engine used in RepeatModeler [abblast|wublast|ncbi] [default: ncbi]
  -rb , --repbase_lib         RepBase library of TEs or additional curated library in FASTA format
  -rl , --repeatmodeler_lib   Pre-computed RepeatModeler library
  --hmms                      Path to directory of TE pfam domain files in HMMER3 format [Default: TransposableELMT/te_hmms]
  --REPEATMODELER_PATH        Path to RepeatModeler exe if not set in $PATH
  --REPEATMASKER_PATH         Path to RepeatMasker exe if not set in $PATH
  --BUILDDATABASE_PATH        Path to BuildDatabase exe if not set in $PATH
  --REPEATCLASSIFIER_PATH     Path to RepeatClassifier exe if not set in $PATH
  --LTRFINDER_PATH            Path to LTR_Finder exe if not set in $PATH
  --GENOMETOOLS_PATH          Path to genometools exe if not set in $PATH
  --USEARCH_PATH              Path to USEARCH exe if not set in $PATH
  --TRANSPOSONPSI_PATH        Path to transposonPSI.pl if not set in $PATH
  --CNV_LTRFINDER2GFF_PATH    Path to cnv_ltrfinder2gff.pl if not set in $PATH
```

## Output files

1. Soft-masked genome assembly in FASTA format
2. RepeatMasker Table file
3. RepeatMasker Out file
