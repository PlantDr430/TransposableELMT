#!/usr/bin/python3 

"""
Written by Stephen A. Wyka (2019)

This is a batch wrapper that uses multiple repeat finding programs including RepeatModeler, 
TransposonPSI, LTR_finder, and LTR_harvest. LTR_harvest is coupled with LTR_digest and an 
HMMsearch against pfam domains associated with LTRs to limit false positive identifications. 
THe constructed libraries are run through RepeatClassifier to classify the LTR's. USEARCH is 
then used on the concatenated library to remove redundantLTR's based on a 80% similarity. 
The non-redundant library is then used with RepeatMasker to soft mask the assembly.

Additional curated libraries are recommended (such as RepBase libraries). Please make sure
to look at the -h /--help menu for all options. 

"""

import os, sys, re, argparse, inspect, shutil, subprocess
rundir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(rundir)
currentdir = os.getcwd()
sys.path.insert(0, parentdir)

class MyFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(
    usage='python3 %(prog)s [options] -in genome_assembly.fasta -o output_basename',
    description = '''    Wrapper script for de novo LTR identification and masking 
    of a genome assembly using RepeatModeler, LTR_Finder, LTR_Harvest + LTR_digest, 
    TransposonPSI, RepeatClassifier, and RepeatMasker.

    Additional curated libraries are recommended (such as RepBase libraries) and can 
    be passed through with the -rl flag''',
    
    epilog = """Written by Stephen A. Wyka (2019)""",
    formatter_class = MyFormatter)

parser.add_argument(
    '-in',
    '--input',
    required=True,
    help = 'Genome assembly in FASTA format',
    metavar=''
)
parser.add_argument(
    '-o',
    '--out',
    required=True,
    help = 'Basename of output directory and file',
    metavar=''
)
parser.add_argument(
    '--cpus',
    default=2,
    type=int,
    help = 'Number of cores to use [default: 2]',
    metavar=''
)
parser.add_argument(
    '-id',
    '--identity',
    default=0.80,
    type=float,
    help = 'Cutoff value for percent identity in USEARCH [default: 0.80]',
    metavar=''
)
parser.add_argument(
    '-en',
    '--engine',
    default='ncbi',
    choices = ['abblast','wublast','ncbi'],
    help = 'Search engine used in RepeatModeler [abblast|wublast|ncbi] [default: ncbi]',
    metavar=''
)
parser.add_argument(
    '-rb',
    '--repbase_lib',
    help = 'RepBase library of TEs or additional curated library in FASTA format',
    metavar=''
)
parser.add_argument(
    '-rl',
    '--repeatmodeler_lib',
    help = 'Pre-computed RepeatModeler library',
    metavar=''
)
parser.add_argument(
    '--hmms',
    default=os.path.join(rundir, 'te_hmms'),
    help = 'Path to directory of TE pfam domain files in HMMER3 format [Default: TransposableELMT/te_hmms]',
    metavar=''
)
parser.add_argument(
    '--REPEATMODELER_PATH',
    help = 'Path to RepeatModeler exe if not set in $PATH',
    metavar=''
)
parser.add_argument(
    '--REPEATMASKER_PATH',
    help = 'Path to RepeatMasker exe if not set in $PATH',
    metavar=''
)
parser.add_argument(
    '--BUILDDATABASE_PATH',
    help = 'Path to BuildDatabase exe if not set in $PATH',
    metavar=''
)
parser.add_argument(
    '--REPEATCLASSIFIER_PATH',
    help = 'Path to RepeatClassifier exe if not set in $PATH',
    metavar=''
)
parser.add_argument(
    '--LTRFINDER_PATH',
    help = 'Path to LTR_Finder exe if not set in $PATH',
    metavar=''
)
parser.add_argument(
    '--GENOMETOOLS_PATH',
    help = 'Path to genometools exe if not set in $PATH',
    metavar=''
)
parser.add_argument(
    '--USEARCH_PATH',
    help = 'Path to USEARCH exe if not set in $PATH',
    metavar=''
)
parser.add_argument(
    '--TRANSPOSONPSI_PATH',
    help = 'Path to transposonPSI.pl if not set in $PATH',
    metavar=''
)
parser.add_argument(
    '--CNV_LTRFINDER2GFF_PATH',
    help = 'Path to cnv_ltrfinder2gff.pl if not set in $PATH',
    metavar=''
)
args=parser.parse_args()

def which_path(file_name):
    for path in os.environ["PATH"].split(os.pathsep):
        full_path = os.path.join(path, file_name)
        if os.path.exists(full_path) and os.access(full_path, os.X_OK):
            return full_path
    return None

# modified os.walk() to only search into one level
def walklevel(path_dir, level=0):
    path_dir = path_dir.rstrip(os.path.sep)
    if os.path.isdir(path_dir):
        assert os.path.isdir(path_dir)
        num_sep = path_dir.count(os.path.sep)
        for root, dirs, files in os.walk(path_dir):
            yield root, dirs, files
            num_sep_this = root.count(os.path.sep)
            if num_sep + level <= num_sep_this:
                del dirs[:]
    else:
        None

# create output folders and define paths to folders
if not os.path.isdir(args.out):
    os.makedirs(os.path.join(args.out, 'working_directory'))
    os.makedirs(os.path.join(args.out, 'repeatmasker'))
result_dir = os.path.abspath(os.path.join(currentdir, args.out))
work_dir = os.path.abspath(os.path.join(args.out, 'working_directory'))
masker_dir = os.path.abspath(os.path.join(args.out, 'repeatmasker'))
if not os.path.isdir(work_dir):
    os.makedirs(os.path.join(work_dir, 'repeatmodeler'))
    os.makedirs(os.path.join(work_dir, 'ltr_finder'))
    os.makedirs(os.path.join(work_dir, 'ltr_harvest'))
    os.makedirs(os.path.join(work_dir, 'usearch'))
    os.makedirs(os.path.join(work_dir, 'TransposonPSI'))
    os.makedirs(os.path.join(work_dir, 'RepeatClassifier'))
else:
    # double check directories exists
    dirs = [os.path.join(args.out, 'working_directory'),os.path.join(args.out, 'repeatmasker'),
    os.path.join(work_dir, 'repeatmodeler'),os.path.join(work_dir, 'ltr_finder'),os.path.join(work_dir, 'ltr_harvest'),
    os.path.join(work_dir, 'usearch'),os.path.join(work_dir, 'TransposonPSI'),
    os.path.join(work_dir, 'RepeatClassifier')]
    for d in dirs:
        if not os.path.isdir(d):
            os.makedirs(d)
modeler_dir = os.path.abspath(os.path.join(work_dir, 'repeatmodeler'))
finder_dir = os.path.abspath(os.path.join(work_dir, 'ltr_finder'))
harvest_dir = os.path.abspath(os.path.join(work_dir, 'ltr_harvest'))
usearch_dir = os.path.abspath(os.path.join(work_dir, 'usearch'))
psi_dir = os.path.abspath(os.path.join(work_dir, 'TransposonPSI'))
class_dir = os.path.abspath(os.path.join(work_dir, 'RepeatClassifier'))

# Checking dependencies
try:
    if which_path('bedtools'):
        BEDTOOLS = 'bedtools'
    else:
        raise
except:
    print('bedtools not found, please make sure parent directory of', \
    'bedtools is located in $PATH')
    sys.exit(1)

try:
    if which_path('samtools'):
        SAMTOOLS = 'samtools'
    else:
        raise
except:
    print('samtools not found, please make sure parent directory of', \
    'bedtools is located in $PATH')
    sys.exit(1)

try:
    if which_path('perl'):
        None
    else:
        raise
except:
    print('perl not found, please make sure parent directory of ', \
    'perl is located in $PATH')
    sys.exit(1)

try:
    if which_path('hmmbuild'):
        None
    else:
        raise
except:
    print('hmmbuild not found, please make sure parent directory of', \
    'hmmbuild is located in $PATH')
    sys.exit(1)

try:
    if which_path('hmmsearch'):
        None
    else:
        raise
except:
    print('hmmsearch not found, please make sure parent directory of', \
    'hmmsearch is located in $PATH')
    sys.exit(1)

# Checking arguments
if args.input:
    input_file = os.path.abspath(os.path.join(work_dir, args.out+'_assembly.fasta'))
    shutil.copy(args.input, input_file)
else:
    print('Error: Please provide a genome assembly, -in or --input')
genome_file = input_file

if args.hmms:
    te_hmms = os.path.abspath(os.path.join(args.hmms, '*.hmm'))

if args.REPEATMODELER_PATH:
    REPEATMODELER = args.REPEATMODELER_PATH
else:
    try:
        if which_path('RepeatModeler'):
            REPEATMODELER = 'RepeatModeler'
        else:
            raise
    except:
        print('ERROR: RepeatModeler not found, please make sure parent directory of', \
        'RepeatModeler is located in $PATH or provide path to executable in', \
        'command with --REPEATMODELER_PATH')
        sys.exit(1)

if args.BUILDDATABASE_PATH:
    BUILDDATABASE = args.BUILDDATABASE_PATH
else:
    try:
        if which_path('BuildDatabase'):
            BUILDDATABASE = 'BuildDatabase'
        else:
            raise
    except:
        print('ERROR: BuildDatabase not found, please make sure parent directory of', \
        'BuildDatabase is located in $PATH or provide path to executable in', \
        'command with --BUILDDATABASE_PATH')
        sys.exit(1)

if args.LTRFINDER_PATH:
    LTRFINDER = args.LTRFINDER_PATH
else:
    try:
        if which_path('ltr_finder'):
            LTRFINDER = 'ltr_finder'
        else:
            raise
    except:
        print('ERROR: ltr_finder not found, please make sure parent directory of', \
        'ltr_finder is located in $PATH or provide path to executable in', \
        'command with --LTRFINDER_PATH')
        sys.exit(1)

if args.CNV_LTRFINDER2GFF_PATH:
    FINDER2GFF = args.CNV_LTRFINDER2GFF_PATH
else:
    try:
        if which_path('cnv_ltrfinder2gff.pl'):
            FINDER2GFF = 'cnv_ltrfinder2gff.pl'
        else:
            raise
    except:
        print('ERROR: cnv_ltrfinder2gff.pl not found, please make sure parent directory of', \
        'cnv_ltrfinder2gff.pl is located in $PATH or provide path to executable in', \
        'command with --CNV_LTRFINDER2GFF_PATH')
        sys.exit(1)

if args.GENOMETOOLS_PATH:
    GENOMETOOLS = args.GENOMETOOLS_PATH
else:
    try:
        if which_path('gt'):
            GENOMETOOLS = 'gt'
        else:
            raise
    except:
        print('ERROR: genometools (gt) not found, please make sure parent directory of', \
        'genometools (gt) is located in $PATH or provide path to executable in', \
        'command with --GENOMETOOLS_PATH')
        sys.exit(1)

if args.TRANSPOSONPSI_PATH:
    TRANSPOSONPSI = args.TRANSPOSONPSI_PATH
else:
    try:
        if which_path('transposonPSI.pl'):
            TRANSPOSONPSI = 'transposonPSI.pl'
        else:
            raise
    except:
        print('ERROR: transposonPSI.pl not found, please make sure parent directory of', \
        'transposonPSI.pl is located in $PATH or provide path to executable in', \
        'command with --TRANSPOSONPSI_PATH')
        sys.exit(1)

if args.USEARCH_PATH:
    USEARCH = args.USEARCH_PATH
try:
    if which_path('usearch'):
        USEARCH = 'usearch'
    else: ## Trying to find USEARCH dependency as naming of executable can be strange
        for path in os.environ['PATH'].split(os.pathsep):
            path_dir = path
            attempts = []
            for files in walklevel(path_dir):
                dir_list = list(files)
                for x in dir_list:
                    for y in x:
                        if 'usearch' in y and 'linux' in y:
                            attempts.append(os.path.join(path,y))
                            if len(attempts) == 1:
                                for program in attempts:
                                    if os.access(program, os.X_OK):
                                        USEARCH = program
                            else:
                                raise
except:
    print('ERROR: Difficulty finding usearch executable, please make sure parent directory of', \
    'usearch is located in $PATH or provide path to executable in', \
    'command with --USEARCH_PATH')
    sys.exit(1)

if args.REPEATCLASSIFIER_PATH:
    REPEATCLASSIFIER = args.REPEATCLASSIFIER_PATH
else:
    try:
        if which_path('RepeatClassifier'):
            REPEATCLASSIFIER = 'RepeatClassifier'
        else:
            raise
    except:
        print('ERROR: RepeatClassifier not found, please make sure parent directory of', \
        'RepeatClassifier is located in $PATH or provide path to executable in', \
        'command with --REPEATCLASSIFIER_PATH')
        sys.exit(1)

if args.REPEATMASKER_PATH:
    REPEATMASKER = args.REPEATMASKER_PATH
else:
    try:
        if which_path('RepeatMasker'):
            REPEATMASKER = 'RepeatMasker'
        else:
            raise
    except:
        print('ERROR: RepeatMasker not found, please make sure parent directory of', \
        'RepeatMasker is located in $PATH or provide path to executable in', \
        'command with --REPEATMASKER_PATH')
        sys.exit(1)

#### START RUN ####

## RepeatModeler run ##

if not args.repeatmodeler_lib: #no pre-computed repeat library is given, so must run RepeatModeler
    print('Running RepeatModeler: This will take a while...')
    modeler_log = os.path.abspath(os.path.join(modeler_dir, 'repeatmodeler.log'))
    if os.path.exists(modeler_log):
        os.remove(modeler_log)
    for x in os.listdir(modeler_dir): # remove old folder if present
            if x.startswith('RM_'):
                shutil.rmtree(os.path.join(modeler_dir,x))
    with open(modeler_log, 'a') as rm_log:
        subprocess.call([BUILDDATABASE, '-name', args.out, '-engine', args.engine, '-dir', work_dir
        ], cwd=modeler_dir, stdout=rm_log, stderr=rm_log)
        subprocess.call([REPEATMODELER, '-e', args.engine , '-database', args.out, '-pa', str(args.cpus)
        ], cwd=modeler_dir, stdout=rm_log, stderr=rm_log)
        for x in os.listdir(modeler_dir):
            if x.startswith('RM_'):
                modeler_out = os.path.abspath(os.path.join(modeler_dir,x))
                modeler_output = os.path.join(modeler_out, 'consensi.fa.classified')
                shutil.copy(modeler_output, os.path.join(work_dir, 'modeler_library.fasta'))

## ltr_finder run ##

print("Running ltr_finder: This won't take long")
finder_output = os.path.abspath(os.path.join(finder_dir, 'ltr_finder.txt'))
finder_gff = os.path.abspath(os.path.join(finder_dir, args.out+'_ltrfinder.gff'))
finder_gff3 = os.path.abspath(os.path.join(finder_dir, args.out+'_ltrfinder.gff3'))
finder_results = os.path.abspath(os.path.join(finder_dir, 'ltrfinder_library.fasta'))
if os.path.exists(finder_output):
    os.remove(finder_output)
with open(finder_output, 'a') as lf_out:
    subprocess.call([LTRFINDER, genome_file], cwd=finder_dir, stdout=lf_out, stderr=lf_out)
    subprocess.call([FINDER2GFF, '-i', finder_output, '-o', args.out+'_ltrfinder.gff', '--gff-ver', 'GFF3'
    ], cwd=finder_dir, stderr=lf_out)
with open(finder_gff, 'r') as gff:
    keep_line = []
    for line in gff:
        if 'LTR_retrotransposon' in line:
            keep_line.append(line)
with open(finder_gff3, 'w') as gff3:
    gff3.write(''.join(keep_line))
subprocess.call([BEDTOOLS, 'getfasta', '-name', '-fi', genome_file, '-bed', finder_gff3, '-fo', 
'ltrfinder_library.fasta'], cwd=finder_dir)
os.remove(genome_file+'.fai')
shutil.copy(finder_results, work_dir)


## ltr_harvest run ##

print("Running ltr_harvest: This shouldn't take long")
harvest_log = os.path.abspath(os.path.join(harvest_dir, 'ltr_harvest.log'))
index_dir = os.path.abspath(os.path.join(harvest_dir, 'index'))
index_name = os.path.abspath(os.path.join(index_dir, args.out))
sorted_gff = os.path.abspath(os.path.join(harvest_dir, args.out+'_sorted_harvest.gff3'))
digest_gff3 = os.path.abspath(os.path.join(harvest_dir, args.out+'_filter_digest.gff3'))
harvest_results = os.path.abspath(os.path.join(harvest_dir,'harvest_library.fasta'))
if os.path.exists(harvest_results):
    os.remove(harvest_results)
if os.path.exists(harvest_log):
    os.remove(harvest_log)
if os.path.exists(index_dir):
    shutil.rmtree(index_dir)
os.makedirs(os.path.join(harvest_dir, 'index'))
with open(harvest_log, 'a') as h_log:
    subprocess.call([GENOMETOOLS, 'suffixerator', '-db', genome_file, '-indexname', os.path.join(index_dir, args.out),
    '-tis', '-suf', '-lcp', '-des', '-ssp', '-sds', '-dna'], cwd=index_dir, stdout=h_log, stderr=h_log)
    subprocess.call([GENOMETOOLS, 'ltrharvest', '-index', index_name, '-out', args.out+'_harvest.fasta',
    '-gff3', args.out+'_harvest.gff3'], cwd=harvest_dir, stdout=h_log, stderr=h_log)
with open(sorted_gff, 'w') as sort_gff, open(harvest_log, 'a') as h_log:
    subprocess.call([GENOMETOOLS, 'gff3', '-sort', args.out+'_harvest.gff3'
    ], cwd=harvest_dir, stdout=sort_gff, stderr=h_log)
os.chdir(harvest_dir)
digest_run = [GENOMETOOLS, 'ltrdigest', '-hmms', te_hmms, '-outfileprefix', args.out, sorted_gff, 
index_name, '>', args.out+'_all_digest.gff3', '2>', 'ltr_digest_errors.log']
os.system(' '.join(digest_run)) #had to use os.system call due to difficulties with subprocess.call
with open(args.out+'_all_digest.gff3', 'r') as gff, open(args.out+'_filter_digest.gff3', 'w') as gff3:
    seqs = []
    contigs = []
    models = []
    all_models = []
    for line in gff:
        if line.startswith('##s'):
            seq_id = re.search(r'(seq\d+)', line).group(1)
            seqs.append(seq_id)
        elif line.startswith('#c'):
            contig_id = re.search(r'(contig_\d+)', line).group(1)
            contigs.append(contig_id)
        elif line.startswith('seq'):
            models.append(line)
        elif line.startswith('###'):
            all_models.append(models)
            models = []
    
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    sorted_seqs = sorted(seqs, key = alphanum_key)
    
    combined = dict(zip(sorted_seqs, contigs))
    
    filter_models = []
    final_models = []
    LTR_models = []
    for feat in all_models:
        for x in feat:
            if 'protein_match' in x:
                filter_models.append(feat)
    for feat in filter_models:
        for x in feat:
            seq_id = re.search(r'(seq\d+)', x).group(1)
            if seq_id in combined.keys():
                y = x.replace(seq_id, combined[seq_id])
                final_models.append(y)
    for feat in final_models:
        if '\trepeat_region\t' in feat:
            LTR_models.append(feat)
    gff3.write(''.join(set(LTR_models)).replace('\trepeat_region\t', '\tLTR_retrotransposon_harvest\t'))
with open(harvest_log, 'a') as h_log:
    subprocess.call([BEDTOOLS, 'getfasta', '-name', '-fi', genome_file, '-bed', digest_gff3, '-fo', 
    'harvest_library.fasta'], cwd=harvest_dir, stderr=h_log)
os.remove(genome_file+'.fai')
shutil.copy(harvest_results, work_dir)

sys.exit()

## TransposonPSI run ##

print("Running TransposonPSI: This will take some time")
psi_log = os.path.abspath(os.path.join(psi_dir, 'TransposonPSI.log'))
psi_gff3 = os.path.abspath(os.path.join(psi_dir, args.out+'_Tpsi.gff3'))
psi_results = os.path.abspath(os.path.join(psi_dir,'Tpsi_library.fasta'))
if os.path.exists(psi_results):
    os.remove(psi_results)
if os.path.exists(psi_log):
    os.remove(psi_log)
try:
    with open(psi_log, 'a') as psi_log:
        subprocess.call([TRANSPOSONPSI, genome_file, 'nuc'], cwd=psi_dir, stdout=psi_log, stderr=psi_log)
except:
    print('There was an error with TransposonPSI, please check the logfile located in {}'.format(psi_dir))
best_hits = [f for f in os.listdir(psi_dir) if os.path.isfile(os.path.join(psi_dir, f)) and 'bestPerLocus.gff3' in f]
best_hits_gff = os.path.abspath(os.path.join(psi_dir, ''.join(best_hits)))
with open(best_hits_gff, 'r') as best_hits, open(psi_gff3, 'w') as gff3:
    for line in best_hits:
        col = line.split('\t')
        target = re.search(r'Target=(\w+)', col[8]).group(1)
        gff3.write(col[0] + '\t' + col[1] + '\t' + target + '\t' + col[3] +
        '\t' + col[4] + '\t'+ col[5] + '\t' + col[6] + '\t' + col[7] + '\t' + col[8])
subprocess.call([BEDTOOLS, 'getfasta', '-name', '-fi', genome_file, '-bed', psi_gff3, '-fo', 
'Tpsi_library.fasta'], cwd=psi_dir, stderr=subprocess.DEVNULL)
os.remove(genome_file+'.fai')
shutil.copy(psi_results, work_dir)

## Preparing files for RepeatClassifier ##

ltr_finder_lib = os.path.abspath(os.path.join(work_dir, 'ltrfinder_library.fasta'))
ltr_harvest_lib = os.path.abspath(os.path.join(work_dir, 'harvest_library.fasta'))
Tpsi_lib = os.path.abspath(os.path.join(work_dir, 'Tpsi_library.fasta'))
concatenated_lib = os.path.abspath(os.path.join(class_dir, 'concat_library.fasta'))
if os.path.exists(concatenated_lib):
    os.remove(concatenated_lib)
with open(concatenated_lib, 'w') as concat_lib:
    subprocess.call(['cat', ltr_finder_lib, ltr_harvest_lib, Tpsi_lib
    ], cwd=class_dir, stdout=concat_lib)

## RepeatClassifier run ##

print("Running RepeatClassifier: This can take some time depending on the size of the library")
unclassified_lib = os.path.abspath(os.path.join(class_dir, 'concat_library.fasta'))
class_log = os.path.abspath(os.path.join(class_dir, 'repeatclassifier.log'))
classified_lib = os.path.abspath(os.path.join(class_dir, 'concat_library.fasta.classified'))
classified_mask = os.path.abspath(os.path.join(class_dir, 'concat_library.fasta.masked'))
if os.path.exists(class_log):
    os.remove(class_log)
with open(class_log, 'a') as cl_log:
    subprocess.call([REPEATCLASSIFIER, '-consensi', unclassified_lib, '-engine', args.engine
    ], cwd=class_dir, stdout=cl_log, stderr=cl_log)

# Preparing files for USEARCH ##

usearch_prep = os.path.abspath(os.path.join(usearch_dir, 'usearch_prep.fasta'))
if args.repbase_lib:
    repbase_lib = os.path.abspath(os.path.join(currentdir, args.repbase_lib))
if args.repeatmodeler_lib:
    repeatmodeler_lib = os.path.abspath(os.path.join(currentdir, args.repeatmodeler_lib))
else:
    repeatmodeler_lib = os.path.abspath(os.path.join(work_dir, 'modeler_library.fasta'))
with open(usearch_prep, 'w') as u_prep:
    subprocess.call(['cat', repeatmodeler_lib, repbase_lib, classified_lib], cwd=usearch_dir, stdout=u_prep)

## USEARCH run ##

print("Running Usearch: This won't take too long.")
usearch_results = os.path.abspath(os.path.join(usearch_dir, 'nonredundant_library.fasta'))
usearch_log = os.path.abspath(os.path.join(usearch_dir, 'usearch.log'))
if os.path.exists(usearch_results):
    os.remove(usearch_results)
if os.path.exists(usearch_log):
    os.remove(usearch_log)
with open(usearch_log, 'a') as u_log:
    subprocess.call([USEARCH, '-cluster_fast', usearch_prep, '-id', str(args.identity), '-threads', 
    str(args.cpus), '-centroids', usearch_results, '-sort', 'length'
    ], cwd=usearch_dir, stdout=subprocess.DEVNULL, stderr=u_log)

## Fixing up some Unknown calls, based on RepBase calls. This is not meant to be extensive ##

final_lib = os.path.abspath(os.path.join(work_dir, 'final_library.fasta'))
LTR_list = ['Copia', 'copia', 'gypsy', 'Gypsy', 'LTR_retrotransposon_harvest', 'CALTR', 'CEN1_SP', 'Roo']
DNA_list = ['mariner', 'cacta', 'hAT', 'EnSpm', 'DNA transposon', 'MuDR', 'Academ', 'Dada', 'Crypton']
RC_list = ['Helitron', 'helitron']
LINE_list = ['Non-LTR Retrotransposon']
SINE_list = ['FOXY']
SIMPLE_list = ['Simple Repeat']
ALL_list = LTR_list + DNA_list + RC_list + LINE_list + SINE_list + SIMPLE_list

with open(usearch_results, 'r') as u_res, open(final_lib, 'w') as final:
    for line in u_res:
        if line.startswith('>') and '#Unknown' in line and any(x in line for x in ALL_list):
            for x in LTR_list:
                if x == 'LTR_retrotransposon_harvest' and x in line:
                    final.write(line.replace('#Unknown', '#LTR/'))
                elif x in line:
                    final.write(line.replace('#Unknown', '#LTR/'+x))
            for x in DNA_list:
                if x == 'DNA transposon' and x in line:
                    final.write(line.replace('#Unknown', '#DNA/'))
                elif x in line:
                    final.write(line.replace('#Unknown', '#DNA/'+x))
            for x in RC_list:
                if x in line:
                    final.write(line.replace('#Unknown', '#RC/'+x))
            for x in LINE_list:
                if x in line:
                    final.write(line.replace('#Unknown', '#LINE'))
            for x in SINE_list:
                if x in line:
                    final.write(line.replace('#Unknown', '#SINE'))
            for x in SIMPLE_list:
                if x in line:
                    final.write(line.replace('#Unknown', '#Simple_repeat'))
        else:
            final.write(line)

## Run RepeatMasker ##

print("Running RepeatMasker: In the home stretch now")
masker_log = os.path.abspath(os.path.join(masker_dir, 'repeatmasker.log'))
masked_genome = os.path.abspath(os.path.join(masker_dir, args.out+'_assembly.fasta.masked'))
repeat_table = os.path.abspath(os.path.join(masker_dir, args.out+'_assembly.fasta.tbl'))
repeat_output = os.path.abspath(os.path.join(masker_dir, args.out+'_assembly.fasta.out'))
if os.path.exists(masker_log):
    os.remove(masker_log)
with open(masker_log, 'a') as mask_log:
    subprocess.call([REPEATMASKER, '-pa', str(args.cpus), '-gff', '-xm', '-lcambig', '-cutoff', '300', '-xsmall', 
    '-gccalc', '-dir', masker_dir, '-lib', final_lib, genome_file], cwd=masker_dir, stdout=mask_log, stderr=mask_log)
shutil.copy(masked_genome, os.path.join(result_dir, args.out+'_masked.fasta'))
shutil.copy(repeat_table, os.path.join(result_dir, args.out+'_repeatmasker.tbl'))
shutil.copy(repeat_output, os.path.join(result_dir, args.out+'_repeatmasker.output'))

print("Job done! Please find your masked genome, repeatmasker table, and repeatmasker output", \
"in the directory {}/".format(args.out))
