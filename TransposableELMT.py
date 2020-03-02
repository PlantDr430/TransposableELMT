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
from collections import defaultdict
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
rundir = os.getcwd()
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
    default=os.path.join(currentdir, 'te_hmms'),
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
args=parser.parse_args()


# Preliminary functions

def which_path(file_name):
    '''Function for checking executable paths'''
    for path in os.environ["PATH"].split(os.pathsep):
        full_path = os.path.join(path, file_name)
        if os.path.exists(full_path) and os.access(full_path, os.X_OK):
            return full_path
    return None

def walklevel(path_dir, level=0):
    '''Modified os.walk() to only search into one level'''
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

# create folders and define paths to folders
if not os.path.isdir(args.out):
    os.makedirs(os.path.join(args.out, 'working_directory'))
    os.makedirs(os.path.join(args.out, 'repeatmasker'))
result_dir = os.path.abspath(os.path.join(rundir, args.out))
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
    os.path.join(work_dir, 'repeatmodeler'),os.path.join(work_dir, 'ltr_finder'),
    os.path.join(work_dir, 'ltr_harvest'),os.path.join(work_dir, 'usearch'),
    os.path.join(work_dir, 'TransposonPSI'),
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

# Checking dependencies and arguments
try:
    if which_path('bedtools'):
        BEDTOOLS = 'bedtools'
    else:
        raise
except:
    print('bedtools not found, please make sure parent directory of', \
    'bedtools is located in $PATH')
    sys.exit()

try:
    if which_path('samtools'):
        SAMTOOLS = 'samtools'
    else:
        raise
except:
    print('samtools not found, please make sure parent directory of', \
    'bedtools is located in $PATH')
    sys.exit()

try:
    if which_path('perl'):
        None
    else:
        raise
except:
    print('perl not found, please make sure parent directory of ', \
    'perl is located in $PATH')
    sys.exit()

try:
    if which_path('hmmbuild'):
        None
    else:
        raise
except:
    print('hmmbuild not found, please make sure parent directory of', \
    'hmmbuild is located in $PATH')
    sys.exit()

try:
    if which_path('hmmsearch'):
        None
    else:
        raise
except:
    print('hmmsearch not found, please make sure parent directory of', \
    'hmmsearch is located in $PATH')
    sys.exit()

if args.input:
    input_file = os.path.abspath(os.path.join(work_dir, args.out+'_assembly.fasta'))
    modeler_copy = os.path.abspath(os.path.join(modeler_dir, args.out+'_assembly.fasta'))
    shutil.copy(args.input, input_file)
    shutil.copy(args.input, modeler_copy)
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
        sys.exit()

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
        sys.exit()

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
        sys.exit()

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
        sys.exit()

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
        sys.exit()

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
    sys.exit()

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
        sys.exit()

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
        sys.exit()

######## START RUN ########

def run_repeatmodeler():
    print('Running RepeatModeler: This will take a while...')
    modeler_log = os.path.abspath(os.path.join(modeler_dir, 'repeatmodeler.log'))
    if os.path.exists(modeler_log):
        os.remove(modeler_log)
    for x in os.listdir(modeler_dir): # remove old folder if present
        if x.startswith('RM_'):
            shutil.rmtree(os.path.join(modeler_dir,x))
    with open(modeler_log, 'a') as rm_log:
        try:
            subprocess.call([BUILDDATABASE, '-name', args.out, '-engine', args.engine, '-dir', modeler_dir
            ], cwd=modeler_dir, stdout=rm_log, stderr=rm_log)
        except:
            print('ERROR: There was an error building a database for RepeatModeler, please check log \
            file located at {}'.format(modeler_log))
            sys.exit()
        try:
            subprocess.call([REPEATMODELER, '-e', args.engine , '-database', args.out, '-pa', str(args.cpus)
            ], cwd=modeler_dir, stdout=rm_log, stderr=rm_log)
            for x in os.listdir(modeler_dir):
                if x.startswith('RM_'):
                    modeler_out = os.path.abspath(os.path.join(modeler_dir,x))
                    modeler_output = os.path.join(modeler_out, 'consensi.fa.classified')
                    shutil.copy(modeler_output, os.path.join(work_dir, 'modeler_library.fasta'))
        except:
            print('ERROR: There was an error running RepeatModeler, please check log \
            file located at {}'.format(modeler_log))
            sys.exit()

def parse_ltr_finder_2gff(fpath, opath):
    '''Parse the LTR_finder output into gff format.'''
    line_breaks = defaultdict(list) # {<contig_ID>-<subset_number> : [Nested list of lines from LTR_output]}
    with open(fpath, 'r') as f_in:
        contig_sub = ''
        for line in f_in:
            match_start = re.match(r'(\[\d+\])\s(.+)\s(.+)', line)
            if match_start != None:
                contig = match_start.group(2)
                sub = ' '.join(re.findall(r'\[(\d+)\]',match_start.group(1)))
                contig_sub = contig+'-'+sub
                continue
            line_breaks[contig_sub].append([line])
    line_breaks.pop('')

    gff_list = []
    for match, lines in line_breaks.items():
        contig = match.split('-')[0]
        subset = match.split('-')[1]
        prog = 'ltr_finder'
        score = ' '.join(re.findall(r'\:\s(\d+)\s\[',lines[1][0]))
        phase = '.'
        strand = re.search(r'(Strand:)(.)',lines[0][0]).group(2)
        # id = contig+'_'+prog+'_model0001'
        parent = contig+'_'+prog+'_model000{}'.format(subset)
        for l in lines:
            l = l[0] # turn list into string
            if l.startswith('Location'):
                sect = re.search(r'(:)\s(\d+)\s-\s(\d+)', l)
                start = sect.group(2)
                stop = sect.group(3)
                attribute = 'ID={}'.format(parent)
                gff_line = [contig, prog, 'LTR_retrotransposon', start, stop, score, strand, phase, attribute]
                gff_list.append(gff_line)
            elif l.startswith("5'-LTR"):
                sect = re.search(r'(:)\s(\d+)\s-\s(\d+)', l)
                start = sect.group(2)
                stop = sect.group(3)
                attribute = 'ID={};Name=Five Prime LTR;Parent={}'.format(parent+'_five_prime_LTR',parent)
                gff_line = [contig, prog, 'five_prime_LTR', start, stop, score, strand, phase, attribute]
                gff_list.append(gff_line)
            elif l.startswith("3'-LTR"):
                sect = re.search(r'(:)\s(\d+)\s-\s(\d+)', l)
                start = sect.group(2)
                stop = sect.group(3)
                attribute = 'ID={};Name=Three Prime LTR;Parent={}'.format(parent+'_three_prime_LTR',parent)
                gff_line = [contig, prog, 'three_prime_LTR', start, stop, score, strand, phase, attribute]
                gff_list.append(gff_line)
            elif l.startswith('TSR'):
                if not 'NOT FOUND' in l: # if target site duplications are found
                    sect = re.search(r'(:)\s(\d+)\s-\s(\d+)\s(,)\s(\d+)\s-\s(\d+)', l)
                    start1 = sect.group(2)
                    stop1 = sect.group(3)
                    start2 = sect.group(5)
                    stop2 = sect.group(6)
                    attribute = 'ID={};Name=Target Site Duplication;Parent={}'.format(parent+'_tsd5',parent)
                    gff_line = [contig, prog, 'target_site_duplication', start1, stop1, score, strand, phase, attribute]
                    gff_list.append(gff_line)
                    attribute = 'ID={};Name=Target Site Duplication;Parent={}'.format(parent+'_tsd3',parent)
                    gff_line = [contig, prog, 'target_site_duplication', start2, stop2, score, strand, phase, attribute]
                    gff_list.append(gff_line)
                else:
                    pass
            elif l.startswith('PPT'):
                sect = re.search(r'(:)\s(\[.+/.+\])\s(\d+)\s-\s(\d+)', l)
                start = sect.group(3)
                stop = sect.group(4)
                attribute = 'ID={};Name=PPT;Parent={}'.format(parent+'_RR_tract',parent)
                gff_line = [contig, prog, 'RR_tract', start, stop, score, strand, phase, attribute]
                gff_list.append(gff_line)
            elif l.startswith('Domain'):
                sect = re.search(r'(:)\s(\d+)\s-\s(\d+)', l)
                start = sect.group(2)
                stop = sect.group(3)
                attribute = 'ID={};Name=Reverse Transcriptase;Parent={}'.format(parent+'_rvt',parent)
                gff_line = [contig, prog, 'transposable_element_gene', start, stop, score, strand, phase, attribute]
                gff_list.append(gff_line)
            else:
                pass

    with open(opath, 'w') as f_out:
        f_out.write('##gff-version 3\n')
        for line in gff_list:
            f_out.write('\t'.join(line) + '\n')

def run_ltr_finder():
    print("Running ltr_finder: This won't take long")
    finder_output = os.path.abspath(os.path.join(finder_dir, 'ltr_finder.txt'))
    finder_log = os.path.abspath(os.path.join(finder_dir, 'ltr_finder.log'))
    finder_gff = os.path.abspath(os.path.join(finder_dir, args.out+'_ltrfinder.gff'))
    finder_gff3 = os.path.abspath(os.path.join(finder_dir, args.out+'_ltrfinder.gff3'))
    finder_results = os.path.abspath(os.path.join(finder_dir, 'ltrfinder_library.fasta'))
    if os.path.exists(finder_log):
        os.remove(finder_log)
    with open(finder_output, 'w') as lf_out, open(finder_log, 'w') as lf_log:
        try:
            subprocess.call([LTRFINDER, genome_file], cwd=finder_dir, stdout=lf_out, stderr=lf_log)
        except:
            print('ERROR: There was an error running LTR_finder, please check the run \
            log located at {}'.format(finder_log))
            sys.exit()
        parse_ltr_finder_2gff(finder_output, finder_gff)
    with open(finder_gff, 'r') as gff: # just keep full length lines (i.e. gene coordinates)
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

def run_ltr_harvest_and_digest():
    print("Running ltr_harvest: This shouldn't take long")
    harvest_log = os.path.abspath(os.path.join(harvest_dir, 'ltr_harvest.log'))
    index_dir = os.path.abspath(os.path.join(harvest_dir, 'index'))
    index_name = os.path.abspath(os.path.join(index_dir, args.out))
    sorted_gff = os.path.abspath(os.path.join(harvest_dir, args.out+'_sorted_harvest.gff3'))
    digest_gff3 = os.path.abspath(os.path.join(harvest_dir, args.out+'_filter_digest.gff3'))
    harvest_results = os.path.abspath(os.path.join(harvest_dir,'harvest_library.fasta'))
    if os.path.exists(harvest_log):
        os.remove(harvest_log)
    os.makedirs(os.path.join(harvest_dir, 'index'))
    with open(harvest_log, 'w') as h_log:
        try:
            subprocess.call([GENOMETOOLS, 'suffixerator', '-db', genome_file, '-indexname', 
            os.path.join(index_dir, args.out),'-tis', '-suf', '-lcp', '-des', '-ssp', '-sds', '-dna'
            ], cwd=index_dir, stdout=h_log, stderr=h_log)
        except:
            print('ERROR: There was an error running suffixerator of genome tools, please check log file \
            located at {}'.format(harvest_log))
            sys.exit()
        try:
            subprocess.call([GENOMETOOLS, 'ltrharvest', '-index', index_name, '-out', args.out+'_harvest.fasta',
            '-gff3', args.out+'_harvest.gff3'], cwd=harvest_dir, stdout=h_log, stderr=h_log)
        except:
            print('ERROR: There was an error running ltr_harvest of genome tools, please check log file \
            located at {}'.format(harvest_log))
            sys.exit()
    with open(sorted_gff, 'w') as sort_gff, open(harvest_log, 'a') as h_log:
        try:
            subprocess.call([GENOMETOOLS, 'gff3', '-sort', args.out+'_harvest.gff3'
            ], cwd=harvest_dir, stdout=sort_gff, stderr=h_log)
        except:
            print('ERROR: There was an error running gff3 (sort) of genome tools, please check log file \
            located at {}'.format(harvest_log))
            sys.exit()
    os.chdir(harvest_dir)
    digest_log = os.path.abspath(os.path.join(harvest_dir, 'ltr_digest_errors.log'))
    digest_gff = os.path.abspath(os.path.join(harvest_dir, args.out+'_all_digest.gff3'))
    digest_filter_gff = os.path.abspath(os.path.join(harvest_dir, args.out+'_filter_digest.gff3'))
    digest_run = [GENOMETOOLS, 'ltrdigest', '-hmms', te_hmms, '-outfileprefix', args.out, sorted_gff, 
    index_name, '>', args.out+'_all_digest.gff3', '2>', 'ltr_digest_errors.log']
    try:
        os.system(' '.join(digest_run)) #had to use os.system call due to difficulties with subprocess.call
    except:
        print('ERROR: There was an error running LTR digest of genome tools, please check log file \
        located at {}'.format(digest_log))
        sys.exit()
    with open(digest_gff, 'r') as gff, open(digest_filter_gff, 'w') as gff3:
        seqs = []
        contigs = []
        models = []
        all_models = []
        for line in gff:
            if line.startswith('##sequence-region'):
                seq_id = re.search(r'(seq\d+)', line).group(1)
                seqs.append(seq_id)
            elif line.startswith('seq'):
                models.append(line)
            elif line.startswith('###'):
                all_models.append(models)
                models = [] # reset model list to nest features together
            elif line.startswith('##gff'):
                pass
            elif line.startswith('#'): # lines with actual contig names
                contig_id = re.search(r'(#)(.*)', line).group(2)
                if ' ' in contig_id:
                    contigs.append(contig_id.split(' ')[0])
                else:
                    contigs.append(contig_id)

        # match up sequence-regions to correct contigs. Contigs are sorted correctly (1 -> last)
        # seq-regions are sorted seq1, seq10, seq10x, etc. BUT seq1 = 1st contig, seq2 = 2nd contig (not seq10)
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
        try:
            subprocess.call([BEDTOOLS, 'getfasta', '-name', '-fi', genome_file, '-bed', digest_gff3, '-fo', 
            'harvest_library.fasta'], cwd=harvest_dir, stderr=h_log)
        except:
            print('ERROR: There was an error running gff3 getfasta of bedtools, please check log file \
            located at {}'.format(harvest_log))
            sys.exit()
    os.remove(genome_file+'.fai')
    shutil.copy(harvest_results, work_dir)

def run_transposonPSI():
    print("Running TransposonPSI: This will take some time")
    psi_log = os.path.abspath(os.path.join(psi_dir, 'TransposonPSI.log'))
    psi_gff3 = os.path.abspath(os.path.join(psi_dir, args.out+'_Tpsi.gff3'))
    psi_results = os.path.abspath(os.path.join(psi_dir,'tpsi_library.fasta'))
    if os.path.exists(psi_log):
        os.remove(psi_log)
    with open(psi_log, 'w') as psi_log:
        try:
            subprocess.call([TRANSPOSONPSI, genome_file, 'nuc'], cwd=psi_dir, stdout=psi_log, stderr=psi_log)
        except:
            print('ERROR: There was an error running TransposonPSI, please check log file \
            located at {}'.format(psi_log))
            sys.exit()
    best_hits=[f for f in os.listdir(psi_dir) if os.path.isfile(os.path.join(psi_dir, f)) and 'bestPerLocus.gff3' in f]
    best_hits_gff = os.path.abspath(os.path.join(psi_dir, ''.join(best_hits)))
    with open(best_hits_gff, 'r') as best_hits, open(psi_gff3, 'w') as gff3:
        for line in best_hits:
            col = line.split('\t')
            target = re.search(r'Target=(\w+)', col[8]).group(1)
            gff3.write(col[0] + '\t' + col[1] + '\t' + target + '\t' + col[3] +
            '\t' + col[4] + '\t'+ col[5] + '\t' + col[6] + '\t' + col[7] + '\t' + col[8])
    subprocess.call([BEDTOOLS, 'getfasta', '-name', '-fi', genome_file, '-bed', psi_gff3, '-fo', 
    'tpsi_library.fasta'], cwd=psi_dir, stderr=subprocess.DEVNULL)
    os.remove(genome_file+'.fai')
    shutil.copy(psi_results, work_dir)

def concatenate_libraries(finder_lib, harvest_lib, tpsi_lib):
    concat_lib = os.path.abspath(os.path.join(class_dir, 'concatenated_library.fasta'))
    with open(concat_lib, 'w') as cat_lib:
        subprocess.call(['cat', finder_lib, harvest_lib, tpsi_lib], cwd=class_dir, stdout=cat_lib)
    shutil.copy(concat_lib, work_dir)

def run_repeatclassifier():
    print("Running RepeatClassifier: This can take some time depending on the size of the library")
    class_log = os.path.abspath(os.path.join(class_dir, 'repeatclassifier.log'))
    if os.path.exists(class_log):
        os.remove(class_log)
    with open(class_log, 'w') as cl_log:
        try:
            subprocess.call([REPEATCLASSIFIER, '-consensi', concatenated_lib, '-engine', args.engine
            ], cwd=class_dir, stdout=cl_log, stderr=cl_log)
        except:
            print('ERROR: There was an error running RepeatClassifier, please check log file \
            located at {}'.format(class_log))
            sys.exit()
    # don't know why, but output is deposited into parent directory, no need to shutil.copy

def prepare_usearch():
    usearch_prep = os.path.abspath(os.path.join(usearch_dir, 'usearch_prep.fasta'))
    # if os.path.exists(usearch_prep):
        # os.remove(usearch_prep)
    if args.repbase_lib:
        repbase_lib = os.path.abspath(os.path.join(rundir, args.repbase_lib))
        with open(usearch_prep, 'w') as u_prep:
            subprocess.call(['cat', repeatmodeler_lib, repbase_lib, classified_lib], cwd=usearch_dir, stdout=u_prep)
    else:
        with open(usearch_prep, 'w') as u_prep:
            subprocess.call(['cat', repeatmodeler_lib, classified_lib], cwd=usearch_dir, stdout=u_prep)
    shutil.copy(usearch_prep, work_dir)

def run_usearch():
    print("Running Usearch: This won't take too long.")
    usearch_run = os.path.abspath(os.path.join(usearch_dir, 'nonredundant_library.fasta'))
    usearch_log = os.path.abspath(os.path.join(usearch_dir, 'usearch.log'))
    if os.path.exists(usearch_log):
        os.remove(usearch_log)
    with open(usearch_log, 'a') as u_log:
        try:
            subprocess.call([USEARCH, '-cluster_fast', usearch_lib, '-id', str(args.identity), '-threads', 
            str(args.cpus), '-centroids', usearch_run, '-sort', 'length'
            ], cwd=usearch_dir, stdout=subprocess.DEVNULL, stderr=u_log)
        except:
            print('ERROR: There was an error running Usearch, please check log file \
            located at {}'.format(usearch_log))
            sys.exit()
    shutil.copy(usearch_run, work_dir)

def fix_calls():
    '''Fix up some Unknown calls, based on RepBase calls. This is not meant to be extensive.'''
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

def run_repeatmasker():
    print("Running RepeatMasker: In the home stretch now")
    masker_log = os.path.abspath(os.path.join(masker_dir, 'repeatmasker.log'))
    masked_g = os.path.abspath(os.path.join(masker_dir, args.out+'_assembly.fasta.masked'))
    repeat_table = os.path.abspath(os.path.join(masker_dir, args.out+'_assembly.fasta.tbl'))
    repeat_output = os.path.abspath(os.path.join(masker_dir, args.out+'_assembly.fasta.out'))
    if os.path.exists(masker_log):
        os.remove(masker_log)
    with open(masker_log, 'w') as mask_log:
        try:
            subprocess.call([REPEATMASKER, '-pa', str(args.cpus), '-gff', '-xm', '-lcambig', '-cutoff', '300', '-xsmall', 
            '-gccalc', '-dir', masker_dir, '-lib', final_lib, genome_file], cwd=masker_dir, stdout=mask_log, stderr=mask_log)
        except:
            print('ERROR: There was an error running RepeatMasker, please check log file \
            located at {}'.format(masker_log))
            sys.exit()
    shutil.copy(masked_g, os.path.join(result_dir, args.out+'_masked.fasta'))
    shutil.copy(repeat_table, os.path.join(result_dir, args.out+'_repeatmasker.tbl'))
    shutil.copy(repeat_output, os.path.join(result_dir, args.out+'_repeatmasker.output'))

if __name__ == "__main__":
    if not args.repeatmodeler_lib: #no pre-computed repeat library is given, so must run RepeatModeler
        repeatmodeler_lib = os.path.abspath(os.path.join(work_dir, 'modeler_library.fasta'))
        if not os.path.exists(repeatmodeler_lib):
            run_repeatmodeler()
        else: # assume repeatmodeler ran successfully as final file is correctly located
            pass
    else:
        repeatmodeler_lib = os.path.abspath(os.path.join(rundir, args.repeatmodeler_lib))

    ltr_finder_lib = os.path.abspath(os.path.join(work_dir, 'ltrfinder_library.fasta'))
    if not os.path.exists(ltr_finder_lib):
        run_ltr_finder()

    ltr_harvest_lib = os.path.abspath(os.path.join(work_dir, 'harvest_library.fasta'))
    if not os.path.exists(ltr_harvest_lib):
        run_ltr_harvest_and_digest()

    tpsi_lib = os.path.abspath(os.path.join(work_dir, 'tpsi_library.fasta'))
    if not os.path.exists(tpsi_lib):
        run_transposonPSI()

    concatenated_lib = os.path.abspath(os.path.join(work_dir, 'concatenated_library.fasta'))
    if not os.path.exists(concatenated_lib):
        concatenate_libraries(ltr_finder_lib, ltr_harvest_lib, tpsi_lib)

    classified_lib = os.path.abspath(os.path.join(work_dir, 'concatenated_library.fasta.classified'))
    if not os.path.exists(classified_lib):
        run_repeatclassifier()

    usearch_lib = os.path.abspath(os.path.join(work_dir, 'usearch_prep.fasta'))
    if not os.path.exists(usearch_lib):
        prepare_usearch()

    usearch_results = os.path.abspath(os.path.join(work_dir, 'nonredundant_library.fasta'))
    if not os.path.exists(usearch_results):
        run_usearch()

    final_lib = os.path.abspath(os.path.join(work_dir, 'final_library.fasta'))
    if not os.path.exists(final_lib):
        fix_calls()

    masked_genome = os.path.abspath(os.path.join(result_dir, args.out+'_masked.fasta'))
    repeatmasker_tbl = os.path.abspath(os.path.join(result_dir, args.out+'_repeatmasker.tbl'))
    repeatmasker_out = os.path.abspath(os.path.join(result_dir, args.out+'_repeatmasker.output'))
    if not os.path.exists(masked_genome) and not os.path.exists(repeatmasker_tbl) and not os.path.exists(repeatmasker_out):
        run_repeatmasker()

    print("Job done! Please find your masked genome, repeatmasker table, and repeatmasker output", \
    "in the directory {}/".format(os.path.abspath(args.out)))