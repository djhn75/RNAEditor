'''
Created on May 13, 2013

@author: david
'''

class Editing(object):
    '''
    classdocs
    '''


    def __init__(selfparams):
        '''
        Constructor
        '''


parser = argparse.ArgumentParser(description='map FastQ Files to the given genome and realigns the reads for SNP-calling.')
parser.add_argument('-i', '--input', metavar='Fastq-File', type=argparse.FileType('r'), help='Input fastq file that should be mapped to the genome', required=True)
parser.add_argument('-o', '--output', metavar='output-prefix', type=str,help='prefix that is written in front of the output files', required=True)
parser.add_argument("-r", "--RefGenome", metavar='Fasta-File', help="File that contains the reference sequences", type=argparse.FileType('r'), default='/media/databases/human/human_g1k_v37.fasta')
parser.add_argument('-s', '--dbsnp', help=' SNP database (dbSNP) in VCF format (downloaded from the GATK homepage)', type=argparse.FileType('r'), default='/media/databases/human/dbsnp_135.b37.vcf')
parser.add_argument('-d', '--sourceDir', help='- Directory to all the tools [default: /bin]', default='/media/', type=readable_dir)
parser.add_argument('-n', '--maxDiff', help=' maximum Number of mismatches in the reads (int) or error rate in percentage (float)[0.04]', type=float, default=0.04)
parser.add_argument('-t', '--threads', help='number of threads', type=int, default=8)
parser.add_argument('--seedDiff', help='maximum Number of mismatches in the seed sequence (int)[2]', type=int, default=2)
parser.add_argument('-p', '--paired', help="Use this paramater if you have paired end reads [false]", action='store_true', default=False)
parser.add_argument('--keepTemp', help='keep the intermediate Files [False]', action='store_true', default=False)
parser.add_argument('--overwrite', help='overwrite existing Files [False]', action='store_true', default=False)
