#!env/bin/python

from optparse import OptionParser
from cartselection.main import main, helper
import os


# Print out info
def print_info(options):
    print '\nInput file (HGNC IDs):      ' + options.hgnc
    print 'APPRIS RefSeq file:         ' + options.appr
    print 'RefSeq transcript db file:  ' + options.refsdb
    print 'RefSeqScan output file:     ' + options.refss
    print 'Gene ID dictionary file:    ' + options.genes
    print '\nCART series: ' + options.series.upper() + '\n'

# Version
ver = 'v1.7.0'

# Script dir
scriptdir = os.path.dirname(os.path.realpath(__file__))

# Command line argument parsing
descr = 'CART selection script'
parser = OptionParser(usage='python cart_selection <options>', version=ver, description=descr)
parser.add_option("--hgnc", dest='hgnc', action='store', help="Input file containing HGNC IDs")
parser.add_option("--appr", dest='appr', action='store', help="APPRIS file")
parser.add_option("--refsdb", dest='refsdb', action='store', help="RefSeq transcript database file")
parser.add_option("--refss", dest='refss', action='store', help="refseq_scan output file")
parser.add_option("--genes", dest='genes', action='store', help="Gene ID dictionary file")
parser.add_option("--build", dest='build', default=None, action='store', help="Genome build")
parser.add_option("--series", dest='series', default=None, action='store', help="CART series")
parser.add_option("--out", dest='out', action='store', help="Output file name")
(options, args) = parser.parse_args()

if options.series is None:
    print '\nUse command line option --series to specify the CART series.\n'
    quit()

print '\n'+'='*100
print '\nCART selection script ' + ver

options = helper.convert_file_names_to_absolute_path(options)

print_info(options)

main(options)

print '- Done'
print '\nOutput files written: '
print ' - ' + options.out + '_selected.txt'
print ' - ' + options.out + '_missing.txt'
print ' - ' + options.out + '_log.txt'
print '\n'+'='*100+'\n'


