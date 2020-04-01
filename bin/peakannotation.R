#!/usr/bin/env Rscript
suppressMessages(library(ChIPseeker))
library(argparser)
suppressMessages(library(GenomicFeatures))

parser <- arg_parser('generates an annotation file for a given ChIP peak file in BED format')
parser <- add_argument(parser, 'peakfile', help='file containing the peak regions in BED format')
parser <- add_argument(parser, 'txdb', help='TxDb annotation database generated from a current annotation')
parser <- add_argument(parser, 'outputFile', help='directory to which the generated annotation should be written')

argv <- parse_args(parser)
#making new txdb for annotation from current ucsc knownGene table
#txdb <-makeTxDbFromUCSC(genome = 'mm9', tablename = 'refGene')
#saving txdb to file
#saveDb(txdb, file = '/Users/daniel.malzl/PycharmProjects/ori_class/mm9_txdb.sqlite')

txdb <- loadDb(argv$txdb)
#tssRegion parameter defines the region for which we still report a peak as Promoter
#read peaks
peaks <- readPeakFile(argv$peakfile)

#annotate peaks
genes <- annotatePeak(peaks, tssRegion = c(0, 0), TxDb = txdb,
                      genomicAnnotationPriority = c("Promoter", "Exon", "Intron", "5UTR", "3UTR", "Intergenic", "Downstream"),
                      ignoreDownstream=FALSE)

#convert to data.frame
geneframe <- as.data.frame(genes)

#writing to files
write.table(geneframe, file = args$outputFile, row.names = FALSE, sep ='\t', quote = FALSE)
