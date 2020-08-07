#!/usr/bin/env nextflow

/*
* MIT License
*
* Copyright (c) 2020 Daniel Malzl
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/
def helpMessage() {
  log.info"""
  ================================================================
   geneIS-nf
  ================================================================
  DESCRIPTION

  Reproducible annotation of initiation sites with genetic features

  Usage:
  nextflow run dmalzl/geneIS-nf

  Options:
    --masterTable     tab-separated file containing the genomic coordinates and
                      names of the mapped initiation sites with header row and
                      first four columns being chr, start, end, name (same as BED)
                      including quantification results

    --txDb            file containing annotation database for a current genetic annotation
                      can be generated with the GenomicFeatures Bioconductor package
                      ideally use refGene table for generation

    --email           email address to use to fetch genesymbols from EntrezIds

    --xCol            column of masterTable holding quantification for x-axis of plots
    --yCol            column of masterTable holding quantification for y-axis of plots

    --axMin           minimum value of axes of final plots (Default: 0)
    --axMax           maximum value of axes of final plots (Default: 8)

    --foldChange      adds diagonal lines in distance of foldChange to plot

    --filePrefix      prefix to use for output files (Default: ISgene)
    --outputDir       name of the directory to save results to (Default: results)

    Authors:
    Daniel Malzl (daniel.malzl@imp.ac.at)
  """.stripIndent()
}

params.help = false

if (params.help) {
  helpMessage()
  exit 0
}

if (params.masterTable) {
  if (!file(params.masterTable).exists()) {
    exit 1, "--masterTable was specified but does not exist"
  }
} else {
  exit 1, "--masterTable is requried but was not set!"
}

if (params.txDb) {
  if (!file(params.txDb).exists()) {
    exit 1, "--txDb was specified but does not exist"
  }
} else {
  exit 1, "--txDb is required but was not set!"
}

if (!params.email) {
  exit 1, "--email needs to be set"
}

log.info ""
log.info " parameters"
log.info " ======================"
log.info " masterTable              : ${params.masterTable}"
log.info " txDb                     : ${params.txDb}"
log.info " email                    : ${params.email}"
log.info " xCol                     : ${params.xCol}"
log.info " yCol                     : ${params.yCol}"
log.info " axMin                    : ${params.axMin}"
log.info " axMax                    : ${params.axMax}"
log.info " filePrefix               : ${params.filePrefix}"
log.info " outputDir                : ${params.outputDir}"
log.info " ======================"
log.info ""

inputChannel = Channel
                  .fromList([[params.filePrefix,
                              file(params.masterTable),
                              file(params.txDb)]])

process computeGeneAnnotation {

  tag { filePrefix }

  input:
  set val(filePrefix), file(masterTable), file(txDb) from inputChannel

  output:
  set val(filePrefix), file("${filePrefix}.chipseeker.tsv"), file(masterTable) into resultsGeneAnnotation

  shell:
  '''
  cut -f 1,2,3,4 !{masterTable} | tail -n +2 > master.tmp.bed
  annotateIniSites.R master.tmp.bed !{txDb} !{filePrefix}.chipseeker.tsv
  '''
}

process mapEntrezIds {

  tag { filePrefix }

  publishDir  path: "${params.outputDir}",
              mode: "copy",
              overwrite: "true",
              pattern: "*.chipseeker.mapped.tsv"
              
  publishDir  path: "${params.outputDir}",
              mode: "copy",
              overwrite: "true",
              pattern: "entrezEntries*"

  input:
  set val(filePrefix), file(annotation), file(masterTable) from resultsGeneAnnotation

  output:
  set val(filePrefix), file("${filePrefix}.chipseeker.mapped.tsv"), file(masterTable) into resultsMapEntrez
  set file("entrezEntries.txt"), file("entrezEntries.txt.ids") into entrezResults

  shell:
  '''
  mapentrez.py -a !{annotation} \
               --email !{params.email} \
               -m entrezEntries.txt \
               -o !{filePrefix}.chipseeker.mapped.tsv
  '''
}

process extendMasterTable {

  tag { filePrefix }

  publishDir  path: "${params.outputDir}",
              mode: "copy",
              overwrite: "true",
              pattern: "${filePrefix}.master.tsv"

  input:
  set val(filePrefix), file(mappedAnnotation), file(masterTable) from resultsMapEntrez

  output:
  set val(filePrefix), file("${filePrefix}.master.tsv"), file("${filePrefix}.master.tsv.plotTable") into resultsExtendMaster

  shell:
  '''
  annotateMasterTable.py -mt !{masterTable} -a !{mappedAnnotation} -o !{filePrefix}.master.tsv
  '''
}

process plotting {

  tag { filePrefix }

  publishDir  path: "${params.outputDir}",
              mode: "copy",
              overwrite: "true",
              pattern: "*.pdf"

  input:
  set val(filePrefix), file(masterTable), file(plotTable) from resultsExtendMaster

  output:
  set file("${filePrefix}.Promoter.pdf"), file("${filePrefix}.Exon.pdf"), file("${filePrefix}.Intron.pdf"), file("${filePrefix}.Intergenic.pdf") into resultsPlotting

  shell:
  '''
  for geneFeature in Promoter Exon Intron Intergenic;
  do
    datashader.py -mt !{plotTable} \
                  --xcol !{params.xCol} --ycol !{params.yCol} \
                  --xmin !{params.axMin} --xmax !{params.axMax} \
                  --ymin !{params.axMin} --ymax !{params.axMax} \
                  --subsetColumn ${geneFeature} \
                  --xlabel "!{params.xCol} log2(RPM)" \
                  --ylabel "!{params.yCol} log2(RPM)" \
                  --density F T \
                  --colormaps Grey,Grey coolwarm \
                  --foldchange !{params.foldChange} \
                  --plotmethod mesh \
                  --figwidth 6 --figheight 6 \
                  -o !{filePrefix}.${geneFeature}.pdf \
                  --xbins 200 --ybins 200 \
                  --labels rest $geneFeature \
                  --plotCounts
    done
  '''
}
