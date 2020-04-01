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

    --txDb            file containing annotation database for a current genetic annotation
                      can be generated with the GenomicFeatures Bioconductor package
                      ideally use refGene table for generation (Default: packaged version)

    --

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
