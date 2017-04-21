#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Li Yao
# @File: prepare-data-for-giremi.py
# @License: MIT
# @Bitbutcket: https://bitbucket.org/li_yao/
# @Github: https://github.com/liyao001

import pysam
import getopt
import sys


def is_snp(table, chrom, pos):
    flag = 0
    try:
        for hit in table.fetch(reference=chrom, start=int(pos) - 1, end=int(pos), parser=pysam.asGTF()):
            flag = 1
            break
    except:
        pass
    return str(flag)


def get_gene_annotation(table, chrom, pos):
    gene_symbol = 'Inte'
    gene_strand = '#'
    try:
        for hit in table.fetch(reference=chrom, start=int(pos) - 1, end=int(pos), parser=pysam.asGTF()):
            gene_symbol = hit.gene_id
            if hit.strand != '.':
                gene_strand = hit.strand
    except:
        pass
    return gene_symbol, gene_strand


def helper(input_file, snp_file, gene_annotation_file):
    snpFile = pysam.Tabixfile(snp_file)
    refSeqFile = pysam.Tabixfile(gene_annotation_file)
    rawInput = open(input_file, mode='r')
    newFile = open(input_file + ".giremi", mode='w')
    new_lines = []
    first_line = 1
    i = 0
    for line in rawInput.readlines():
        i+=1
        if first_line == 1:
            first_line = 0
            continue
        cols = line.replace('\n', '').replace('\r', '').split('\t')
        coordinate = cols[0]
        position = cols[1]
        gene_symbol, gene_strand = get_gene_annotation(refSeqFile, coordinate, position)
        snp_mark = is_snp(snpFile, coordinate, position)
        record = (coordinate, str(int(position) - 1), position, gene_symbol, snp_mark, gene_strand + '\n')
        new_lines.append('\t'.join(record))
    newFile.writelines(new_lines)
    newFile.close()
    rawInput.close()

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:s:g:", ["input_file=", "snp_file=", "gene_file="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit()
    input_file = None; snp_file = None; gene_annot = None;
    if len(opts)==0:
        help()
        sys.exit()
    for o, a in opts:
        if o in ("-i","--input_file"):
            input_file = a
        elif o in ("-s", "--snp_file"):
            snp_file = a
        elif o in ("-g", "--gene_file"):
            gene_annot = a
    print 'Please wait......'
    helper(input_file, snp_file, gene_annot)
    print 'Done!'

if __name__ == '__main__':
    main()
