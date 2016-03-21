#! /usr/bin/python

import os
from pandas import DataFrame
import pandas as pd
from math import log
mapFile = ''
tcgaMap = None

def getNormalizedFileList(baseDir):
    result = []
    os.chdir(baseDir)
    for fileOrFolder in os.listdir(os.curdir):
        if os.path.isdir(fileOrFolder):
            print 'Please check your folder, folder exists!'
        else:
            if fileOrFolder.find("genes.normalized_results") != -1:
                result.append(os.getcwd()+os.sep+fileOrFolder)
    return result

def getGeneExpression(normalizedFiles, geneSymbol, output):
    result = []
    
    for file in normalizedFiles:
        table = pd.read_table(file)
        condition = table.gene_id.str.contains(geneSymbol)
        where = table[condition]
        if len(where) == 1:
            tmp = [list(where['gene_id'])[0], str(list(where['normalized_count'])[0]), getBarcode(mapFile, os.path.basename(file)), str(log(float(list(where['normalized_count'])[0]), 2))]
            result.append(tmp)
        else:
            print 'err'
        #break
    sep = '\t'
    toWrite = open(output, 'w')
    
    for tmp in result:
        line = '\t'.join(tmp)
        toWrite.write(line+'\n')
    toWrite.close()

def getBarcode(mapFile, currentFile):
    global tcgaMap
    if tcgaMap is None:
        tcgaMap = pd.read_table(mapFile)
    condition = tcgaMap.filename.str.contains(currentFile)
    result = tcgaMap[condition]
    if len(result) > 0:
        return str(list(result['barcode(s)'])[0])
    else:
        return 'err'
if __name__ == '__main__':
    mapFile = raw_input("Where is the map file: ")
    data = raw_input("Where is the data folder: ")
    fl = getNormalizedFileList(data)
    gene = raw_input("Gene symbol: ")
    store = raw_input("Path for storing the result: ")
    getGeneExpression(fl, 'YY1\|', store)
