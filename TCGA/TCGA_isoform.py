#! /usr/bin/python
# Copyright (C) 2016 Li Yao <yaoli95@outlook.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import os
from pandas import DataFrame
import pandas as pd
from math import log
mapFile = ''
tcgaMap = None
isoformId = []

def getNormalizedFileList(baseDir):
    result = []
    os.chdir(baseDir)
    for fileOrFolder in os.listdir(os.curdir):
        if os.path.isdir(fileOrFolder):
            print 'Please check your folder, folder exists!'
        else:
            if fileOrFolder.find("isoforms.normalized_results") != -1:
                result.append(os.getcwd()+os.sep+fileOrFolder)
    return result

def getIsoformExpression(normalizedFiles, isos, output):
    result = []
    
    for file in normalizedFiles:
        for iso in isos:
            table = pd.read_table(file)
            condition = table.isoform_id.str.contains(iso)
            where = table[condition]
            if len(where) == 1:
                if float(list(where['normalized_count'])[0]) == 0.0:
                    logValue = 'NaN'
                else:
                    tmp = [list(where['isoform_id'])[0], str(list(where['normalized_count'])[0]), getBarcode(mapFile, os.path.basename(file)), str(log(float(list(where['normalized_count'])[0]), 2))]
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
    while True:
        id = raw_input("Isoform ID(Enter q to quit): ")
        if id != 'q':
            isoformId.append(id)
        else:
            break
    store = raw_input("Path for storing the result: ")
    getIsoformExpression(fl, isoformId, store)
