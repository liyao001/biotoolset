import os, getopt, sys

"""
PBS Script Builder
-f RNA-seq files' storing folder
-t template file
-p prefix
-e file extension
-c cut length
&copy; Li Yao 2016.
yaoli95@outlook.com
"""

tree = {}

def GetFileList(dir, fileList):
    newDir = dir
    if os.path.isfile(dir):
        fileList.append(dir.decode('gbk'))
    elif os.path.isdir(dir):  
        for s in os.listdir(dir):
            newDir=os.path.join(dir,s)
            GetFileList(newDir, fileList)  
    return fileList

try:
    opts, args = getopt.getopt(sys.argv[1:], "f:t:p:e:c:", ["folder=", "template=", "prefix=", "extension=", "cut="])
except getopt.GetoptError, err:
    print str(err) 
    sys.exit()

flag = 0
fd = None; tp = None; pre = None; ext = None; cut = 0;

if len(opts)==0:
    sys.exit()
for o, a in opts:
    if o in ("-f","--folder"):
        fd = a
    elif o in ("-t", "--template"):
        tp = a
    elif o in ("-p", "--prefix"):
        pre = a
    elif o in ("-e", "--extension"):
        ext = a
    elif o in ("-c", "--cut"):
        cut = int(a)

if fd == None or tp == None or pre == None or ext == None or cut == 0:
    sys.exit("Lack arguments!")

allFiles = GetFileList(fd, [])
#print allFiles
tf = open(tp, 'r')
pbsScript = tf.read()
used = []
experiment = {}
for file in allFiles:
    #if file.find('.fastq') != -1:
    fname, fextension = os.path.splitext(file)
    if fextension == '.'+ext:
        err = os.path.basename(file)
        exp = err[:-1*cut]
        if experiment.has_key(exp):
            experiment[exp].append(file)
        else:
            experiment[exp] = [file]
ekeys = experiment.keys()

for exp in ekeys:
    tmp = open(pre+"-"+exp+".pbs", 'w')
    template = pbsScript.replace('{PREFIX}', pre).replace('{ERRCODE}', exp)
    for k, v in enumerate(experiment[exp]):        
        template = template.replace('{FILE'+str(k+1)+'}', v)
    tmp.write(template)
    tmp.close()
    os.system("qsub "+pre+"-"+exp+".pbs")
