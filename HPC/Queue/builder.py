#this python script is designed for building pbs files
import os

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

allFiles = GetFileList(r'YOUR PATH', [])
tf = open('YOUR TEMPLATE', 'r')
template = tf.read()
used = []
for file in allFiles:
    tissue = os.path.basename(os.path.dirname(file))
    err = os.path.basename(file)
    exp = err[-20:-11]
    if exp not in used:
        used.append(exp)
    else:
        continue
    tmp = open(err+".pbs", 'w')
    tmp.write(template.replace('{ERRCODE}', err).replace('{TISSUE}', tissue).replace('{EXP}', exp))
    tmp.close()
    
tf.close()
