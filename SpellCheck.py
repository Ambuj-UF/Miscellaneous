################################################################################################################                                #
# This program checks taxon name spelling mistakes in your alignment files.                                    #
# Put this program in the directory that contains all your alignment files.                                    #
# Copyright (C) {2014}  {Ambuj Kumar, Kimball-Brain lab group, Biology Department, University of Florida}      #
#                                                                                                              #
# This program is free software: you can redistribute it and/or modify                                         #
# it under the terms of the GNU General Public License as published by                                         #
# the Free Software Foundation, either version 3 of the License, or                                            #
# (at your option) any later version.                                                                          #
#                                                                                                              #
# This program is distributed in the hope that it will be useful,                                              #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                #
# GNU General Public License for more details.                                                                 #
#                                                                                                              #
# This program comes with ABSOLUTELY NO WARRANTY;                                                              #
# This is free software, and you are welcome to redistribute it                                                #
# under certain conditions;                                                                                    #
#                                                                                                              #
################################################################################################################




import glob
import sys
from Bio import SeqIO
from Bio.AlignIO import MultipleSeqAlignment


print "Select the input file format \n 1. fas \n 2. nex \n 3. phy \n"

file_format = input('\n')
    
extList = ["*.fas", "*.nex", "*.phy", "*.phy", "*.phy"]
typeList = ["fasta", "nexus", "phylip", "phylip-sequential", "phylip-relaxed"]

try:
    fileList = glob.glob(extList[file_format - 1])
except NameError:
    sys.exit("Files with extension %s not found \n" %extList[file_format - 1])

idDict = {}
    
for filename in fileList:
    handle = open(filename, "rU")
    idList = []
    for record in SeqIO.parse(handle, typeList[file_format - 1]):
        idList.append(record.id)
    gene = filename.split(".")[0]
    idDict[gene] = idList
    handle.close()

idList = []
keyList = []
    
for key, val in idDict.items():
    try:
        idList.append([x.split('|')[0] for x in val])
    except KeyError:
        idList.append(val)
        
    keyList.append(key)

counter = 1
    
for value in idList:
    n = counter
    while n < len(idList):
        for i, j in zip(value, idList[n]):
            if 1.0 > float([x == y for (x, y) in zip(i, j)].count(True))/len(j) > 0.8:
                print "Found " + i + " in gene " + keyList[counter-1] + " but gene " + keyList[n] + " (has " + j+")\n"
            
        n = n + 1
    counter = counter + 1

print "I can correct these spelling mistakes...If you want me to..;) \n Do you?? \n 1. Yes \n 2. No \n"

choice = input('\n')
while choice == 1:
    print "Enter the correct taxon name \n"
    usrInp = raw_input('\n')

    for filename in fileList:
        handle = open(filename, 'rU')
        record = list(SeqIO.parse(handle, 'nexus'))
        msa = MultipleSeqAlignment(record)
        for i, val in enumerate(msa):
            if 1.0 > float([x == y for (x, y) in zip(usrInp, msa[i].id)].count(True))/len(msa[i].id) > 0.8:
                msa[i].id = usrInp
            else:
                print "No Spelling mistakes found in file %s \n" % filename

        fp = open(filename, 'w')
        SeqIO.write(msa, fp, typeList[file_format-1])
        handle.close()
        fp.close()
    print "Want me to do some more editing?? \n Well! I can do this whole day \n What about you?? \n 1. Yes \n 2. No \n"
    choice = input('\n')
