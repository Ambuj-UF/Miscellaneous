################################################################################################################
# Tool for removing indels and stop codons from the alignment file                                             #
# It also performes species centric alignment correction                                                       #
#                                                                                                              #
# Copyright (C) {2014}  {Ambuj Kumar, Kimball-Braun lab group, Biology Department, University of Florida}      #
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


import sys
import argparse
import textwrap

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna
    from Bio.Align import MultipleSequenceAlignment as msa
except:
    sys.exit("Biopython not found on the system")

parser = argparse.ArgumentParser(prog='Remove-Indel',
                                 version= '1.0',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
    ----------------------------------------------------------------------------------------------------------
    \t\t\t Designed at Kimbal-Braun Lab Group, University of Florida
    
    ----------------------------------------------------------------------------------------------------------
    
    '''))

parser.add_argument('-i', type=str, required = True, help='Enter input alignment filename')
parser.add_argument('-o', type=str, required = True, help='Enter output alignment filename')
parser.add_argument('-itype', type=str, required = True, choices=['fasta', 'nexus', 'phylip', 'phylip-interleived', 'phylip-relaxed'], help='Enter input alignment file format')
parser.add_argument('-otype', type=str, required = True, choices=['fasta', 'nexus', 'phylip', 'phylip-interleived', 'phylip-relaxed'], help='Enter output alignment file format')
parser.add_argument('-c', action='store_true', default=False, help='Use if you want species centric alignment correction')
parser.add_argument('-name', type=str, default = None, required = True, help='Enter species name')


args = parser.parse_args()

if args.c == True and args.name == None:
    print("-name required in -c mode\n")


handle = open(args.i, 'rU')
records = list(SeqIO.parse(handle, args.itype))

def split(str, num):
    return [ str[start:start+num] for start in range(0, len(str), num) ]

def group(L):
    first = last = L[0]
    for n in L[1:]:
        if n - 1 == last:
            last = n
        else:
            yield first, last
            first = last = n
    yield first, last


def checkSpecies(records):
    posList = list()
    for rec in records():
        if rec.id == args.name:
            dataSplit = split(rec.seq)
            for i, data in enumerate(dataSplit):
                if '-' in data:
                    posList.append(i)

    groupData = [x for x in group(posList)]
    msaObject = msa(records)
    for val in groupData:
        if dataSplit[val[0]].count('-') != 3 and dataSplit[val[1]].count('-') != 3:
            msaObject = msaObject[0:(val[0]-1)*3] + msaObject[(val[1]*3)+1:len(msaObject)]

    for i, rec in enumerate(records):
        records[i].seq = Seq(msaObject[rec.id], generic_dna)

    return records



def main():
    
    if args.c == True:
        records = checkSpecies(records)
    
    for i, rec in enumerate(records):
        newSeq = Seq("", generic_dna)
        seqData = split(rec.seq, 3)
        for j, data in enumerate(seqData):
            if '-' in data and data.count('-') != 3 or 'TGA' in data or 'TAG' in data or 'TAA' in data:
                seqData[j] = Seq("---", generic_dna)
        for newData in seqData:
            newSeq = newSeq + newData

        records[i].seq = newSeq

    with open(args.o, 'w') as fp:
        SeqIO.write(records, fp, args.otype)



if __name__ == "__main__":
    main()

