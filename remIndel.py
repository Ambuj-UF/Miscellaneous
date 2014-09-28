################################################################################################################
# Tool for removing indels and stop codons from the alignment file                                             #
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


import argparse
import textwrap
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

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

args = parser.parse_args()

handle = open(args.i, 'rU')
records = list(SeqIO.parse(handle, args.itype))

def split(str, num):
    return [ str[start:start+num] for start in range(0, len(str), num) ]


def main():
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

