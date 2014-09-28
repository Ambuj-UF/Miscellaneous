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



from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

handle = open('ASPM.phy', 'rU')
records = list(SeqIO.parse(handle, 'phylip-relaxed'))

def split(str, num):
    return [ str[start:start+num] for start in range(0, len(str), num) ]

for i, rec in enumerate(records):
    newSeq = Seq("", generic_dna)
    seqData = split(rec.seq, 3)
    for j, data in enumerate(seqData):
        if '-' in data and data.count('-') != 3:
            seqData[j] = Seq("---", generic_dna)
        elif 'TGA' in data or 'TAG' in data or 'TAA' in data:
            seqData[j] = Seq("---", generic_dna)
    for newData in seqData:
        newSeq = newSeq + newData

    records[i].seq = newSeq

with open('output.phy', 'w') as fp:
    SeqIO.write(records, fp, 'phylip-relaxed')

