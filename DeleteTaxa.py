################################################################################################################
# This script deletes user suplied taxons from all the alignment files present in the directory.               #
# Run this script in a folder that contains all your alignment files and enter taxaon name which               #
# you want to delete.                                                                                          #
#                                                                                                              #
# To run this script you require biopython installed on your system                                            #
#                                                                                                              #
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
from Bio import SeqIO
from Bio.AlignIO import MultipleSeqAlignment


files = glob.glob("*.nex")
print files
print "Enter the taxon name \n"
usrID = raw_input('\n')
for filename in files:
    handle = open(filename, 'rU')
    record = list(SeqIO.parse(handle, "nexus"))
    handle.close()
    fp = open(filename, 'w')

    msa = MultipleSeqAlignment(record)
    for i, val in enumerate(msa):
        if val.id.split('|')[0] == usrID:
            align1 = msa[ :i, :]
            align2 = msa[i+1:len(msa), :]
            align1.extend(align2)
            msa = align1

    SeqIO.write(msa, fp, "nexus")
    fp.close()
    print "%s taxon deleted from all the files \n" % usrID
