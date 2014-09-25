################################################################################################################
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

import csv
from Bio import SeqIO

with open('trait.csv', 'rb') as csvfile:
  csvData = csv.reader(csvfile, delimiter=' ', quotechar='|')
  sp_name = [row for row in csvData]
    
handle = open('alignment.fasta', 'rU')
records = list(SeqIO.parse(handle, 'fasta'))

newRec = list()
for rec in records:
  if rec.id in sp_name:
    newRec.append(rec)
    
with open('output.fasta', 'w') as fp:
  SeqIO.write(newRec, fp, 'fasta')
