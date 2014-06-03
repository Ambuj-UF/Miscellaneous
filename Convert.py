'''Use this program to convert your alignment file to any other alignement file
    format. Run this program in the directory where your alignment files are 
    present '''


from Bio import SeqIO
from Bio.Alphabet import IUPAC, Gapped
import glob

print "Select the input file format \n 1. Fasta[*.fas] \n 2. Nexus[*.nex] \n 3. Phylip-Interlaced[*.phy] \n 4) Phylip-sequential[*.phy] \n 5) Relaxed Phylip[*.phy] \n"

usrInp = input('\n')

print "Select the output file format \n 1. Fasta[*.fas] \n 2. Nexus[*.nex] \n 3. Phylip-Interlaced[*.phy] \n 4) Phylip-sequential[*.phy] \n 5) Relaxed Phylip[*.phy] \n"

usrInpOut = input('\n')

file_format = ["*.fas", "*.nex", "*.phy", "*.phy", "*.phy"]
extType = ["fasta", "nexus", "phylip", "phylip-sequential", "phylip-relaxed"]

filenames = glob.glob(file_format[usrInp-1])

if usrInp == 1 and usrInpOut == 2:
    for files in filenames:
        alignment = AlignIO.read(open(files), "fasta", alphabet=Gapped(IUPAC.protein))
        g = open(files.split(".")[0] + file_format[usrInpOut - 1].lstrip('*'), 'w')
        g.write (alignment.format("nexus"))
        g.close()

else:
    try:
        for files in filenames:
            handle = open(files, 'rU')
            record = list(SeqIO.parse(handle, extType[usrInp - 1]))
            fp = open(files.split(".")[0] + file_format[usrInpOut - 1].lstrip('*'), 'w')
            SeqIO.write(record, fp, extType[usrInpOut - 1])
            fp.close()
            handle.close()
    except:
        print "Bad Bad alignments \n Program terminated"
