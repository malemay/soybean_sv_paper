# MIT License
# 
# Copyright (c) 2019 Benedict Paten
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
    # The above copyright notice and this permission notice shall be included in all
    # copies or substantial portions of the Software.
# 
    # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    # SOFTWARE.
# 
from pyfaidx import Fasta
import fileinput

# Load reference genome
fa = Fasta('../../refgenome/Gmax_508_v4.0_mit_chlp.fasta')

# stream VCF
for line in fileinput.input():
    # if header, print and got to next line
    if line[0] == '#':
        line = line.rstrip()
        print line
        continue
    # else parse variant record
    line = line.rstrip().split('\t')
    chrom = line[0]
    pos = int(line[1])
    ref = line[3]
    alt = line[4]
    if ref[0] == alt[0]:
        # if variant and padding present, print
        print '\t'.join(line)
    else:
        # else missing padding
        # find padding base
        nuc = fa[str(chrom)][pos - 2]
        line[3] = nuc.seq + line[3]
        line[4] = nuc.seq + line[4]
        # shift position
        line[1] = str(pos - 1)
        # print padded record
        print '\t'.join(line)

