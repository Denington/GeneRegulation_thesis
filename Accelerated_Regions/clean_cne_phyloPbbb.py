#a quick script for taking phyloP wigs spit out into phylop_multiAcc_wigs_fullcons and similar, concatenating

import glob
import os


import sys


folder = sys.argv[1]
outname = sys.argv[2]
outfolder = sys.argv[3]

os.chdir(folder)
#no_pruning
#all_but


for i in range(14):
	excluding_accs = sorted(glob.glob('%s.*no_pruning*'%i))
	including_accs = sorted(glob.glob('%s.*all_but*'%i))

	with open('%s/%s.accels.all_but.%s.wig'%(outfolder,i,outname),'w') as out:
		for block in excluding_accs:
			o = open(block,'r')
			header = o.readline() #remove wig header
			goodheader = header[:header.find(':')] + header[header.find(' start'):]
			out.write(goodheader)
			for line in o:
				if ':' in line:
					line = line[:line.find(':')] + line[line.find(' start'):]
				out.write(line)
			#out.write(o.read())
		out.close()

	with open('%s/%s.accels.including.accs.%s.wig'%(outfolder,i,outname),'w') as out:
		for block in including_accs:
			o = open(block,'r')
			header = o.readline() #remove wig header
			goodheader = header[:header.find(':')] + header[header.find(' start'):]
			out.write(goodheader)
			for line in o:
				if ':' in line:
					line = line[:line.find(':')] + line[line.find(' start'):]
				out.write(line)
			#out.write(o.read())
		out.close()





#			for line in o:
#				if line[0] == 'f':
#					continue
#				out.write(line)
#			o.close()

#	with open('%s/%s.accels.no_pruning.%s'%(outfolder,i,outname),'w') as out:
#		for block in including.accs:
#			o = open(block,'r')
#			o.readline() #remove wig header
#			for line in o:
#				if line[0] == 'f':
#					continue
#				out.write(line)
#			o.close()
