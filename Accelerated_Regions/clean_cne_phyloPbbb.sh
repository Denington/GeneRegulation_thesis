#!/bin/sh


infolder=/Users/jk2014/Documents/Projects/Accelerated_Regions/phylop_multiAcc_wigs_mammalonly
outname=woofboy
outfolder=/Users/jk2014/Documents/Projects/Accelerated_Regions/
python /Users/jk2014/Documents/Projects/Accelerated_Regions/scripts/clean_cne_phyloPbbb.py $infolder $outname $outfolder

wigToBigWig $outfolder/0.accels.all_but.$outname.wig -clip /Users/jk2014/Downloads/Moving_stuff_to_laptop/Presentations_documents/hg19.chrom.sizes  $outfolder/0.accels.$outname.bw
wigToBigWig $outfolder/0.accels.including.$outname.wig -clip /Users/jk2014/Downloads/Moving_stuff_to_laptop/Presentations_documents/hg19.chrom.sizes  $outfolder/0.accels.including.$outname.bw

wigToBigWig $outfolder/1.accels.all_but.$outname.wig -clip /Users/jk2014/Downloads/Moving_stuff_to_laptop/Presentations_documents/hg19.chrom.sizes  $outfolder/1.accels.excluding.$outname.bw
wigToBigWig $outfolder/1.accels.including.$outname.wig -clip /Users/jk2014/Downloads/Moving_stuff_to_laptop/Presentations_documents/hg19.chrom.sizes  $outfolder/1.accels.including.$outname.bw

wigToBigWig $outfolder/2.accels.all_but.$outname.wig -clip /Users/jk2014/Downloads/Moving_stuff_to_laptop/Presentations_documents/hg19.chrom.sizes  $outfolder/2.accels.excluding.$outname.bw
wigToBigWig $outfolder/2.accels.including.$outname.wig -clip /Users/jk2014/Downloads/Moving_stuff_to_laptop/Presentations_documents/hg19.chrom.sizes  $outfolder/2.accels.including.$outname.bw

wigToBigWig $outfolder/3.accels.all_but.$outname.wig -clip /Users/jk2014/Downloads/Moving_stuff_to_laptop/Presentations_documents/hg19.chrom.sizes  $outfolder/3.accels.excluding.$outname.bw
wigToBigWig $outfolder/3.accels.including.$outname.wig -clip /Users/jk2014/Downloads/Moving_stuff_to_laptop/Presentations_documents/hg19.chrom.sizes  $outfolder/3.accels.including.$outname.bw

wigToBigWig $outfolder/4.accels.all_but.$outname.wig -clip /Users/jk2014/Downloads/Moving_stuff_to_laptop/Presentations_documents/hg19.chrom.sizes  $outfolder/4.accels.excluding.$outname.bw
wigToBigWig $outfolder/4.accels.including.$outname.wig -clip /Users/jk2014/Downloads/Moving_stuff_to_laptop/Presentations_documents/hg19.chrom.sizes  $outfolder/4.accels.including.$outname.bw

wigToBigWig $outfolder/4.accels.all_but.$outname.wig -clip /Users/jk2014/Downloads/Moving_stuff_to_laptop/Presentations_documents/hg19.chrom.sizes  $outfolder/4.accels.excluding.$outname.bw
wigToBigWig $outfolder/4.accels.including.$outname.wig -clip /Users/jk2014/Downloads/Moving_stuff_to_laptop/Presentations_documents/hg19.chrom.sizes  $outfolder/4.accels.including.$outname.bw
