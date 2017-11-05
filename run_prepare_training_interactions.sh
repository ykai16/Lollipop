chip=/media/Data/data/hg19/ChIA-PET/HeLa/ChIA-PET2/GSE72816_HeLa_summits.bed
#chip=/media/Data/data/hg19/ENCODE_CTCF/peak/HeLa_ENCFF002CDW.bed
chia=/media/Data/Dropbox/Research/Ongoing/Lollipop2.0/data/HeLa/HeLa_MICC_intrachrom_loops.bed
hic=/media/Data/Dropbox/Research/Ongoing/Lollipop2.0/data/HeLa/HeLa_HIC_loops.txt
outfile=/media/Data/Dropbox/Research/Ongoing/Lollipop2.0/data/HeLa/HeLa_training_interactions.txt
echo "python prepare_training_interactions.py -p $chip -a $chia -c $hic -o $outfile"
#python prepare_training_interactions.py -p $chip -a $chia -c $hic -o $outfile


chip=/media/Data/data/hg19/ChIA-PET/GM12878/output_1.0/GSE72816_GM12878_summits.bed
#chia=/media/Data/Dropbox/Research/Ongoing/Lollipop2.0/data/GM12878/GM12878_MICC_intrachrom_loops.bed
chia=/media/Data/Dropbox/Research/Ongoing/Lollipop2.0/data/GM12878/GM12878_MICC_intrachrom_loops_0.15.bed
hic=/media/Data/Dropbox/Research/Ongoing/Lollipop2.0/data/GM12878/GM12878_HIC_loops.txt
outfile=/media/Data/Dropbox/Research/Ongoing/Lollipop2.0/data/GM12878/GM12878_training_interactions_0.15.txt
echo "python prepare_training_interactions.py -p $chip -a $chia -c $hic -o $outfile"
python prepare_training_interactions.py -p $chip -a $chia -c $hic -o $outfile

chip=/media/Data/data/hg19/ChIA-PET/K562/ChIA-PET2/K562_CTCF_summits.bed
chia=/media/Data/Dropbox/Research/Ongoing/Lollipop2.0/data/K562/K562_MICC_intrachrom_loops.bed
hic=/media/Data/Dropbox/Research/Ongoing/Lollipop2.0/data/K562/K562_HIC_loops.txt
outfile=/media/Data/Dropbox/Research/Ongoing/Lollipop2.0/data/K562/K562_training_interactions.txt
echo "python prepare_training_interactions.py -p $chip -a $chia -c $hic -o $outfile"
#python prepare_training_interactions.py -p $chip -a $chia -c $hic -o $outfile
