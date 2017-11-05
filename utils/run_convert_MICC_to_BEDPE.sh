MICC=data/GSE72816_HeLa.interactions.MICC
cutoff=0
fdr=1
outfile=/media/Data/Dropbox/Research/Ongoing/CLIP/data/HeLa/HeLa_MICC_intrachrom_loops.bed

echo "python convert_MICC_to_BEDPE.py -i $MICC -c $cutoff -f $fdr -o $outfile"
python convert_MICC_to_BEDPE.py -i $MICC -c $cutoff -f $fdr -o $outfile


MICC=data/GSE72816_GM12878.interactions.MICC
cutoff=0
fdr=1
outfile=/media/Data/Dropbox/Research/Ongoing/CLIP/data/GM12878/GM12878_MICC_intrachrom_loops.bed

echo "python convert_MICC_to_BEDPE.py -i $MICC -c $cutoff -f $fdr -o $outfile"
python convert_MICC_to_BEDPE.py -i $MICC -c $cutoff -f $fdr -o $outfile

MICC=data/K562_CTCF.interactions.MICC
cutoff=0
fdr=1
outfile=/media/Data/Dropbox/Research/Ongoing/CLIP/data/K562/K562_MICC_intrachrom_loops.bed

echo "python convert_MICC_to_BEDPE.py -i $MICC -c $cutoff -f $fdr -o $outfile"
python convert_MICC_to_BEDPE.py -i $MICC -c $cutoff -f $fdr -o $outfile
