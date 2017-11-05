motif=/media/Data/Dropbox/Research/Ongoing/3D_Genome/CTCF_motif/data/CTCF_motif_with_phastCon.txt

bs=/media/Data/Dropbox/Research/Ongoing/3D_Genome/CTCF_motif/data/K562/K562_CTCF_ENCFF002DDJ_renamed_summits.bed
table=/media/Data/Dropbox/Research/Ongoing/Lollipop/data/K562/K562_signals.txt
clf=/media/Data/Dropbox/Research/Ongoing//Lollipop/data/K562/K562_model.joblib.pkl
outdir=/media/Data/Dropbox/Research/Ongoing/3D_Genome/CTCF_motif/results/K562
echo "python predict_CTCF_interactions.py -m $motif -b $bs -t $table -c $clf -o $outdir"
#python predict_CTCF_interactions.py -m $motif -b $bs -t $table -c $clf -o $outdir



bs=/media/Data/Dropbox/Research/Ongoing/3D_Genome/CTCF_motif/data/GM12878/GM12878_CTCF_ENCFF002CDP_renamed_summits.bed
table=/media/Data/Dropbox/Research/Ongoing/Lollipop/data/GM12878/GM12878_signals.txt
clf=/media/Data/Dropbox/Research/Ongoing//Lollipop/data/GM12878/GM12878_model.joblib.pkl
outdir=/media/Data/Dropbox/Research/Ongoing//Lollipop/predictions/GM12878
echo "python predict_CTCF_interactions.py -m $motif -b $bs -t $table -c $clf -o $outdir"
python predict_CTCF_interactions.py -m $motif -b $bs -t $table -c $clf -o $outdir





motif=/media/Data/Dropbox/Research/Ongoing/3D_Genome/CTCF_motif/data/CTCF_motif_with_phastCon.txt
mode=2
clf=/media/Data/Dropbox/Research/Ongoing/3D_Genome/CTCF_motif/data/HeLa/HeLa-S3_CTCF_ENCFF002CDW_renamed_summits.bed
table=/media/Data/Dropbox/Research/Ongoing/3D_Genome/CTCF_motif/data/HeLa/HeLa_signals.txt
predictor=/media/Data/Dropbox/Research/Ongoing/3D_Genome/CTCF_motif/motif_interactions/features/HeLa_RF_model.joblib.pkl
outdir=/media/Data/Dropbox/Research/Ongoing/3D_Genome/CTCF_motif/results/HeLa

#python predict_CTCF_interactions.py -m $motif -i $mode -c $clf -t $table -p $predictor -o $outdir
