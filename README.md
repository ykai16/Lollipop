# Lollipop

Lollipop is a machine-learning-based framework for predicting the CTCF-mediated interactome by integrating genetic, epigenetic and gene expression data. In our paper *Predicting CTCF-mediated chromatin interactions by integrating genomic and epigenomic features*（`https://www.biorxiv.org/content/early/2017/12/01/215871`）, it was used for:

* Creating positive and negative training data.
* Training a model that distinguishes positive loops from negative loops.
* Applying the trained model to a cell-type of interest to make *de novo* predictions of CTCF-mediated loops. 

## Dependencies
Lollipop requires the following packages:


* Numpy `http://www.numpy.org`
* Pandas `http://pandas.pydata.org`
* Scikit-learn `http://scikit-learn.org/stable/`
* HTSeq `https://htseq.readthedocs.io/en/release_0.9.1/`

We recommend to use [Anaconda python distribution](https://www.anaconda.com/what-is-anaconda/) for installation of the above packages.

## Input Data

The input data used by *Lollipop* are available in `input_data`. For a summary of used data set, please see Supplementary Methods and Table 1 in the paper.

## Generating Training Data

Pre-generated training data used in the paper are available in `training_data`. The data format is:

| chrom   |      start1      |  start2 |     response      | ...features...     |
|----------|:-------------:|------:|:-------------:|:-------------:|

One can also generate training data for any cell-type of interest, as long as experimental data for CTCF-mediated loops, such as CTCF ChIA-PET and Hi-ChIP data, are available. It takes two steps to do so:

### Step 1. Preparing positive and negative loops for training purpose.

Usage:

`python prepare_training_interactions.py -p $CTCF_peak -a $CTCF_ChIA-PET_interactions -c $CTCF_HiC_interactions -o $training_interactions`

Parameters:  

`-p $CTCF_peak:`CTCF peak file in BED format.

`-a $CTCF_ChIA-PET_interactions:`CTCF-mediated interactions identified by ChIA-PET or other methods. The file format is `chrom1 start1 end1 chrom2 start2 end2 IAB FDR strand1 strand2`, where `IAB` is the number of PETs connecting the anchors and `FDR` is the statistical significance.

`-o $training_interactions:`Output file with positive and negative loops in the following format: `chrom anchor1 anchor2 response loop-length`, where `anchor1/2` is the genomic coordinate of the middle point of left/right anchor.

### Step2. Characterizing prepared loops.

Usage: 

`python add_features.py -i $training_interactions -t $information_table -o $training_data`

Parameters:

`-i $training_interactions:` Output file from step1.

`-t $information_table:` A table containing the paths of genomic and epigenomic datasets to derive features. An example of this table can be seen in `input_data`. 

`-o $training_data:` Output file with positive and negative loops characterized by a set of features. 


## Building a model
A model can be generated from the prepared training data, by using `train_model.py`.

Usage:

`python train_model.py -t $training_data -o $output_folder`

Parameters:

`-t $training_data:`The file with positive and negative loops characterized by features.

`-o $output_folder:`The path of the folder where you want to put the resulting model and cross-validation results. ROC and PR curves are generated.

## Making *De Novo* Predictions

Lollipop employs a random forest classifier to distinguish positive from negative loops. The classifier trained from three cell-lines (in `.pkl` format) and the *de novo* predictions made by each classicier are available in `denovo_predictions`. The format of predicted loops is:

| chrom   |      start1      |  start2 |     probability      |    yes\_or_no      |  
|----------|:-------------:|------:|:-------------:|:-------------:|

Predicted loops that can be visualized in genome browsers, including UCSC genome browser, IGV and Washington U genome browser, are also available in the same folder.

One can also apply the trained models to make *de novo* predictions in a cell-type of interest by running `make_denovo_predictions.py`.

Usage:

`python make_denovo_predictions.py -b $CTCF_Peaks -t $information_table -c $classifier -o $output_folder`

Parameters:

`-b $CTCF_Peaks:`CTCF peak file in BED format.

`-t $information_table:` A table containing the paths of genomic and epigenomic datasets to derive features.

`-c $classifier:`The trained classifier used for making predictions.

`-o $output_folder:`The output folder for results. Output files include predicted loops in different formats that can be visulized in genome browser.















