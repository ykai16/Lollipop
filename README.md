# Lollipop

Lollipop is a machine-learning-based pipeline for constructing the CTCF-mediated interactome by integrating genetic, epigenetic and gene expression data. In our paper *Predicting CTCF-mediated long-range interactions by integrating genetic, epigenetic and gene expression data*（`link`）, it was used for:

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

The input data used by *Lollipop* are available in `input_data`. For a complete list of used data set, please see Supplementary Methods and Table 1 in the paper.

## Generation of Training Data

Pre-generated training data used in the paper are available in `training_data`. The data format is:

| chrom   |      start1      |  start2 |     response      | ...features...     |
|----------|:-------------:|------:|:-------------:|:-------------:|

You can also generate training data for any cell-type of interest, as long as you have experimental data for CTCF-mediated loops, such as CTCF ChIA-PET and Hi-ChIP data. It takes two steps to do so:

1. Use `prepare_training_interactions.py` to prepare the positive and negative loops for training purpose.

2. Run `add_features.py` to characterize the prepared loops. 

## Making *De Novo* Predictions

Lollipop uses random forest classifier to distinguish positive loops from negative loops. The classifier trained from the three cell-lines (in `.pkl` format) and the *de novo* predictions made by each classicier are available in `denovo_predictions`. The format of predicted loops is:

| chrom   |      start1      |  start2 |     probability      |    yes\_or_no      |  
|----------|:-------------:|------:|:-------------:|:-------------:|

Predicted loops that can be visualized in genome browsers, such as UCSC genome browser, IGV and Washington U genome browser, are also available in the same folder.

We encourage you to apply the trained models to make *de novo* predictions in a cell-type of your interest. After preparing the necessary input data as discussed above, you can run `predict_CTCF_interactions.py` to obtain the predicted loops, in multiple formats that are compatible with genome browsers. 













