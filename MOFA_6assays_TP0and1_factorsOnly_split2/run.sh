#!/bin/bash


cd /mnt/bioadhoc/Groups/Peters/nthrupp/CMI-PB3/FINAL-test/submissions/models_splitTestTrain/MOFA_6assays_TP0and1_factorsOnly_split2
pwd

set -e
set -o pipefail

echo 'Creating data directories'

mkdir -p data/MOFA_models
mkdir -p data/split1
mkdir -p data/regression_models
mkdir -p data/predictions


echo 'Pre-processing'
cd pre-processing


Rscript prep_train_data.R &> prep_train_data.out

cd ../

echo 'MOFA'
cd MOFA

Rscript -e "rmarkdown::render('first_pass_noScale_train.Rmd')"

Rscript prep_testData.R &>  prep_testData.out
Rscript prep_trainData.R &>  prep_trainData.out

cd ../

echo 'Regression models'
cd regression_models

Rscript predict_CCL3.R &> predict_CCL3.out

cd ../

echo 'Predictions'
cd predict_2023

Rscript prep_challengeData.R &> Rscript prep_challengeData.out
Rscript CCL3_imputeMedian.R &> CCL3_imputeMedian.out

cd ../

echo 'RUN COMPLETE'
