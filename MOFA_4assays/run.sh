#!/bin/bash

cd /mnt/BioAdHoc/Groups/Peters/nthrupp/CMI-PB3/sampleTestTrain/MOFA_4assays
pwd

set -e
set -o pipefail

echo 'Creating data directories'

mkdir -p data/MOFA_models
mkdir -p data/split1 
mkdir -p data/regression_models

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

Rscript predict_TcellPol_bl.R &> predict_TcellPol_bl.out 


Rscript predict_CCL3.R &> predict_CCL3.out
Rscript predict_IgG_PT.R &> predict_IgG_PT.out 
Rscript predict_Mn.R &> predict_Mn.out 
Rscript predict_TcellPol.R &> predict_TcellPol.out 
 

Rscript predict_IgG_PT_1.R &> predict_IgG_PT_1.out 
Rscript predict_IgG_PT_2.R &> predict_IgG_PT_2.out 

cd ../

echo 'RUN COMPLETE'
