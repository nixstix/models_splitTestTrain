source('../scripts/libs.R')
source('../scripts/run_glmnet.R')
detach("package:MOFA2", unload=TRUE) # need to detach MOFA in order to use the predict function from glmnet

load('../data/MOFA_models/predictors.RData')
load('../data/MOFA_models/test.RData')


# test predictions on 2022 data
# -----------------------------
# -----------------------------


task = 'tcellPol'

## baseline model
## --------------

# use demographic data, as well as assay data that should be highly predictive of the prediction task



features = c("Day0.TcellPol.pol", 'Day0.Meta.age')

a = predictors.baseline[,c(features)]
summary(a)

b = test_data.baseline[,c(features)]
summary(b)

model1 = run_glmnet(p = features, task.name = task, alpha=1)
model1$model_coef
model1$model_cor


save(model1,  file='../data/regression_models/regression_models_TcellPol_bl.RData')
