 # REFORMATTING OF TRAINING DATA FOR INPUT INTO MOFA
# -------------------------------------------------

# This data has already been normalised by the CMI-PB team. 
# We will renormalise the gene expression data using vst 

## LOAD DATA
## ---------

source('../scripts/libs.R')


# load data provided by CMI-PB3
dat = readRDS('../../../../../data//processed_datasets/master_allData_batchCorrected.RDS')
names(dat)

# processed data doesn't contain DOB info, let's get this from some other CMI-PB3 data
d2 = readRDS('../../../../data/processed_datasets/master_harmonized_data.RDS')

## ADD AGE DATA TO METADATA
## ------------------------

# get DOB and convert to age (for both challenge and training data)
dob = c(d2$training$subject_specimen$year_of_birth,
  d2$challenge$subject_specimen$year_of_birth)
boost_date = c(d2$training$subject_specimen$date_of_boost, 
        d2$challenge$subject_specimen$date_of_boost)
age = as.numeric(difftime(boost_date,  dob, units = 'weeks')/52)
age

# add age to meta data 
id = c(d2$training$subject_specimen$specimen_id, d2$challenge$subject_specimen$specimen_id)

setdiff(dat$subject_specimen$specimen_id, id)
setdiff(id, dat$subject_specimen$specimen_id)

df = data.frame(specimen_id = id, age = age)
head(df)

dat$subject_specimen = left_join(dat$subject_specimen, df)

## SELECT DATASETS TO USE
## ----------------------

# in this case we will use all datasets together, and separate them later
ds=c('2020_dataset', '2021_dataset', '2022_dataset', '2023_dataset')
dat$subject_specimen = dat$subject_specimen[dat$subject_specimen$dataset %in% ds, ]
table(dat$subject_specimen$dataset)
dim(dat$subject_specimen)

# we are having issues with R not knowing whether to treat sample IDs as numbers or characters
# let's convert to characters to be sure
dat$subject_specimen$specimen_id = paste('x', dat$subject_specimen$specimen_id, sep='')
dat$subject_specimen$subject_id = paste('x', dat$subject_specimen$subject_id, sep='')
head(dat$subject_specimen)

## MERGE DATA
## -----------


# MOFA requires the same number of samples for all assays it analyses 
# merge all matrices to replace empty values with NAs

meta = dat$subject_specimen
colnames(meta) = paste("Meta.", colnames(meta), sep='')
dim(meta)
head(meta)

# ab
x = as.data.frame(t(dat$plasma_ab_titer$batchCorrected_data))
dim(x)
rownames(x) = paste('x', rownames(x), sep='')
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))

# cytokines
x = as.data.frame(t(dat$plasma_cytokine_concentrations_by_olink$batchCorrected_data))
colnames(x) = paste("CytokineO.", colnames(x), sep='')
rownames(x) = paste('x', rownames(x), sep='')
head(x)
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))


# cytokine legendplex
x = as.data.frame(t(dat$plasma_cytokine_concentrations_by_legendplex$normalized_data))

colnames(x) = paste("CytokineL.", colnames(x), sep='')
rownames(x) = paste('x', rownames(x), sep='')
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))

# T cell activation
x = as.data.frame(t(dat$t_cell_activation$raw_data))
colnames(x) = paste("TcellActivation.", colnames(x), sep='')
head(x)
rownames(x) = paste('x', rownames(x), sep='')
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))

# t cell polarisation
# take pre-processed data provided by CMI-PB team
tcellpol = read_tsv('../../../../data/processed_datasets/Th1(IFN-γ)_Th2(IL-5)_polarization_ratio.tsv')
tcellpol$specimen_id = as.character(tcellpol$specimen_id)
head(tcellpol)

x = tcellpol[,c(1,2)]
head(x)
x$specimen_id = paste('x', x$specimen_id, sep='')
colnames(x) = c("TcellPol.specimen_id", "TcellPol.pol")
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "TcellPol.specimen_id"))


# cell_freq
x = as.data.frame(t(dat$pbmc_cell_frequency$batchCorrected_data))
colnames(x) = paste("Freq.", colnames(x), sep='')
rownames(x) = paste('x', rownames(x), sep='')
head(x)
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))

# expression
# note: VST followed by batch correction
counts = dat$pbmc_gene_expression$raw_count$raw_data
colnames(counts) = paste('x', colnames(counts), sep='')
counts = counts[, which(colnames(counts) %in% meta$Meta.specimen_id)]
counts = apply(counts,2, as.integer)
rownames(counts) = rownames( dat$pbmc_gene_expression$raw_count$raw_data)
dim(counts)
counts[1:5, 1:5]

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = dat$subject_specimen[dat$subject_specimen$specimen_id %in% colnames(counts) , ],
                              design = ~ as.factor(infancy_vac))
dds
dds = estimateSizeFactors(dds)
vsd <- vst(dds, blind=TRUE)
dim(assay(vsd))

vsd$subject_id = as.factor(vsd$subject_id)

bef_noHVF1 = plotPCA(vsd, "dataset") # PCA before
bef_noHVF2 = plotPCA(vsd, "subject_id") # PCA before
bef_noHVF1

assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch=as.vector(dds$dataset))
aft_noHVF1 = plotPCA(vsd, "dataset")
aft_noHVF2 = plotPCA(vsd, "subject_id")
aft_noHVF1


x = as.data.frame(t(assay(vsd)))
colnames(x) = paste("Expr.", colnames(x), sep='')
x[1:5, 1:5]
meta = left_join(meta, rownames_to_column(x), by=c(  "Meta.specimen_id" = "rowname"))

# separate DF into list again
dat = list(meta = as.data.frame(meta[,grep('Meta', colnames(meta))]),
                  geneExp = as.data.frame(t(meta[,grep('Expr.', colnames(meta))])) ,
                  ab = as.data.frame(t(meta[,grep('^IgG', colnames(meta))])) ,
                  cytokineO = as.data.frame(t(meta[,grep('CytokineO', colnames(meta))])), 
                  cytokineL = as.data.frame(t(meta[,grep('CytokineL', colnames(meta))])), 
                  tcellPol = as.data.frame(t(meta[,grep('TcellPol', colnames(meta))])),
                  tcellAct = as.data.frame(t(meta[,grep('TcellAct', colnames(meta))])), 
                  cell_freq = as.data.frame(t(meta[,grep('Freq', colnames(meta))])
))

colnames(dat$geneExp) = dat$meta$Meta.specimen_id
colnames(dat$ab) = dat$meta$Meta.specimen_id
colnames(dat$cytokineO) = dat$meta$Meta.specimen_id
colnames(dat$cytokineL) = dat$meta$Meta.specimen_id
colnames(dat$cell_freq) = dat$meta$Meta.specimen_id
colnames(dat$tcellPol) = dat$meta$Meta.specimen_id
colnames(dat$tcellAct) = dat$meta$Meta.specimen_id

dat$ab[1:5, 1:13]

lapply(dat, dim)

# split into test and train
# 35-65%

idx = dat$meta %>%
  filter(Meta.dataset %in% c('2022_dataset', '2021_dataset', '2020_dataset')) %>%
  select('Meta.subject_id')  %>%
  unique() 
idx = idx$Meta.subject_id
length(idx)

set.seed(426)
idx.sample = sample(idx, size = round(length(idx)*0.65), replace = F)
idx.test = idx[!idx %in% idx.sample]


train.meta = dat$meta %>%
  filter(Meta.subject_id %in% idx.sample)

test.meta = dat$meta %>%
  filter(Meta.subject_id %in% idx.test)


dim(train.meta)
dim(test.meta)

# keep challenge data separate
challenge.idx = dat$meta[dat$meta$Meta.dataset %in% c('2023_dataset'), 'Meta.specimen_id']
length(challenge.idx)

train_data = list(
  meta = train.meta,
  geneExp = dat$geneExp[, train.meta$Meta.specimen_id], 
  ab = dat$ab[, train.meta$Meta.specimen_id], 
  cytokineO = dat$cytokineO[, train.meta$Meta.specimen_id],
  cytokineL = dat$cytokineL[, train.meta$Meta.specimen_id],
  cell_freq = dat$cell_freq[, train.meta$Meta.specimen_id],
  tcellPol =dat$tcellPol[, train.meta$Meta.specimen_id],
  tcellAct = dat$tcellAct[, train.meta$Meta.specimen_id]
  
)
train_data$cytokineO
train_data$tcellPol
lapply(train_data, dim)
lapply(train_data, class)

test_data = list(
  meta = test.meta,
  geneExp = dat$geneExp[, test.meta$Meta.specimen_id], 
  ab = dat$ab[, test.meta$Meta.specimen_id], 
  cytokineO = dat$cytokineO[, test.meta$Meta.specimen_id],
  cytokineL = dat$cytokineL[, test.meta$Meta.specimen_id],
  cell_freq = dat$cell_freq[, test.meta$Meta.specimen_id],
  tcellPol =dat$tcellPol[, test.meta$Meta.specimen_id],
  tcellAct = dat$tcellAct[, test.meta$Meta.specimen_id]
)
test_data$ab[1:5, 1:13]
test_data$tcellPol
lapply(test_data, dim)

challenge_data = list(
  meta = dat$meta[dat$meta$Meta.specimen_id %in% challenge.idx, ],
  geneExp = dat$geneExp[, challenge.idx], 
  ab = dat$ab[, challenge.idx], 
  cytokineO = dat$cytokineO[, challenge.idx],
  cytokineL = dat$cytokineL[, challenge.idx],
  cell_freq = dat$cell_freq[, challenge.idx],
  tcellPol =dat$tcellPol[, challenge.idx],
  tcellAct = dat$tcellAct[, challenge.idx]
)
challenge_data$ab[1:5, 1:13]
challenge_data$tcellPol
lapply(challenge_data, dim)



save(train_data, file='../data/split1/TRAIN.RData')
save(test_data, file='../data/split1/TEST.RData')
save(challenge_data, file='../data/split1/CHALLENGE.RData')

