# MetabolicSurv
R package : A biomarker validation approach for predicting survival using metabolic signature, this package develope biomarker signature for metabolic data. It contains a set of functions and cross validation methods  to validate and select biomarkers when the outcome of interest is survival. The package can handle prognostic factors and mainly metabolite matrix as input, the package can served as biomarker validation tool.

## Why use the package
* It can be used with any form of high dimensional/omics data such as: Metabolic data, Gene expression matrix, incase you dont have a data it can simulate hypothetical scinerio of a high dimensional data based on the desired biological parameters
* It developed any form of signature from the high dimensional data to be used for other purpose
* It also employs data reduction techniques such as PCA, PLS and Lasso 
* It classifies subjects based on the signatures into Low and high risk group
* It incorporate the use of subject prognostic information for the to enhance the biomarker for classification
* It gives information about the surival rate of subjects depending on the classification



## Installation

You can install the released version of MetabolicSurv from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("MetabolicSurv")

```
## Illustrations to simulate a Metabolomic profile matrix
Apart from the survival prediction and classification, \pkg{MetabolicSurv} can also be used to generate an artificial Metabolomic profile matrix, survival data (Survival time and censoring indiicator) and clinical covariates which will be referred to as prognostic factors to be used for further analysis or for other pursoses. Since there a few publicly available metabolic profile matrix this package can be used to firstly simulate each of this respective dataset which is required to evaluate the other basic and advance function in the package.


``` r
	library(MetabolicSurv)
	Data <- MSData(nPatients = 200, nMet = 3000, Prop = 0.5)
	Metdata <- Data$Mdata
	Survdata <- Data$Survival
	Censordata <- Data$Censor
	Progdata <- Data$Prognostic
	
``` 

The code above was used to simulate a metabolomic, survival and prognostic data with a total of 200 patients with 3000  metabolites in the metabolomic profile matriix  assuming that the proportion of patients having low risk is 0.5 . The proportion can be adjusted depending on how strict one need to be in assuming equal or unequal proportion of classification based on biological findings or intelligent guess. The Metabolomic profile matrix is stored in Metdata, the survival time is stored in Survdata, Censoring information in Censordata and the Prognosticfactor/clinical covariates in Progdata.


## A quick Demostration to solve a problem

``` r
"Problem of interest"

"Given a set of subjects with known riskscores and prognostic features how can we use this information to obtain their risk of surving and what group does each respective subject belongs to?"

```

``` r
##  Loading the package
library("MetabolicSurv")

##  Loading one of the inbuilt data
data(DataHR)
names(DataHR)

##  This function does Classification, Survival Estimation and Visualization
Result = EstimateHR(Risk.Scores=DataHR[,1],Data.Survival=DataHR[,2:3]
,Prognostic=DataHR[,4:5],Plots=FALSE,Quantile=0.50)

## Survival information
Result$SurvResult


## Group information
Result$Riskgroup
```

## Functions in the package

| Category	|	Functions	|	Description	                                |
| --------- | --------- | ------------------------------------------- |
| Basic	|	MSpecificCoxPh	|	Metabolite by metabolite Cox proportional hazard analysis
|  | SurvPcaClass	| Classifier based on first PCA
| | SurvPlsClass| Classifier based on first PLS 
| | Majorityvotes |Classifiction for Majority Votes	
| | Lasoelacox | Wrapper function for glmnet
| | MSData | Generate Artificial Metabolic Survival Data
| Advance	| CVLasoelacox |Cross Validations for Lasso Elastic Net predictive models and Classification
|  | CVSim	| Cross-validation for Top $K_{1}, \ldots, K_{n}$ metabolites
| | CVPcaPls |	Cross-validations for PCA and PLS based methods
| | CvMajorityvotes | Cross-validation for majority votes	
| | MetFreq |Frequency of Selected Metabolites from the Metabolite specific Cross Validation
| | QuantileAnalysis |Sensitivity of the quantile used for classification
| | Icvlasoel |	Inner and outer cross-validations for shrinkage methods
| | DistHR | Null distribution of the estimated HR
| | SIMet | Sequentially increase the number of top $K$ metabolites
