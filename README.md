# dynamic_influences_on_static_measures
This repo contains code for the paper "Dynamic influences on static measures of metacognition" authored by Kobe Desender, Luc Vermeylen, and Tom Verguts.

#Experiment code:
RDK_SATS contains the code that was used to run Experiment 1. Due to covid-19 this study was run online (using jspych). 

#Analysis, 5 R files that implement the following: 
1. DDM Simulations demontrating how meta-d'/d' (=mratio) depends on decision boundary
2. Parameter recovery study on the evidence accumulation model with post-decision accumulation
3A. Fitting an evidence accumulation model with post-decision accumulation to Experiment 1 (explicit SATO instructions) to estimate m-ratio and v-ratio
3B. Fitting Experiment 1 using a Pleskac & Busemeyer version of the model in which confidence only depends on post-decisional evidence 
4A. Fitting the same model to Experiment 2 (data of Drescher et al.)
4B. Fitting Experiment 2 using a Pleskac & Busemeyer version of the model in which confidence only depends on post-decisional evidence 
5A. Fitting the same model to Experiment 3 (data of Prieto et al.)
5B. Fitting Experiment 3 using a Pleskac & Busemeyer version of the model in which confidence only depends on post-decisional evidence 

#Data:
The data of Experiment 1 can be found in the file "data_exp1.csv"
The data of Experiment 2 (Drescher et al., 2018, PLoS ONE), can be found here: https://github.com/l-drescher/raw_data_MW_MC_CC
The data of Experiment 3 (Prieto et al., in preparation), can be found in the Confidence Database: https://osf.io/s46pr/
