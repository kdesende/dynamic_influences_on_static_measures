# dynamic_influences_on_static_measures
This repo contains code for the paper "Dynamic influences on static measures of metacognition" authored by Kobe Desender, Luc Vermeylen, and Tom Verguts.

#Experiment code:
RDK_SATS contains the code that was used to run Experiment 1. Due to covid-19 this study was run online (using jspych). 

#Analysis, 5 R files that implement the following: 
1. DDM Simulations demontrating how meta-d'/d' (=mratio) depends on decision boundary
2. Parameter recovery study on the evidence accumulation model with post-decision accumulation
3A. Fitting an evidence accumulation model with post-decision accumulation to Experiment 1 (explicit SATO instructions, with 6-choice confidence ratings) to estimate m-ratio and v-ratio ("data_exp1_categoricalconf.csv")
3B. Replication of Experiment 1 with continuous confidence ratings ("data_exp1.csv")
4. Fitting the same model to Experiment 2A (data of Drescher et al.)
5. Fitting the same model to Experiment 2C (data of Prieto et al.)

#Data:
The data of Experiment 1 can be found in the file "data_exp1_categoricalconf.csv"
The data of Experiment 1S reported in the Supplementary Materials can be found in the file "data_exp1.csv"
The data of Experiment 2A (Drescher et al., 2018, PLoS ONE), can be found here: https://github.com/l-drescher/raw_data_MW_MC_CC
The data of Experiment 2B can be found in the file "data_exp2B.csv"
The data of Experiment 2C (Prieto et al., in preparation), can be found in the Confidence Database: https://osf.io/s46pr/
