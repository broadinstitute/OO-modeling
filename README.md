# OO-modeling
Welcome to the code and data repository for Analyzing the Impact of a Real-life Outbreak Simulator on Pandemic Mitigation: an Epidemiological Modeling Study.

A quick overview of the code present:
- Small_Model_CMU.Rmd and Small_Model_BYU.Rmd simulate the spread of a pathogen through their respective OO contact networks. The output from these scripts are stored in final_CMU.Rdata and final_CMU.Rdata, where it feautures in the main regression analysis used to generate FIgure 5.
- Model_BYU.Rmd contains the code used for the stochastic epidemiological model of the full college campus, as well as scripts for all figures pertaining to BYU.
- Model_CMU.Rmd contains the scripts for all figures pertaining to CMU.
- histories_CMU.csv and histories_BYU.csv contain data on all events recorded by OO at each campus, respectively.
- leger_CMU and leger_BYU match user ID's and P2P kit codes for each simulation, respectively.
- out (folder) contains the graphics output by Model_CMU.Rmd and Model_BYU.Rmd.

