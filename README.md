NEAT performance on datasets using limited tumour samples per patient
=====================================================================

Code written by Alex Lubbock <code@alexlubbock.com> in Ian Overton's
group <ian.overton@ed.ac.uk>. Licensed under CC-BY-NC-SA v4 (please
see LICENSE.txt also available at: 
http://creativecommons.org/licenses/by-nc-sa/4.0/legalcode)

Introduction
------------

This code loads protein expression data from reverse-phase protein
arrays along with patient age, and examines the effect
of including differing numbers of tumour samples per patient on
validation of the 'NEAT' algorithm (N-cadherin, EPCAM, Age, mTOR).

The sampling regime uses Sobol sampling, a quasi-random technique
which ensures low discrepancy, together with a mixed radix
indexing system which allows us to apply the sampling consistently
even where patients have differing numbers of samples.

One master combinadic (combination identifier) is generated for each
sample combination. For a given dataset, this is then mapped to a
set of combinadics, one for each patient being sampled. Finally,
this identifier is mapped onto the set of samples to include in
the dataset under consideration.

The approach has the advantages of sampling the available range of
samples evenly (Sobol sequences are quasi-random, low discrepancy),
deterministically and with traceability (from the master combinadic
and the raw data, we can easily determine which patients and
which samples were included). The deterministic nature can also
support parallel execution with guaranteed non-overlapping
sample sets.

Requirements
------------

This code has been tested on R 3.2.2 using Linux and Mac OSX, but
is expected to work on Windows and other recent versions of
R.

On our test machine (2014 Macbook Pro with Core i7 "Haswell"
2.5GHz), memory usage (RAM) stayed under 200Mb. The output files
require around 150Mb of free disk space. Execution time was
approximately 9 hours.

Instructions
------------

The following R packages need to be installed:

 * randtoolbox
 * gmp
 * plyr
 * survival

These can be installed from within R with the following command:

    install.packages(c("randtoolbox", "gmp", "plyr", "survival"))

The command below runs the analysis:

    Rscript dosampling.R

Three files will be output: 1-sample-HR.tsv,
2-sample-HR.tsv and 3-sample-HR.tsv. These files contain NEAT model
predictive results using 1, 2 or 3 tumour samples per patient
respectively.

Each file contains a three column table with one million rows.
Each row is the NEAT model predictive result for a given combinadic
- a particular combination of tumour samples from the validation cohort. 
The three columns are: 
1. Combinadic ID (a unique identifier for the combinadic)
2. Log hazard ratio for NEAT on the samples in that combinadic
3. Log-rank p-value

To override the default of one million samples per number of tumour
samples, open the dosampling.R script and change the NUM.SAMP.COMBN
variable near the top. Please note that the NUM.SAMP.COMBN value 
must be an integer multiple of 1000.

*If the script is re-run, output files will be overwritten
without warning.*

List of Files
-------------

NEAT-model.RData - An RData object containing the NEAT model itself ("coxph" class object from R's "survival" library)
combinadics.R - R functions for combinadic and mixed radix arithmetic
dosampling.R - R script for performing the tumour sampling procedure described above
expressiondata.csv - Anonymised per-tumour sample RPPA expression data for the three NEAT protein markers (N-Cadherin, EPCAM and mTOR)
patientdata.csv - Anonymised ages, survival time and survival status
README.md - This README file
LICENSE.txt - The CC-BY-NC-SA license text
