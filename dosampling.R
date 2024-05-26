##' dosampling.R, written by Alex Lubbock <code@alexlubbock.com>
##'
##' Copyright (c) The University of Edinburgh, 2016. Licensed under the
##' Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International
##' License. See attached LICENCE file or view online at
##' http://creativecommons.org/licenses/by-nc-sa/4.0

##' This code loads protein expression data from reverse-phase protein
##' arrays (RPPAs) along with patients age, to examine the effect
##' of including differing numbers of tumour samples per patient on
##' the effectiveness of a classifier called NEAT, using N-Cadherin,
##' EpCAM and mTOR expression, together with patients' ages.
##'
##' The sampling regime uses Sobol sampling, a quasi-random technique
##' which ensures low discrepancy, together with a mixed-radix
##' indexing system which allows us to apply the sampling consistently
##' even where patients have differing numbers of samples.

##' Settings
NUM.SAMP.COMBN <- 1e6 # number of sample combinations
WRITE.BUFFER.SIZE <- 1000 # number of iterations before writing to disk

##' NUM.SAMP.COMBN should be integer multiple of WRITE.BUFFER.SIZE
stopifnot(NUM.SAMP.COMBN %% WRITE.BUFFER.SIZE == 0)
writes.per.num.samp <- floor(NUM.SAMP.COMBN / WRITE.BUFFER.SIZE)
total.writes <- writes.per.num.samp * 3

##' Start the timer
ptm <- proc.time()

##' Load libraries and combinadic code
library(survival)
library(plyr)
source("combinadics.R")

##' Load the NEAT model (object name is neat.mdl)
load("NEAT-model.RData")

##' Load expression and patient data
expdata <- read.csv("expressiondata.csv", row.names=1)
patdata <- read.csv("patientdata.csv", row.names=1)

##' Create Survival object from survival times and statuses
patdata <- data.frame(Survival=Surv(time=patdata$Survival.time,
                      event=patdata$Survival.status), Age=patdata$Age,
                      row.names=rownames(patdata))

##' Extract the IDs for samples and patients
sampnm <- rownames(expdata)
patnm <- sub("_[A-Z0-9_/]*$", "", sampnm)
##' Attach the patient ID to the samples data frame
expdata <- data.frame(Patient=patnm, expdata)
##' Table of number of samples per patient
sampfreq <- table(patnm)

for(num.samp in 1:3) {
    ##' Patients we can sample from (as opposed to use as-is,
    ##' because they have too few samples)
    patients.to.sample <- names(which(sampfreq > num.samp))
    ##' Number of sample combinations for each patient
    samp.combn <- choose(sampfreq[patients.to.sample], num.samp)

    buffer <- NULL
    writes <- 0
    ##' Sobol facilitates even coverage of possible combinations,
    ##' which we can refer to uniquely with an identifier called
    ##' a combinadic
    combind <- floor(sobol(NUM.SAMP.COMBN) * prod(samp.combn))
    for(i in combind) {
        ##' Include all samples from patients with <= max # samples
        samples.to.include <- patnm %in% names(which(sampfreq <= num.samp))
        ##' Generate a combinadic for each patient, for their sample combination,
        ##' derived from the master combinadic using mixed radix arithmetic
        indexes <- as.mixed.radix(samp.combn, i - 1) + 1
        if(length(patients.to.sample)>0) {
            for(ind in 1:length(patients.to.sample)) {
                smp.idx <- combinadic(samp.combn[ind], num.samp, indexes[ind]) # list sample indices
                smp <- sampnm[patnm == patients.to.sample[ind]][smp.idx] # the actual sample ids
                samples.to.include <- samples.to.include | sampnm %in% smp
            }
        }
        ##' build the data frame
        sub.dat <- data.frame(ddply(expdata[samples.to.include,], ~Patient, summarise, mTOR_med=median(mTOR), EpCAM_med=median(EpCAM), N_Cad_med=median(N_Cad))[,-1], patdata[,c("Age","Survival")])
        
        ##' run the test
        mdl.pred <- predict(neat.mdl, sub.dat)
        pval <- pchisq(survdiff(Survival~mdl.pred>0, data=sub.dat)$chisq, 1, lower.tail=FALSE)
        hr <- exp(coef(coxph(Survival~mdl.pred>0, data=sub.dat)))
        buffer <- rbind(buffer, list(Combinadic=i, HR=hr, LRpval=pval))
        if(nrow(buffer) == WRITE.BUFFER.SIZE) {
            writes <- writes + 1
            write.table(buffer, file=paste0(num.samp,"-sample-HR.tsv"), append=ifelse(writes == 1, FALSE, TRUE), quote=FALSE, col.names=FALSE, row.names=FALSE)
            percent.complete <- (((num.samp - 1) * writes.per.num.samp + writes) / total.writes) * 100
            cat(sprintf("%.2f%% complete; elapsed time=%.1fs\n", percent.complete, (proc.time() - ptm)[3]))
            buffer <- NULL
        }
    }

    ##' Write buffer should be empty
    stopifnot(is.null(buffer))
}
