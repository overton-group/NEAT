##' combinadics.R, written by Alex Lubbock <code@alexlubbock.com>

##' Copyright (c) The University of Edinburgh, 2016. Licensed under the
##' Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International
##' License. See attached LICENCE file or view online at
##' http://creativecommons.org/licenses/by-nc-sa/4.0

##' Combinadics are individual combination instances of n choose k. R provides 
##' the function combn for listing these in their entirety, but this soon gets
##' intractable for large n. This code can be used for for iterating over all
##' instances or, more likely, selecting a subset of the instances for sampling.

##' Uses the gmp library to provide bigint support, avoiding integer overflow

library(gmp) #' for arbitrary precision maths
library(randtoolbox) #' for sobol()

##' Returns the m1-position combinadic from
##' n choose k. A combinadic indicates the
##' specific elements to choose for a
##' combination. A vector of length k is
##' returned which specifies the indexes of
##' the selected elements.
combinadic <- function(n, k, m1) {
    mychoose <- choose
    if(n > 100) mychoose <- chooseZ

    m <- m1 - 1
    nck <- mychoose(n, k)
    
    if(nck == 0) {
        stop(n, " choose ", k, "is zero! No combinadic is possible.")
    }
    if(m < 0 || m >= nck) {
        stop("m1 must be between 0 and ", nck)
    }

    ans <- mat.or.vec(k, 1)
    x <- nck - m - 1 # x is the dual of i

    a <- n
    b <- k
    for(i in (1:k) - 1) {
        a <- a - 1
        while(mychoose(a, b) > x) a <- a - 1
        x <- x - mychoose(a, b)
        b <- b - 1
        ans[i + 1] <- n - a
    }

    ans
}

##' Takes a set of scalar bases and converts a single integer
##' value into 'buckets' for each of these radices.
##'
##' Useful for converting time e.g. from seconds to hh:mm:ss
##' But here mainly for sampling from buckets with mixed size
##' using a single master sample index (the input value)
##'
##' Make sure the input value scale starts from 0, not 1!
as.mixed.radix <- function(bases, input) {
    output <- c()
    for(b in rev(bases)) {
        output <- c(input %% b, output)
        input <- input %/% b
    }
    if(input!=0) warning("Overflow - need more radices!")
    if(output[length(output)] %% 1 != 0) {
        warning("Input is not an integer - rounding down")
        output[length(output)] <- floor(output[length(output)])
    }
    output
}

##' Samples the i'th combinadic of k elements from x
##' Allows deterministic, low discrepancy sampling from
##' any vector.
sample.combinadic <- function(x, k, i) {
    stopifnot(k <= length(x))
    numcom <- bigchoose(length(x), k)
    comb <- sobol(i)[i] * numcom
    x[combinadic(length(x), k, comb)]
}
