//
//  StochHMM_usage.h
//  StochHMM
//
//  Created by Paul Lott on 2/12/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#ifndef StochHMM_StochHMM_usage_h
#define StochHMM_StochHMM_usage_h

const char usage[]  = "\n\
StochHMM - Stochastic Hidden Markov Model Framework (version 0.37) Date: December 12, 2013\n\
\n\
Usage: StochHMM -model <model_file>  -seq <seq_file> [options]\n\
\n\
Command Line options:\n\
\t-help or -h		print usage statement\n\
\n\
Files: reqires a sequence file and a model file\n\
\t-model <model file>\t\timport model file\n\
\t-seq <sequence file>\t\t\timport sequence file in fasta format\n\
\n\
Non-stochastic Decoding:  Different algorithms available for decoding\n\
\t-viterbi\t\t\tperforms viterbi traceback\n\
\t-posterior\t\tCalculates posterior probabilities\n\
\t\t\tIf no output options are supplied, this will return the posterior scores\n\
\t\t\tfor all of the states.\n\n\
\t\t-threshold <score>: Return only the States with a GFF_DESC, if they are\n\
\t\t\tgreater than or equal to the threshold amount.\n\n\
\t-nbest <number of paths> \t\tperforms n-best viterbi algorithm\n\
\n\
Stochastic Decoding:\n\
\t-stochastic <Type of stochastic algorithm to use> -repetitions <number of tracebacks to sample>\n\
\t\tTypes:\n\
\t\tforward\t\t\tperforms stochastic traceback using forward algorithm\n\
\t\tviterbi\t\t\tperforms stochastic traceback using modified-viterbi algorithm\n\
\t\tposterior\t\t\tperforms stochastic traceback using posterior algorithm\n\
\n\
Output options:\n\
\t-gff\t\t\tprints path in GFF format\n\
\t-path\t\t\tprints state path according to state number\n\
\t-label\t\t\tprints state path as labels\n\
\t-hits\t\tprints hit table from stochastic sampling for each position and state\n\
\n\
Written by Paul Lott at University of California, Davis\n\
Please direct any questions, suggestions or bugs reports to Paul Lott at plott@ucdavis.edu\n\
\n\
";


#endif
