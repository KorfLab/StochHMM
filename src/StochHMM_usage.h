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
StochHMM - Stochastic Hidden Markov Model Framework (version 0.25) Date: February 14, 2013\n\
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
Posterior:\n\
\t-posterior\t\tCalculates posterior probabilities\n\
\n\
Output options:\n\
\t-gff\t\t\tprints path in GFF format\n\
\t-path\t\t\tprints state path according to state number\n\
\t-label\t\t\tprints state path as labels\n\
\n\
Written by Paul Lott at University of California, Davis\n\
Please direct any questions, suggestions or bugs reports to Paul Lott at plott@ucdavis.edu\n\
\n\
";


#endif
