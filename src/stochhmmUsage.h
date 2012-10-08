//
//  StochHMM_usage.h
//  StochHMMme_rc
//
//  Created by Paul C Lott on 7/1/11.
//  Copyright 2011 University of California, Davis. All rights reserved.
//

#pragma mark Usage
const char usage[]  = "\n\
StochHMM- Stochastic Hidden Markov Model Framework (version 0.2)\n\
usage: StochHMM -model <model_file>  -seq <seq_file> [options]\n\
\n\
example: StochHMM -model Arabidopsis.hmm -seq Arab412.fasta -viterbi -gff -label -path \n\
\n\
Command Line options:\n\
\t-help or -h		print usage statement\n\
\n\
Files: reqires a sequence file and a model file\n\
\t-model file\t\timport model file\n\
\t-seq file\t\t\timport sequence file in fasta format\n\
\n\
Non-stochastic Decoding:  Different algorithms available for decoding\n\
\t-viterbi\t\t\tperforms viterbi traceback\n\
\n\
Output options:\n\
\t-gff\t\t\t\tprints path in GFF format\n\
\t-path\t\t\tprints state path according to state number\n\
\t-label\t\t\tprints state path as labels\n\
\t-trellis\t\t\tprint the trellis values\n\
\n\
Written by Paul Lott at University of California, Davis\n\
Please direct any questions, suggestions or bugs reports to Paul Lott at plott@ucdavis.edu\n\
\n\
";

