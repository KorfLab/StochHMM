//
//  TM_usage.h
//  StochHMMme_rc
//
//  Created by Paul C Lott on 7/1/11.
//  Copyright 2011 University of California, Davis. All rights reserved.
//

#pragma mark Usage
const char usage[]  = "\n\
StochHMMmv- Stochastic-based HMM multivariable package (version 0.1β)\n\
usage: StochHMMmv [options]\n\n\
example: StochHMMmv -model Arabidopsis.hmm -seq Arab412.fasta -stoch forward -rep 10000 -rpt 100 -gff -label -path -posterior\n\n\
Command Line options:\n\
\t-help or -h		print usage statement\n\
\n\
Files: reqires a sequence file and a model file\n\
\t-model file\t\timport model file\n\
\t-track file\t\t\timport tracks file in fasta format, comma separated list,tab-delimited list,\n\
\n\
Non-stochastic Decoding:  Different algorithms available for decoding\n\
\t-viterbi\t\t\tperforms viterbi traceback\n\
\t-nbest #\t\t\tperforms nbest traceback returning # of paths\n\
\n\
Stochastic Decoding:\n\
\t-stoch type\t\tperforms stochastic traceback\n\
\t\tTypes:\n\
\t\t\tforward\tStochastic traceback based on forward probabilities\n\
\t\t\tviterbi\tStochastic traceback based on viterbi probabilities\n\
\n\
\t-rep #\t\t\t# of stochastic traceback repetitions to perform\n\
\t-rpt #\t\t\t# top of stochastic tracebacks to report \n\
\n\
Output options:\n\
\t-gff\t\t\t\tprints path in GFF format\n\
\t-path\t\t\tprints state path according to state number\n\
\t-label\t\t\tprints state path as labels\n\
\t-posterior\t\tprints posterior probability table for labels and states\n\
\n\
Written by Paul Lott at University of California, Davis\n\
Please direct any questions, suggestions or bugs reports to Paul Lott at plott@ucdavis.edu\n\
\n\
";

