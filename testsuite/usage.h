//
//  usage.h
//  StochHMM
//
//  Created by Paul Lott on 3/5/12.
//  Copyright (c) 2012 University of California, Davis. All rights reserved.
//

#ifndef StochHMM_usage_h
#define StochHMM_usage_h

const char usage[]  = "Usage: HMM_Counter  -seq <Fasta File> [Options]\n\
Options:\n\
\t-help:  Print the usage statement\n\
\t-order: Define what order of dependence\n\
\t-direction (forward|reverse|both): Direction to read sequences. Forward or Reverse complement or both.\n\
\t-periodic  <INT>:  Defines the period to use.  Coding sequence =3\n\
\t-pwm:  Read aligned sequence and create a PWM\n\
\t-upstream <INT> :  Amount of upstream sequence in each sequence\n\
\t-downstream <INT> : Amount of downstream sequence in each sequence\n\
\t-pseudocount <INT> : Amount of pseudocounts to add to table\n\
\t-withWords : Prints out tables with the Words corresponding to row and column alphabet\n\
\t-masked <INT> :  the supplied value identifies how many states are identified in your sequence.\n\n\
Created by Paul Lott,  UCDavis Genome Center\n\
Contact: lottpaul@gmail.\n";


#endif
