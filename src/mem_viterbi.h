//
//  mem_viterbi.h
//  StochHMM
//
//  Created by Paul Lott on 11/2/12.
//  Copyright (c) 2012 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#ifndef __StochHMM__mem_viterbi__
#define __StochHMM__mem_viterbi__

#include <iostream>
#include <stdint.h>
#include "mem_trellis.h"
#include "hmm.h"
#include "state.h"
#include "sequences.h"
#include <time.h>
#include <bitset>

namespace StochHMM {
    
    void viterbi(mem_trellis* trell, model* hmm, sequences* seqs);
    
    void viterbi_two(mem_trellis* trell, model* hmm, sequences* seqs);
    
    std::vector<std::bitset<400> >* process_to_trans(model *hmm);
    std::vector<std::bitset<400> >* process_from_trans(model *hmm);
    void process_initial(model* hmm, std::bitset<400>& ending);
    
}

#endif /* defined(__StochHMM__mem_viterbi__) */
