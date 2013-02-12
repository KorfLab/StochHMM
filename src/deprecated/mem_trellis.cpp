//
//  mem_trellis.cpp
//  StochHMM
//
//  Created by Paul Lott on 11/2/12.
//  Copyright (c) 2012 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include "mem_trellis.h"

namespace StochHMM {
    
    mem_trellis::mem_trellis(){
        
    }
    
    mem_trellis::~mem_trellis(){
        delete traceback;
        delete stoch_traceback;
        delete nth_traceback;
        delete viterbi_score;
        delete forward_score;
        delete backward_score;
        delete ending_stoch_tb;
        delete ending_nth_viterbi;
    }
    
}