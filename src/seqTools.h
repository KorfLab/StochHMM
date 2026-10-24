//
//  seqTools.h
//  StochHMM
//
//  Created by Paul Lott on 5/18/12.
//  Copyright (c) 2012 University of California, Davis. All rights reserved.
//

#ifndef StochHMM_seqTools_h
#define StochHMM_seqTools_h
#include <iostream>
#include <string>
#include <algorithm>
#include "sequence.h"
#include "hmm.h"
namespace StochHMM{
    sequence shuffle(sequence* seq);
    
    sequence random_sequence(std::vector<double> const& freq, size_t , track*);
    sequence random_sequence(emm*);
    
    sequence reverseComplement();
    sequence translate();
    
    
    void motifScoring();
    
    void markov_length_distribution(model*);
    
    void markov_generate_sequence(model*);
    
//Put in PWM scoring with options for set threshold or determine threshold through shuffling and calculation of FDR(give an FDR threshold)
    
    
}



#endif
