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
#include <algorithm>
#include <stdint.h>
#include "mem_trellis.h"
#include "hmm.h"
#include "state.h"
#include "sequences.h"
#include <time.h>
#include <bitset>

namespace StochHMM {
    
    typedef std::vector<double> zero_order;
    typedef std::vector<zero_order> first_order;
    typedef std::vector<first_order> second_order;
    typedef std::vector<second_order> third_order;
    typedef std::vector<third_order> fourth_order;
    typedef std::vector<fourth_order> fifth_order;
    
    void viterbi(mem_trellis* trell, model* hmm, sequences* seqs);
    void viterbi_two(mem_trellis* trell, model* hmm, sequences* seqs);
    void viterbi_three(mem_trellis* trell, model* hmm, sequences* seqs);
	void viterbi_four(mem_trellis* trell, model* hmm, sequences* seqs);

    
    size_t* indexToindices(size_t index, size_t length);
    
    std::vector<std::bitset<200> >* process_to_trans(model *hmm);
    std::vector<std::bitset<200> >* process_from_trans(model *hmm);
    void process_initial(model* hmm, std::bitset<200>& ending);
	void permute(size_t word_length, std::vector<std::vector<size_t> >& words);
    //void expand (size_t value, std::vector<size_t>& ambiguous, std::vector<std::vector<size_t> >* templ, std::vector<std::vector<size_t> >* expanded);
	void expand_amb(size_t value, std::vector<size_t>& ambiguous, std::vector<size_t>& word, std::vector<std::vector<size_t> >& ret_val);
	void _expand_amb(size_t value, size_t position, std::vector<size_t>& ambiguous,std::vector<size_t>& word, std::vector<std::vector<size_t> >& ret_val);
	
	
	std::vector<std::vector<bool> >* process_to_trans_new(model *hmm);
    std::vector<std::vector<bool> >* process_from_trans_new(model *hmm);
    void process_initial_new(model* hmm, std::vector<bool>& ending);
	
	
	class emissions{
    public:
        emissions(emm* emiss,std::vector<short>* sq);
        void parse_emission(emm* emiss);
		void calc_ambig();
        double get_emission(size_t position);
		void add_word(std::vector<size_t>& ambig, std::vector<std::vector<size_t> >& words);
        
        size_t order;
		std::vector<short>* seq;
        
        zero_order*     zero;
        first_order*    first;
        second_order*   second;
        third_order*    third;
        fourth_order*   fourth;
        fifth_order*    fifth;
        
    };
        
    
}

#endif /* defined(__StochHMM__mem_viterbi__) */
