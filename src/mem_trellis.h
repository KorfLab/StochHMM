//
//  mem_trellis.h
//  StochHMM
//
//  Created by Paul Lott on 11/2/12.
//  Copyright (c) 2012 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#ifndef __StochHMM__mem_trellis__
#define __StochHMM__mem_trellis__

#include <iostream>
#include <stdint.h>
#include <vector>


namespace StochHMM {
    
    struct tb_score{
        float score;
        int16_t tb_ptr;
    };
    
    
    typedef std::vector<std::vector<uint16_t> > two_int_table;
    typedef std::vector<std::vector<float> > two_float_table;
    typedef std::vector<std::vector<std::vector<uint16_t> > > three_int_table;
    typedef std::vector<std::vector<std::vector<float> > > three_float_table;
    
    
    class mem_trellis{
    public:
        mem_trellis();
        ~mem_trellis();
                
        //Traceback Tables
        two_int_table*      traceback;          //Simple traceback table
        three_float_table*  stoch_traceback;    //Stochastic traceback table
        three_int_table*    nth_traceback;      //Nth-Viterbi traceback table
        
        //Score Tables
        two_float_table*     viterbi_score;      //Storing viterbi scores
        two_float_table*    forward_score;      //Storing Forward scores
        two_float_table*    backward_score;     //Storing Backward scores
        
        //Ending Cells
        double   ending_viterbi_score;
        uint16_t ending_viterbi_tb;
        double   ending_forward_score;
        
        std::vector<tb_score>*    ending_stoch_tb;
        std::vector<tb_score>*    ending_nth_viterbi;
    };
}


#endif /* defined(__StochHMM__mem_trellis__) */
