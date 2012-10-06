//
//  seqTools.cpp
//  StochHMM
//
//  Created by Paul Lott on 5/18/12.
//  Copyright (c) 2012 University of California, Davis. All rights reserved.
//

#include "seqTools.h"

namespace StochHMM{
    
    
    sequence shuffle(sequence *seq){
        sequence shuffled_seq(ALPHA_NUM);
        std::string strSeq=seq->getUndigitized();
        for(size_t i=0;i<strSeq.size();++i){
            size_t k = rand() % (strSeq.size());
            if (k!=i){
                iter_swap(strSeq.begin()+i, strSeq.begin()+k);
            }
        }
        
        shuffled_seq.setSeq(strSeq, seq->getTrack());
        
        return shuffled_seq;
    }
    
    sequence random_sequence(std::vector<double>& freq, size_t length, track* tr){
        sequence random_seq(ALPHA_NUM);
        
        if (tr==NULL){
            std::cerr << "Track is not defined" << std::endl;
            return random_seq;
        }
        
        size_t alphaSize=tr->getAlphaSize();
        size_t freqSize=freq.size();
        if (alphaSize!=freqSize){
            std::cerr << "Frequency distribution size and Alphabet size must be the same." << std::endl;
            return random_seq;
        }
        
        //Create CDF of frequency distribution
        std::vector<std::pair<double,std::string> > cdf;
        double sum = 0.0;
        for(size_t i=0;i<freqSize;++i){
            sum+=freq[i];
            std::pair<double,std::string> val (sum, tr->getAlpha(i));
            cdf.push_back(val);
        }
        
        //Generate random sequence
        std::string random_string;
        for(size_t j=0;j<length;++j){
            double val = ((double)rand()/((double)(RAND_MAX)+(double)(1))); //Generate random 
            for (size_t m=0;m<freqSize;++m){ //Check to see which value is 
                if (cdf[m].first>=val){
                    random_string+=cdf[m].second;
                    break;
                }
            }
        }
        
        random_seq.setSeq(random_string, tr);
        return random_seq;
    }
}