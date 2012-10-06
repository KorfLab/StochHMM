//
//  main.cpp
//  TestAlignment
//
//  Created by Paul Lott on 5/18/12.
//  Copyright (c) 2012 University of California, Davis. All rights reserved.
//

#include <iostream>
#include "alignment.h"
#include "seqTools.h"

int main(int argc, const char * argv[])
{
    srand(time(NULL));
    StochHMM::track first;
    first.setAlphaType(StochHMM::ALPHA_NUM);
    first.addAlphabetChar("A");
    first.addAlphabetChar("C");
    first.addAlphabetChar("G");
    first.addAlphabetChar("T");
    
    char* targetSeq = "ACGTTCGTACGT";
    char* querySeq  = "TAGTCGTGGCTA";
    
    StochHMM::sequence target(targetSeq,&first);
    StochHMM::sequence query(querySeq,&first);
    
    StochHMM::sequence rand = StochHMM::shuffle(&target);
    std::cout << rand.getUndigitized() << std::endl;
    
    StochHMM::alignment test;
    
    double score = test.align(target,query,StochHMM::cLocal,5.0,-1.0,-5.0,-2.5);
    //double pvalue = test.calc_pvalue();
    
    
    // insert code here...
    std::cout << "Hello, World!\n";
    return 0;
}

