//
//  main.cpp
//  ModelParsing
//
//  Created by Paul Lott on 10/14/11.
//  Copyright 2011 University of California, Davis. All rights reserved.
//

#include <iostream>
#include <string>
#include "hmm.h"
#include "seqTracks.h"
using namespace StochHMM;

//typedef double  (*transitionFunc) (const std::string*, const size_t, const std::string*, const size_t);
//typedef double  (*emissionFunc) (const std::string*, const size_t);


double transFunc_Test(const std::string* seq, const size_t pos, const std::string* trb_seq, const size_t trb){
    return -2.0;
}


double emissFunc_PWM(const std::string* seq , const size_t pos){
    return -1.5;
}

StateFuncs stFuncs;



int main (int argc, const char * argv[])
{
    std::string pwm="PWM";
    stFuncs.assignEmmissionFunction(pwm, *emissFunc_PWM);
    
    
    
    std::string files[]={"TestModel1.hmm",
        "TestModel2.hmm","TestModel3.hmm",
        "TestModel4.hmm","TestModel5.hmm",
        "TestModel6.hmm","TestModel7.hmm",
        "TestModel8.hmm","TestModel9.hmm",
        "TestModel10.hmm","TestModel11.hmm",
        "TestModel12.hmm","TestModel13.hmm",
        "TestModel14.hmm","TestModel15.hmm",
        "TestModel16.hmm"
    };
    
    
    for(int i=0;i<10;i++){
        
        std::cout << files[i] <<std::endl;
        
        std::string file = "TestModels/ValidModels/"  + files[i];
        model test;
        if (test.import(file,&stFuncs)){
            //test.print();
        }
        else{
            std::cerr << "Error reading in model "<< file << std::endl;
        }
    }
    
    for(int i=0;i<16;i++){
        
        std::cout << files[i] <<std::endl;
        
        std::string file = "TestModels/InvalidModels/"  + files[i];
        model test;
        if (test.import(file,&stFuncs)){
            test.print();
        }
        else{
            std::cerr << "Error reading in model "<< file << std::endl;
        }
    }
    
    return 0;
}

