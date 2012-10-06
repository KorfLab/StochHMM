//
//  main.cpp
//  TestNewStochHMM
//
//  Created by Paul Lott on 12/16/11.
//  Copyright (c) 2011 University of California, Davis. All rights reserved.
//

#include <iostream>
#include <string>
#include "hmm.h"
#include "seqTracks.h"

using namespace StochHMM;


int main (int argc, const char * argv[])
{
    //std::string file("newFormat_woTemplate.hmm");
    std::string file("E2.hmm");
    
    model ModelParsing(file,NULL);
    //model ModelParsing;
    ModelParsing.print();
    ModelParsing.writeGraphViz("Test.dot",false);
    
    seqTracks trks;
    std::string seqFile("Test.fa");
    trks.loadFasta(ModelParsing, seqFile);
    
    
    return 0;
}
