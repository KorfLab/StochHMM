//
//  main.cpp
//  TestTrellis
//
//  Created by Paul Lott on 5/11/12.
//  Copyright (c) 2012 University of California, Davis. All rights reserved.
//

#include <iostream>
#include <string>
#include "hmm.h"
#include "sequence.h"
#include "trellis.h"
#include <time.h>
using namespace StochHMM;


int main(int argc, const char * argv[])
{
    std::string file = "3_16Eddy.hmm";
    
    std::cout << "Importing Dice Model" << std::endl;
    clock_t start = clock();
    model ModelParsing;
    ModelParsing.import(file,NULL);
    ModelParsing.print();
    clock_t stop = clock();
    
    std::cout << "Import Model Time: " << (double)(stop-start)/(double)CLOCKS_PER_SEC << std::endl;
    
    tracks* trk = ModelParsing.getTracks();
    
    sequences sqs(trk->size());
    
    char *str = "GGCA";
    std::cout << "Import Track/Sequence Information" << std::endl;
    start = clock();

    track* tr=trk->getTrack("TRACK1");
    tr->setAmbiguous();
    sequence* sq = new(std::nothrow) sequence(str,trk->getTrack("TRACK1"));
    //sq->print();
    sqs.addSeq(sq);
   
    
    stop = clock();
    
    std::cout << "Import Track/Sequence Time: " << (double)(stop-start)/(double)CLOCKS_PER_SEC << std::endl;
    
    
    std::cout << "Creating Trellis " << std::endl;
    start = clock();
    trellis trell(&ModelParsing,&sqs,SIMPLE);
    stop = clock();
    
    std::cout << "Creating Trellis Time: " << (double)(stop-start)/(double)CLOCKS_PER_SEC << std::endl;
    
    
    std::cout << "Calculating Viterbi " << std::endl;
    start = clock();
    //trell.nthViterbi(5);
    trell.posterior();
    stop = clock();
    std::cout << "Calculating Viterbi Time: " << (double)(stop-start)/(double)CLOCKS_PER_SEC << std::endl;
    
    
    //trell.print();
    std::cout << "Get Traceback " << std::endl;
    start = clock();
    traceback_path& path = trell.posteriorTraceback();
    stop = clock();
    std::cout << "Get Traceback Time: " << (double)(stop-start)/(double)CLOCKS_PER_SEC << std::endl;
    
    
    std::string pth;
    
    std::cout << "Get Path from Traceback " << std::endl;
    start = clock();
    path.label(pth);
    stop = clock();
    std::cout << "Get Path from Traceback: " << (double)(stop-start)/(double)CLOCKS_PER_SEC << std::endl;
    std::cout << "Score: " << exp(path.getScore())   << std::endl;
    std::cout << pth << std::endl;
    
    std::map<std::string,int> paths;
    std::map<std::string,int>::iterator iter;
    for(size_t i=0;i<100000;++i){
        path = trell.stochasticPosterior();
        std::string tmp;
        path.label(tmp);
        paths[tmp]++;
    }
    
    for(iter=paths.begin();iter!=paths.end();++iter){
        std::cout << "Path:\t" << (*iter).first  << "\t" << (*iter).second << std::endl;
    }
    
    
    trell.print();
    
    return 0;
    
}

