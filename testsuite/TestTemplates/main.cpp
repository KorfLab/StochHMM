//
//  main.cpp
//  TestTemplates
//
//  Created by Keith Dunaway on 10/21/11.
//  Copyright 2011 UC Davis. All rights reserved.
//

#include <vector>
#include <string>
#include <iostream>
#include <map>
#include "text.h"
#include "modelTemplate.h"
//using namespace StochHMM;

int main (int argc, const char * argv[])
{
    std::map<std::string,std::string> ReplacedUserVals;
    ReplacedUserVals["TABLE"] = "0,0,0,0,\n0,0,0,0,\n0,0,0,0,\n";
    ReplacedUserVals["SS_1"] = "1,1,1,1,\n1,1,1,1,\n1,1,1,1,\n";
    ReplacedUserVals["VAL"] = "2,2,2,2,\n2,2,2,2,\n2,2,2,2,\n";
    ReplacedUserVals["SS_2"] = "3,3,3,3,\n3,3,3,3,\n3,3,3,3,\n";
    ReplacedUserVals["VAL2"] = "4,4,4,4,\n4,4,4,4,\n4,4,4,4,\n";
    std::string filestring= "/Users/keithdunaway/TestTemplate.txt";
    std::string modelstring = StochHMM::slurpFile(filestring);
    //std::cout << modelstring << "\n";
    StochHMM::templates Tester(modelstring);
    std::string printmodelname = "START_SITE";
    std::string IDnumber = "001";
    
    std::string printablestring = Tester.getTemplate(printmodelname, IDnumber,ReplacedUserVals);
    
    std::cout << printablestring << "\n";
    return 0;
}

