//
//  main.cpp
//  DupHMM 1.2
//
//  Created by Ravi Dandekar on 8/20/13.
//  Copyright (c) 2013 Ravi Dandekar. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdint.h>
#include <algorithm>
#include <bitset>
#include <time.h>
#include <bitset>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string>
#include <vector>
#include <stdlib.h>
#include <iomanip>
#include <stdio.h>
#include <cstring>
#include <map>
#include <utility>

#include "StochHMMlib.h"
#include "text.h"
#include "track.h"
#include "stochMath.h"
#include "PDF.h"

static std::vector<std::vector<double>> emission_table;

static void initialize_emission_vectors(std::string& filename){
    
    std::cout << "reading " << filename << "\n";
    
    std::ifstream emissions_fh;
    emissions_fh.open(filename.c_str());
    if (!emissions_fh.is_open()){
        std::cerr << "Could not open Emission Probability File\n";
    }
    
    /* create the emission table to be the size of the emission table */
    StochHMM::stringList header_row;
    std::string header;
    getline(emissions_fh,header,'\n');
    header_row.splitString(header, "\t");
    emission_table.resize(header_row.size(), std::vector<double>(100,-INFINITY));
    std::cout << "size is " << header_row.size() << "\n";
    
    std::string line;
    size_t n(0);
    while (getline(emissions_fh, line, '\n')){
        StochHMM::stringList row;
        row.splitString(line, "\t");
        
        for (size_t model = 0; model < row.size(); model++){
            double value;
            if (!StochHMM::stringToDouble(row[model], value)) {
                std::cerr << "Couldn't convert value to double" << std::endl;
            }
            emission_table[model][n] = log(value);
        }
    }
    
    std::cout << "loaded " << filename << "\n";
}

double emission (const double val, const std::vector<double>* param) {
    int model = (*param)[0];
    return emission_table[model][val];
}



/*=============*/
/*  MAIN CODE  */
/*=============*/

// Usage Statement
std::string usage = "usage: DupHMM_1.2 <Seq File> <emission prob file> <hmm file>\n";

int main (int argc, const char * argv[]) {
	
    /* command line */
    if (argc != 4) {
        std::cout << usage << std::endl;
        exit(2);
    }
    std::string seq_file  =   argv[1];
    std::string emissions =   argv[2];
    std::string hmm_file  =   argv[3];
    
    /* init */
    initialize_emission_vectors(emissions);
    
    StochHMM::StateFuncs hmm_functions;
    hmm_functions.assignPDFFunction("P05x", *emission);
    hmm_functions.assignPDFFunction("P1x",  *emission);
    hmm_functions.assignPDFFunction("P2x",  *emission);
    hmm_functions.assignPDFFunction("P3x",  *emission);
    hmm_functions.assignPDFFunction("P4x",  *emission);
    hmm_functions.assignPDFFunction("P5x",  *emission);
    
    StochHMM::model hmm;
    std::cout << "attempting to import " << hmm_file << "\n";
    hmm.import(hmm_file, &hmm_functions);
    hmm.print();
    
	//    StochHMM::seqTracks jobs;
    
    /*
	 
	 jobs.loadSeqs(hmm, seq_file);
	 StochHMM::seqJob *job=jobs.getJob();
	 
	 while (job != NULL) {
	 
	 std::cout << job << std::endl;
	 
	 job = jobs.getJob();
	 }
     */
	
    return 0;
}
