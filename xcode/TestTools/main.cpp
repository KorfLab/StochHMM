//
//  main.cpp
//  TestTools
//
//  Created by Paul Lott on 2/20/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include <iostream>
#include "StochHMMlib.h"


int main(int argc, const char * argv[])
{
	const char * characters[4] = {"A","C","G","T"};

	StochHMM::track tr;
	tr.addAlphabetChar(4,characters);
	
	StochHMM::sequence seq;
	std::string sq("ACGTACGTACGTACGTACGT");
	seq.setSeq(sq, &tr);
	
	std::cout << seq.undigitize() << std::endl;
	seq.shuffle();
	std::cout << seq.undigitize() << std::endl;
	
	
	std::vector<double> freq(4,0.25);
	
	seq = StochHMM::random_sequence(freq,100,&tr);
	
	std::cout << seq.undigitize() << std::endl;
	
	
	return 0;
}

