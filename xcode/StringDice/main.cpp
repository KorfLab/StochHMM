//
//  main.cpp
//  StringDice
//
//  Created by Paul Lott on 11/9/15.
//  Copyright © 2015 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include <iostream>
#include "StochHMMlib.h"

int main(int argc, const char * argv[]) {
	
	//Import HMM
	std::string model_file = "Dice.hmm";
	StochHMM::model hmm;
	hmm.import(model_file);
	hmm.random_walk();
	
	//Get Track from Model
	StochHMM::track* tr = hmm.getTrack(0);
	
	//Convert String to sequence
	std::string sq("123456123456");
	StochHMM::sequence* seq = new StochHMM::sequence(sq, tr);
	
	//Add sequence to sequences
	StochHMM::sequences seqs;
	seqs.addSeq(seq, tr);
	
	//Process Model
	StochHMM::trellis trell(&hmm,&seqs);
	trell.viterbi();
	
	StochHMM::traceback_path* path = new(std::nothrow) StochHMM::traceback_path(trell.get_model());
	trell.traceback(*path);
	std::cout << "Score: " << path->getScore() << std::endl;
	path->print_label();
	

    return 0;
}
