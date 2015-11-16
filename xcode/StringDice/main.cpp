//
//  main.cpp
//  StringDice
//
//  Created by Paul Lott on 11/16/15.
//  Copyright © 2015 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

//
//  main.cpp
//  StringDice
//
//  Copyright © 2015 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include <iostream>
#include "StochHMMlib.h"

int main(int argc, const char * argv[]) {
	
	//Import HMM
	std::string model_file = "Dice_multi.hmm";
	StochHMM::model hmm;
	hmm.import(model_file);
	
	hmm.print();
	
	//Get Track from Model
	StochHMM::track* dice_tr = hmm.getTrack(0);
	StochHMM::track* color_tr = hmm.getTrack(1);
	
	
	//Convert String to sequence
	std::string sq("123456664");
	std::string col("RGBBBRRRB");
	StochHMM::sequence* seq = new StochHMM::sequence(sq, dice_tr);
	StochHMM::sequence* col_seq = new StochHMM::sequence(col, color_tr);
	
	//Add sequence to sequences
	StochHMM::sequences seqs;
	seqs.addSeq(seq, dice_tr);
	seqs.addSeq(col_seq, color_tr);
	
	seqs.print();
	
	//Process Model
	StochHMM::trellis trell(&hmm,&seqs);
	trell.viterbi();
	
	StochHMM::traceback_path* path = new(std::nothrow) StochHMM::traceback_path(trell.get_model());
	trell.traceback(*path);
	std::cout << "Score: " << path->getScore() << std::endl;
	path->print_label();
	
	
	return 0;
}

