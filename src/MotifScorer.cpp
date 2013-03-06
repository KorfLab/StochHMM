//
//  MotifScorer.cpp
//  StochHMM
//
//  Created by Paul Lott on 2/14/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include <iostream>
#include <fstream>
#include "emm.h"
#include "sequence.h"
#include "sequences.h"
#include "seqTracks.h"
#include "pwm.h"

#include "MotifScorer.h"

using namespace StochHMM;

int main(int argc, const char * argv[])
{
	PWM mat;
	//mat.import("RSS23standard.txt");
	mat.import("Heptamer.txt");
	//std::cout << mat.stringify() << std::endl;
	
	sequences seqs;
	seqs.getFastas("chr22.fa",mat.getTrack());
	
	mat.score(&seqs);
	
//	seqs.addSeq(sq, mat.getTrack());
//	
//	mat.score(&seqs);
	
	//mat.score(sq);
	//std::cout << sq->getHeader() << std::endl;
	
}

