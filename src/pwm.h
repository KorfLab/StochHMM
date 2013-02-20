//
//  pwm.h
//Copyright (c) 2007-2012 Paul C Lott 
//University of California, Davis
//Genome and Biomedical Sciences Facility
//UC Davis Genome Center
//Ian Korf Lab
//Website: www.korflab.ucdavis.edu
//Email: lottpaul@gmail.com
//
//Permission is hereby granted, free of charge, to any person obtaining a copy of
//this software and associated documentation files (the "Software"), to deal in
//the Software without restriction, including without limitation the rights to
//use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
//the Software, and to permit persons to whom the Software is furnished to do so,
//subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
//FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
//COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
//IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
//CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef StochHMM_pwm_cpp
#define StochHMM_pwm_cpp    
    
#include <vector>
#include <math.h>
#include <algorithm>
#include <vector>
#include <list>
#include "track.h"
#include "sequences.h"
#include "stochTypes.h"
#include "emm.h"
namespace StochHMM{
	
	
	
	
	class matrixPosition;
	
	//Position Weight Matrix Class
	//Before using we need to calculate the pValues for the matrix
	class PWM{
	public:
		PWM();
		void import(std::string& file);
		void import(const char* file);
		bool parse(const std::string& matrix);
		
		void score(sequences* seqs);
		void scoreSimple(sequences* seqs);
		void scoreUndefSpacer(sequences* seqs);
		void scoreVariableSpacer(sequences* seqs);
		
		void score(sequence* seq);
		void scoreSimple(sequence* seq);
		void scoreUndefSpacer(sequence* seq);
		void scoreVariableSpacer(sequence* seq);
		
		inline void setBackground(emm* bg){bgWeight = bg;}
		inline void setTrack(track* tr){trk = tr;}
		inline void setSimpleThreshold(float thresh){simpleThreshold = thresh;};
		
		inline void setCurrentThreshold(float* thresh){currentThreshold = thresh;}
		inline track* getTrack(){return trk;}

		std::string stringify();
		void print();
	private:
		
		bool _parseTrack(std::string& txt);
		bool _parseAmbiguous(std::string& txt);
		bool _parsePositions(std::string& txt);
		bool _parseThreshold(std::string& txt);
		bool _parseBackground(std::string& txt);
		bool _parseSpacer(std::string& txt);
		bool _splitPositions(std::string& txt ,stringList& sts);
		bool _getOrderedPositionNames(stringList& states, stringList& names);
		void _finalizeTransitions();
		float calculateBack(sequences *seqs, size_t position, float sum);
		float calculateBack(sequence *seq, size_t position, float sum);
		
		track* trk;
		valueType valType;
		
		bool simple;
		bool variableSpacer;
		bool undefinedSpacer;
		float simpleThreshold;
		float* currentThreshold;
		
		std::vector<matrixPosition*> weightMatrix;
		std::map<std::string,size_t> positionNames;
		
		std::vector<matrixPosition*> frontWeightMatrix;
		std::vector<matrixPosition*> backWeightMatrix;
		std::list<float>* backScores;
		std::bitset<1024>* backScored;
		
		std::vector<size_t> undefSpacerSizes;
		std::vector<matrixPosition*> variableSpacerMatrix;
		std::bitset<1024>* variableTransition;
		std::string frontWeightName;
		std::string backWeightName;
		size_t min_spacer;
		size_t max_spacer;
		
		emm* bgWeight; //Background weight
	};
	
	
	/*! \class matrixPosition
	 Stores weight information for a position in the position weight matrix
	 */
	class matrixPosition{
	public:
		matrixPosition();
		~matrixPosition();
		bool parse(std::string& txt, track* trk, stringList& names);
		float getEmissionValue(sequences*, size_t);
		float getEmissionValue(sequence*, size_t);
		inline emm* getEmission(){return positionMatrix;};
		inline void addTransition(emm* trans){transitions.push_back(trans);}
		inline std::vector<std::string>& getTransitionNames(){return transition_names;}
		inline void setThreshold(float thresh){threshold = thresh;}
		inline float* getThresholdPtr(){return &threshold;}
		inline float getThreshold(){return threshold;}
		inline size_t transitionsSize(){return transitions.size();}
		std::string stringify();
		
	private:
		emm* positionMatrix;	//Weight information for position
		bool thresholdSet;		//True if a threshold is set for this position(global threshold)
		float threshold;		//Global threshold for position
		std::vector<emm*> transitions;	//The next weights after this position (Transitions)
		std::vector<std::string> transition_names; //Next weights after this position
		std::string name;		//Name assigned to this position
	};

	

    
}
#endif
