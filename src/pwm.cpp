//
//  pwm.cpp
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

#include "pwm.h"
namespace StochHMM{

	//! Constructor for PWM class
	PWM::PWM(){
		valType=PROBABILITY;
		simple = true;
		simpleThreshold = -INFINITY;
		currentThreshold = &simpleThreshold;
		bgWeight = NULL;
		min_spacer = SIZE_MAX;
		max_spacer = 0;
		return;
	}
	
	/*! Import a Position Weight Matrices file
	 Imports file and parses the file position weight matrix file
	 \param[in] std::string file
	 */
	void PWM::import(std::string& file){
		std::string modelString=slurpFile(file);
		parse(modelString);
	}
	
	/*! Import a Position Weight Matrices file
	 Imports file and parses the file position weight matrix file
	 \param[in] const char* file
	 */
	void PWM::import(const char* file){
		std::string fl(file);
		import(fl);
	}
	
	/*! Parses a Position Weight Matrices string
	 Parses string containing position weight matrix definitions
	 \param[in] std::string file
	 */
	bool PWM::parse(const std::string& matrix){
		size_t thresh = matrix.find("THRESHOLD DEFINITION");
        size_t track  = matrix.find("TRACK SYMBOL DEFINITIONS");
        size_t ambig  = matrix.find("AMBIGUOUS SYMBOL DEFINITIONS");
        size_t pwm = matrix.find("POSITION WEIGHT DEFINITIONS");
		size_t back = matrix.find("BACKGROUND DEFINITION");
		size_t space = matrix.find("SPACER DEFINITIONS");
		size_t blank;
        size_t nlChar;
		
		if (thresh != std::string::npos){
			blank=matrix.find("\n\n",thresh);
			
			size_t nlCharEq = matrix.rfind("####\n",blank);
			size_t nlCharNum= matrix.rfind("====\n",blank);
			//Check for optional dividing line
			if (nlCharEq!=std::string::npos){
				nlChar=nlCharEq+5;
			}
			else if (nlCharNum!=std::string::npos){
				nlChar=nlCharNum+5;
			}
			else{  //No divider line
				nlChar=matrix.find("\n",thresh);
				nlChar++;
			}
			
			
			std::string thr (matrix.substr(nlChar,blank-nlChar));
			
			if (!_parseThreshold(thr)){
				return false;
			}
			
		}
		
		if (track != std::string::npos){
			blank=matrix.find("\n\n",track);

            size_t nlCharEq = matrix.rfind("####\n",blank);
            size_t nlCharNum= matrix.rfind("====\n",blank);
            //Check for optional dividing line
            if (nlCharEq!=std::string::npos){
                nlChar=nlCharEq+5;
            }
            else if (nlCharNum!=std::string::npos){
                nlChar=nlCharNum+5;
            }
            else{  //No divider line
                nlChar=matrix.find("\n",track);
                nlChar++;
            }


            std::string trck (matrix.substr(nlChar,blank-nlChar));

            if (!_parseTrack(trck)){
                return false;
            }
			
        }
        else{
            std::cerr << "Required section: TRACK SYMBOL DEFINITIONS missing from the model" << std::endl;
            return false;
        }
		
		if (ambig != std::string::npos){
			blank=matrix.find("\n\n",ambig);

            size_t nlCharEq = matrix.rfind("####\n",blank);
            size_t nlCharNum= matrix.rfind("====\n",blank);
            //Check for optional dividing line
            if (nlCharEq!=std::string::npos){
                nlChar=nlCharEq+5;
            }
            else if (nlCharNum!=std::string::npos){
                nlChar=nlCharNum+5;
            }
            else{  //No divider line
                nlChar=matrix.find("\n",ambig);
                nlChar++;
            }

            std::string amb(matrix.substr(nlChar,blank-nlChar));

            if (!_parseAmbiguous(amb)){
                return false;
            }

		}
		
		if (back!= std::string::npos){
			blank=matrix.find("\n\n",back);
			
            size_t nlCharEq = matrix.rfind("####\n",blank);
            size_t nlCharNum= matrix.rfind("====\n",blank);
            //Check for optional dividing line
            if (nlCharEq!=std::string::npos){
                nlChar=nlCharEq+5;
            }
            else if (nlCharNum!=std::string::npos){
                nlChar=nlCharNum+5;
            }
            else{  //No divider line
                nlChar=matrix.find("\n",back);
                nlChar++;
            }
			
            std::string background(matrix.substr(nlChar,blank-nlChar));
			
            if (!_parseBackground(background)){
                return false;
            }
		}
		
		//Parse the positions
		if (pwm != std::string::npos){
			std::string positions = matrix.substr(pwm);
			_parsePositions(positions);
		}
		
		
		//Parse the Spacer Information
		if (space != std::string::npos){
			blank=matrix.find("\n\n",space);
			
            size_t nlCharEq = matrix.rfind("####\n",blank);
            size_t nlCharNum= matrix.rfind("====\n",blank);
            //Check for optional dividing line
            if (nlCharEq!=std::string::npos){
                nlChar=nlCharEq+5;
            }
            else if (nlCharNum!=std::string::npos){
                nlChar=nlCharNum+5;
            }
            else{  //No divider line
                nlChar=matrix.find("\n",space);
                nlChar++;
            }
			
            std::string spacer(matrix.substr(nlChar,blank-nlChar));
			
            if (!_parseSpacer(spacer)){
                return false;
            }
		}
		
		_finalizeTransitions();
		
		return true;
	}
	
	matrixPosition::matrixPosition(){
		positionMatrix	= NULL;
		thresholdSet	= false;
		threshold		= -INFINITY;
	}
	
	matrixPosition::~matrixPosition(){
		delete positionMatrix;
		positionMatrix = NULL;
	}
	
	
	//!Parses the emission for each position from a string
    //! \param txt String representation of emissions
    //! \param names stringList of all state names defined in the model
    //! \param trks Tracks defined in the model
    //! \param wts Weight defined of the model
    //! \param funcs StateFunction defined for the model
    bool matrixPosition::parse(std::string& txt, track* trk, stringList& names){
		
		size_t nameHeader = txt.find("NAME:");
		size_t transHeader = txt.find("TRANSITION:");
		size_t thresholdHeader = txt.find("THRESHOLD:");
		//size_t emmHeader = txt.find("EMISSION:");
		size_t end = txt.find("//END");
		
		if (end != std::string::npos){
			txt = txt.substr(0,end);
		}
		
		//Parse Name
		if (nameHeader != std::string::npos){
			size_t endline = txt.find("\n",nameHeader);
			std::string temp_name = txt.substr(nameHeader,endline-nameHeader);
			stringList tmp;
			tmp.splitString(temp_name, ":");
			clear_whitespace(tmp[1], "\t\n ");
			name = tmp[1];
		}
		else{
			std::cerr << "Position is missing NAME definition\n"<< txt << std::endl;
			exit(2);
		}
		
		
		//Parse Transitions
		if (transHeader != std::string::npos){
			std::string temp;
			size_t endline = txt.find("\n",transHeader);
			temp = txt.substr(transHeader, endline - transHeader);
			stringList transi;
			transi.splitString(temp,":,");
			for(size_t i = 1; i < transi.size(); i++){
				clear_whitespace(transi[i], " \t\n");
				transition_names.push_back(transi[i]);
			}
		}
		
		
		//Parse Thresholds
		if (thresholdHeader != std::string::npos){
			size_t endline = txt.find("\n",thresholdHeader);
			std::string temp = txt.substr(thresholdHeader, endline-thresholdHeader);
			stringList tmp;
			tmp.splitString(temp, ":");
			clear_whitespace(tmp[1], "\t\n ");
			threshold = atof(tmp[1].c_str());
			thresholdSet = true;
		}
		
		
        stringList lst;
        lst.splitND(txt,"EMISSION:");
        
		emm* temp_emm = new(std::nothrow) emm;
		
		if (temp_emm==NULL){
			std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
			exit(1);
		}
		
		if (!temp_emm->parse(lst[0],trk)){
			return false;
		}
		
		positionMatrix = temp_emm;
        
        return true;
    }
	
	
	std::string PWM::stringify(){
		std::string pwm;
		std::string split("#################################\n");
		pwm += "#POSITION WEIGHT MATRIX\n\n";
		
		if ( simpleThreshold!= -INFINITY){
			pwm +="<THRESHOLD DEFINITION>\n";
			pwm += split;
			pwm += double_to_string(simpleThreshold);
			pwm += "\n\n";
		}
		
		pwm +="<TRACK SYMBOL DEFINITIONS>\n";
		pwm += split;
		pwm += trk->stringify();
		pwm += "\n";
		
		
		if (bgWeight != NULL){
			pwm += "<BACKGROUND DEFINITION>\n";
			pwm += split;
			pwm += bgWeight->stringify();
		}
		
		if (undefinedSpacer){
			pwm += "<SPACER DEFINITIONS>\n";
			pwm += split;
			pwm += frontWeightName + "-SPACER-" + backWeightName + ":";
			pwm += join(undefSpacerSizes, ',');
			pwm += "\n\n";
		}
		
		
		pwm +="<POSITION WEIGHT DEFINITIONS>\n";
		pwm += split;
		
		
		for(size_t i=0; i< weightMatrix.size(); i++){
			pwm += weightMatrix[i]->stringify();
			pwm += split;
		}
		pwm += "//END";
		return pwm;
	}
	
	std::string matrixPosition::stringify(){
		std::string mat;
		mat = "NAME:\t" + name + "\n";
		
		if (thresholdSet){
			mat+= "THRESHOLD:\t" + double_to_string(threshold) + "\n";
		}
		
		if (transition_names.size() > 0){
			mat += "TRANSITIONS:" + join(transition_names, ',') + "\n";
		}
		
		mat += positionMatrix->stringify();
		return mat;
	}
	
	
	void PWM::score(sequences* seqs){
		if (simple){
			scoreSimple(seqs);
		}
		
		if (undefinedSpacer){
			scoreUndefSpacer(seqs);
		}
		
		if (variableSpacer){
			scoreVariableSpacer(seqs);
		}
		return;
	}
	
	void PWM::score(sequence* seq){
		if (simple){
			scoreSimple(seq);
		}
		
		if (undefinedSpacer){
			scoreUndefSpacer(seq);
		}
		
		if (variableSpacer){
			scoreVariableSpacer(seq);
		}
		return;
	}
	
	
	float matrixPosition::getEmissionValue(sequences* seqs, size_t pos){
		return positionMatrix->get_emission(*seqs, pos);
	}
	
	float matrixPosition::getEmissionValue(sequence* seq, size_t pos){
		return positionMatrix->get_emission(*seq, pos);
	}
	
	
	
	void PWM::scoreSimple(sequences* seqs){
		size_t seq_size = seqs->getLength();
		size_t motif_size = weightMatrix.size();
		
		if (seq_size < motif_size){
			return;
		}
		
		float score(0);
		for(size_t position = 0; position < seq_size - motif_size; position++){
			score = 0;
			for (size_t motif_pos = 0; motif_pos < motif_size; motif_pos ++){
				score += weightMatrix[motif_pos]->getEmissionValue(seqs,position+motif_pos);
				if (score <= *currentThreshold){
					break;
				}
			}
			
			if (score <= *currentThreshold){
				continue;
			}
			
			std::cout << position << "\t" << score << "\n";
		}
	}
	
	void PWM::scoreSimple(sequence* seq){
		size_t seq_size = seq->getLength();
		size_t motif_size = weightMatrix.size();
		
		if (seq_size < motif_size){
			return;
		}
		
		float score(0);
		for(size_t position = 0; position < seq_size - motif_size; position++){
			score = 0;
			for (size_t motif_pos = 0; motif_pos < motif_size; motif_pos ++){
				score += weightMatrix[motif_pos]->getEmissionValue(seq,position+motif_pos);
				if (score <= *currentThreshold){
					break;
				}
			}
			
			if (score <= *currentThreshold){
				continue;
			}
			
			std::cout << position << "\t" << score << "\n";
		}
	}
	
	void PWM::scoreUndefSpacer(sequences* seqs){
		size_t seq_size = seqs->getLength();
		size_t motif_size = weightMatrix.size();
		size_t front_size = frontWeightMatrix.size();
		size_t back_size = backWeightMatrix.size();
		
		if (seq_size < motif_size + max_spacer){
			return;
		}
		
		float front_score(0);
		float back_score(0);
		float sumScore(0);
		size_t spacerSize(0);
		for(size_t position = 0; position < seq_size - (motif_size + max_spacer); position++){
			front_score = 0;
			back_score = 0;
			sumScore = 0;
			//Calculate the Front motif score
			for (size_t front_pos = 0; front_pos < front_size; front_pos++){
				front_score += frontWeightMatrix[front_pos]->getEmissionValue(seqs,position+front_pos);
				if (front_score <= *currentThreshold){
					break;
				}
			}
			
			if (front_score <= *currentThreshold){
				continue;
			}
			
			//Calculate the back scores
			for(size_t sp = 0; sp < undefSpacerSizes.size(); sp++){
				spacerSize = undefSpacerSizes[sp];
				
				//Check to see if back_motif was previously calculated.
				//If it was then get score and assign score.
				
				//If not then calculate
				size_t back_pos(0);
				for(; back_pos < back_size; back_pos++){
					back_score += backWeightMatrix[back_pos]->getEmissionValue(seqs, position+front_size+spacerSize+back_pos);
					sumScore = front_score+back_score;
					if (front_score+back_score <= *currentThreshold){
						break;
					}
				}
				
				//TODO: Assign score to list
				//Assign score of backward
				if (back_pos == back_size){
					//std::cout << "We've made it through" << std::endl;
				}
				
				if (sumScore >= *currentThreshold){
					std::cout << position << "\t" <<spacerSize << "\t" << sumScore << "\n";
				}
			}
		}
		return;
	}
	
	
	void PWM::scoreUndefSpacer(sequence* seq){
		size_t seq_size = seq->getLength();
		size_t motif_size = weightMatrix.size();
		size_t front_size = frontWeightMatrix.size();
		size_t back_size = backWeightMatrix.size();
		
		if (seq_size < motif_size + max_spacer){
			return;
		}
		
		float front_score(0);
		float back_score(0);
		float sumScore(0);
		size_t spacerSize(0);
		for(size_t position = 0; position < seq_size - (motif_size + max_spacer); position++){
			front_score = 0;
			back_score = 0;
			sumScore = 0;
			//Calculate the Front motif score
			for (size_t front_pos = 0; front_pos < front_size; front_pos++){
				front_score += frontWeightMatrix[front_pos]->getEmissionValue(seq,position+front_pos);
				if (front_score <= *currentThreshold){
					break;
				}
			}
			
			if (front_score <= *currentThreshold){
				continue;
			}
			
			//Calculate the back scores
			for(size_t sp = 0; sp < undefSpacerSizes.size(); sp++){
				spacerSize = undefSpacerSizes[sp];
				
				//Check to see if back_motif was previously calculated.
				//If it was then get score and assign score.
				
				//If not then calculate
				size_t back_pos(0);
				for(; back_pos < back_size; back_pos++){
					back_score += backWeightMatrix[back_pos]->getEmissionValue(seq, position+front_size+spacerSize+back_pos);
					sumScore = front_score+back_score;
					if (front_score+back_score <= *currentThreshold){
						break;
					}
				}
				
				//TODO: Assign score to list
				//Assign score of backward
				if (back_pos == back_size){
					//std::cout << "We've made it through" << std::endl;
				}
				
				if (sumScore >= *currentThreshold){
					std::cout << position << "\t" <<spacerSize << "\t" << sumScore << "\n";
				}
			}
		}
		return;
	}
	
	
	void PWM::scoreVariableSpacer(sequences* seqs){
		size_t seq_size = seqs->getLength();
		size_t motif_size = weightMatrix.size();
		size_t front_size = frontWeightMatrix.size();
		//size_t back_size = backWeightMatrix.size();
		
		if (seq_size < motif_size + max_spacer){
			return;
		}
		
		float score(0);
		float sumScore(0);
		size_t spacerSize(0);
		for(size_t position = 0; position < seq_size - (motif_size + max_spacer); position++){
			score = 0;
			sumScore = 0;
			//Calculate the Front motif score
			for (size_t front_pos = 0; front_pos < front_size; front_pos++){
				score += frontWeightMatrix[front_pos]->getEmissionValue(seqs,position+front_pos);
				if (score <= *currentThreshold){
					break;
				}
			}
			
			if (score <= *currentThreshold){
				continue;
			}
			
			//Calculate the back scores
			for(size_t sp = 0; sp < variableSpacerMatrix.size(); sp++){
				
				//Check to see if back_motif was previously calculated.
				//If it was then get score and assign score.
				
				//If not then calculate
				
				score+=variableSpacerMatrix[sp]->getEmissionValue(seqs, position+front_size+sp);
				
				if ((*variableTransition)[sp]){
					sumScore = calculateBack(seqs, position+front_size+sp+1, score);
					if (sumScore >= *currentThreshold){
						std::cout << position << "\t" <<spacerSize << "\t" << sumScore << "\n";
					}
				}
				
				
			}
		}
		return;
	}
	
	void PWM::scoreVariableSpacer(sequence* seq){
		size_t seq_size = seq->getLength();
		size_t motif_size = weightMatrix.size();
		size_t front_size = frontWeightMatrix.size();
		//size_t back_size = backWeightMatrix.size();
		
		if (seq_size < motif_size + max_spacer){
			return;
		}
		
		float score(0);
		float sumScore(0);
		size_t spacerSize(0);
		for(size_t position = 0; position < seq_size - (motif_size + max_spacer); position++){
			score = 0;
			sumScore = 0;
			//Calculate the Front motif score
			for (size_t front_pos = 0; front_pos < front_size; front_pos++){
				score += frontWeightMatrix[front_pos]->getEmissionValue(seq,position+front_pos);
				if (score <= *currentThreshold){
					break;
				}
			}
			
			if (score <= *currentThreshold){
				continue;
			}
			
			//Calculate the back scores
			for(size_t sp = 0; sp < variableSpacerMatrix.size(); sp++){
				
				//Check to see if back_motif was previously calculated.
				//If it was then get score and assign score.
				
				//If not then calculate
				
				score+=variableSpacerMatrix[sp]->getEmissionValue(seq, position+front_size+sp);
				
				if ((*variableTransition)[sp]){
					sumScore = calculateBack(seq, position+front_size+sp+1, score);
					if (sumScore >= *currentThreshold){
						std::cout << position << "\t" <<spacerSize << "\t" << sumScore << "\n";
					}
				}
				
				
			}
		}
		return;
	}
	
	
	float PWM::calculateBack(sequences *seqs, size_t position, float sum){
		float score = sum;
		for (size_t pos = 0; pos < backWeightMatrix.size(); pos++){
			sum += backWeightMatrix[pos]->getEmissionValue(seqs, position+pos);
			if (sum < *currentThreshold){
				return -INFINITY;
			}
		}
		return score;
	}
	
	
	
	float PWM::calculateBack(sequence *seq, size_t position, float sum){
		float score = sum;
		for (size_t pos = 0; pos < backWeightMatrix.size(); pos++){
			sum += backWeightMatrix[pos]->getEmissionValue(seq, position+pos);
			if (sum < *currentThreshold){
				return -INFINITY;
			}
		}
		return score;
	}
	
	
	
	bool PWM::_parseSpacer(std::string& txt){
		undefinedSpacer = true;
		simple = false;
		variableSpacer = false;
	
		stringList lst;
		lst.splitString(txt,":");
		
		
		//Split header to see where to insert spacer
		stringList order;
		order.splitString(lst[0],"-");
		
		frontWeightName = order[0];
		backWeightName = order[2];
		
		//Push Frontside Weights into Matrix
		if (positionNames.count(frontWeightName) == 1){
			size_t frontIndex = positionNames[frontWeightName];
			for(size_t i=0;i <= frontIndex; i++){
				frontWeightMatrix.push_back(weightMatrix[i]);
			}
		}
		else{
			std::cerr << "Couldn't find " << frontWeightName << " in the POSITION WEIGHT DEFINITIONS" << std::endl;
			exit(2);
		}
		
		//Push backside weights into Matrix
		if (positionNames.count(backWeightName) == 1){
			size_t backIndex = positionNames[backWeightName];
			for(size_t i=backIndex;i < weightMatrix.size() ; i++){
				backWeightMatrix.push_back(weightMatrix[i]);
			}
		}
		else{
			std::cerr << "Couldn't find " << backWeightName << " in the POSITION WEIGHT DEFINITIONS" << std::endl;
			exit(2);
		}
		
		//Parse Spacer Sizes
		stringList spc;
		spc.splitString(lst[1],",");
		for(size_t i=0;i<spc.size();i++){
			clear_whitespace(spc[i], " \t\n");
			size_t val = atoi(spc[i].c_str());
			if (val > max_spacer){
				max_spacer = val;
			}
			
			if (val < min_spacer){
				min_spacer = val;
			}
			
			undefSpacerSizes.push_back(val);
		}
		
		return true;
	}
	
	
	bool PWM::_parseThreshold(std::string& txt){
		simpleThreshold = atof(txt.c_str());
		return true;
	}
	
	bool PWM::_parseBackground(std::string& txt){
		emm* temp_emm = new(std::nothrow) emm;
		
		if (temp_emm==NULL){
			std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
			exit(1);
		}
		
		if (!temp_emm->parse(txt,trk)){
			return false;
		}
		
		bgWeight = temp_emm;
		
		return true;
	}

	
	
	void PWM::_finalizeTransitions(){
		
		//Assign emission transitions
		for(size_t i=0; i < weightMatrix.size(); i++){
			std::vector<std::string>& tmp= weightMatrix[i]->getTransitionNames();
			for (size_t j = 0; j < tmp.size(); j++) {
				size_t itera = positionNames[tmp[j]];
				weightMatrix[i]->addTransition(weightMatrix[itera]->getEmission());
			}
			
			if (tmp.size() > 1){
				simple = false;
				variableSpacer = true;
				undefinedSpacer = false;
			}
		}
		
		
		//TODO: DETERMINE CORE REGIONS AND THEN ADAPT SCORING TO THOSE
		if (variableSpacer){
			
			//Determine Front Core
			//Check the Transitions and determine core regions
			size_t firstMultiTrans(SIZE_MAX);
			size_t backWeightStart(SIZE_MAX);
			for(size_t i=0 ; i < weightMatrix.size(); i++){
				
				
				if (weightMatrix[i]->transitionsSize() > 1){
					firstMultiTrans = i;
					std::vector<std::string>& tmp = weightMatrix[i]->getTransitionNames();
					for( size_t j = 0; j < tmp.size(); j++){
						size_t indx = positionNames[tmp[i]];
						if (indx > i+1){
							backWeightStart = indx;
							break;
						}
					}
					break;
				}
				else{
					frontWeightMatrix.push_back(weightMatrix[i]);
				}
			}
			
			//Determine Spacer and transitions
			variableTransition = new (std::nothrow) std::bitset<1024>;
			for(size_t i=firstMultiTrans; i< backWeightStart; i++){
				variableSpacerMatrix.push_back(weightMatrix[i]);
				if (weightMatrix[i]->transitionsSize() > 1){
					(*variableTransition)[i]=1;
				}
			}
			
			
			//Determine Back Core
			for(size_t i = backWeightStart; i < weightMatrix.size() ; i++){
				backWeightMatrix.push_back(weightMatrix[i]);
			}
		}
		return;
	}
	

	bool PWM::_parseTrack(std::string& txt){
		stringList lst;
		lst.splitString(txt, "\n");
		
		trk=new(std::nothrow) track();
			
		if (trk==NULL){
			std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
			exit(1);
		}
		
		if (trk->parse(lst[0])){
			trk->setIndex(0);
		}
		else {
			std::string info = "Couldn't parse new track line.  Please check formatting of : " + lst[0];
			std::cerr << info << std::endl;
			exit(1);
		}
			
		return true;
	}
	
	bool PWM::_parseAmbiguous(std::string& txt){
        stringList lst;
        lst.splitString(txt, "\n");
        for(size_t i=0;i<lst.size();i++){
            stringList ln;
            ln.splitString(lst[i],":\t ");
			
            if (trk!=NULL){
                if (!trk->parseAmbiguous(ln[1])){
                    std::cerr << "Couldn't parse the Ambiguous section for " << ln[0] << std::endl;
                    return false;
                }
            }
            else{
                std::string info = "Ambiguous Characters Section\nSupplied track name doesn't correspond to any previously parsed tracks.\nPlease check the formatting and names.\n Unfound Name: " + ln[0];
                
                std::cerr << info << std::endl;
                return false;
            }
            
        }
        return true;
    }
	
	bool PWM::_parsePositions(std::string& txt){
		//3. create and parse states
        
        stringList positions;
        _splitPositions(txt,positions);
		
		stringList names;
        _getOrderedPositionNames(positions, names);
		for(size_t i=0;i<names.size();i++){
			positionNames[names[i]] = i ;
		}
        
        for(size_t iter=1;iter<positions.size();iter++){
			
			matrixPosition* temp = new (std::nothrow) matrixPosition;
			temp->parse(positions[iter], trk, names);
			weightMatrix.push_back(temp);
			
			//std::cout << positions[iter] << std::endl;
			
		}
        
        return true;
    }
	
	//Split states into individual state strings
    bool PWM::_splitPositions(std::string& txt ,stringList& sts){
        
        size_t start=0;
        size_t end=0;
        
        while(start!=std::string::npos){
            end=txt.find("NAME:",start+1);
            
            std::string st = txt.substr(start,end-start);
            
            clear_whitespace(st, "#><");
            
			sts.push_back(st);
			
            start=txt.find("NAME:",end);
            
        }
        return true;
    }
	
	
	//Get all Name in order
	bool PWM::_getOrderedPositionNames(stringList& pos, stringList& names){
        
		for(size_t i=0;i<pos.size();i++){
            size_t nameHeader=pos[i].find("NAME:");
			if (nameHeader == std::string::npos){
				continue;
			}
            size_t nameLineEnding=pos[i].find_first_of("\n",nameHeader);
            std::string name = pos[i].substr(nameHeader+5,nameLineEnding-(nameHeader+5));
            clear_whitespace(name, " \t\n");
            if (names.contains(name)){
                std::cerr << "Position with name of: " << name << " is defined twice in the model\n";
                return false;
            }
            else{
                names.push_back(name);
            }
        }
        return true;
    }
	
	
	
	
				
}