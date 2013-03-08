//
//  externDefinitions.cpp
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

#include "externDefinitions.h"


namespace StochHMM{
    
    //!Copy constructor for ExDefSequence
//    ExDefSequence::ExDefSequence(const ExDefSequence& rhs){
//		defs = rhs.defs;
//    }
    
//    ExDefSequence& ExDefSequence::operator=(const ExDefSequence& rhs){
//		for(size_t i=0; i < defs.size() ;++i){
//			if (defs.defined(i)){
//				defs[i] = new (std::nothrow) 
//			}
//		}
//        return *this;
//    }
    
    
    //!Print the External definitions to stdout 
    void ExDefSequence::print(){
        for(size_t i = 0;i<defs.size();i++){
			if (defs.defined(i)){
				std::cout << i  << "\t" << defs[i]->stringify() << std::endl;
			}
        }
        return;
    }
    
    //!Defined if there is a external definition defined for the position
    //! \param position Position of the sequence to check for external definition
    bool ExDefSequence::defined(size_t position){
        
        if (defs.defined(position)){
            return true;
        }
        else{
            return false;
        }
    }
    
    //!External definitions are either absolute or weighted.  Checks to see if the 
    //!external definition at position in sequence is absolute
    //! \param position Position in teh sequence
    //! \return true if the external definition is absolute
    bool ExDefSequence::isAbsolute(size_t position){
		if (defs.defined(position) && defs[position]->isAbsolute()){
			return true;
		}
	
		return false;
    }
    
    //!Get the absolute state defined for the position in the sequence
    //!\param position Position of the sequence
    //!\return integer indice to state
    size_t ExDefSequence::getAbsState(size_t position){
		
        if(defs.defined(position) && defs[position]->isAbsolute()){
            return defs[position]->getState();
        }
        else{
            std::cerr << "Calling getAbsState on weighted state" << std::endl;
            return -1;
        }
    }
    
    //!Check to see if the external definition at the position is weighted and not absolute
    //! \param position Position in the sequence
    //! \return true if the external definition is weighted
    //
    bool ExDefSequence::isWeighted(size_t position){
		if (defs.defined(position)){
			return !defs[position]->isAbsolute();
		}
        return false;
    }

    //! Get the weight for the external definition at a position in the sequence
    //! \param position  Position in the sequence
    //! \param stateIter Index of state to get weight
    //! \return double value of weight to apply to state at the position
    double ExDefSequence::getWeight(size_t position, size_t stateIter){
		if (defs.defined(position)){
            return defs[position]->getWeight(stateIter);
        }
        else{
            return 0.0;
        }
    }
    
    
    //!Create a ExDef type
    //!Default: absolute=false, weightedState = -2
    ExDef::ExDef(){
        absolute=false;
        weightedState=-2;
        //st=NULL;
    }

    //!Create a weightDef type
    weightDef::weightDef(size_t state_size):ExDef(), weights(state_size){
        absolute=false;
        weightedState=SIZE_MAX;
		//weights.assign(state_size, 0);
        //st=NULL;
    }
    
    //!Assign a weight to a particular state
    //! \param stateIter  integer Iterator to the state
    //! \param logValue  Log value of weight to apply 
    void weightDef::assignWeight(size_t stateIter, double logValue){
		if (!weights.defined(stateIter)){
			weights[stateIter] = logValue;
		}
		else{
			weights[stateIter]+=logValue;
		}
        
        return;
    }
    
    
    //!Get the string representation of the ExDef
    // \return std::string representation of the ExDef
    std::string ExDef::stringify(){
        
        std::string output;
        output+="STATE: ";
        output+= int_to_string(weightedState) + "\n";
        return output;
    }

    //! Get the string representation of the weightDef
    // \return std::string representation of the weightDef
    std::string weightDef::stringify(){
        std::string test;
        
        test+="STATES: ";
        
        for(size_t i=0;i<weights.size();i++){
            if (weights[i]!=0.0){
                test+= "\t" + int_to_string(i) + "\t" + double_to_string(weights[i]) + "\n";
            }
        }
        return test;
    }
    
    //TODO: Change return so it will return false if not able to parse from file
    
    //! Parses the ExDefSequence from a file stream
    //! \param file File stream to be used to parse the External definitions from
    //! \return true if parsing was successful
    bool ExDefSequence::parse(std::ifstream& file,stateInfo& info){
        //use getDefs to parse the lines
        //and create external definition
        
        stringList ln;
        
        std::string comment_char="#";
        std::string ws=" []";
        std::string split="\t:";
        
        while(ln.fromDef(file, ws, split)){

            if (ln[0].compare("EXDEF")==0){

                if (ln[1].compare("ABSOLUTE")==0){
                    _parseAbsDef(ln,info);
                }
                else if (ln[1].compare("WEIGHTED")==0){
                    _parseWeightDef(ln,info);
                }
                else{
                    std::cerr << ln[0] << " not a valid external definition\n";
                    continue;
                }
            }
            
            char nl_peek=file.peek();
            if (nl_peek!='['){
                break;
            }
        }
        
        //defs->print();
        
        return true;
    }

    //! Parse the absolute def from stringList give the stateInfo
    bool ExDefSequence::_parseAbsDef(stringList& ln, stateInfo& info){
        bool start=false;
        bool stop=false;
        bool state=false;
        
        size_t startPosition(SIZE_MAX);
        size_t stopPosition(SIZE_MAX);
        std::vector<size_t> path;
        
        for (size_t i=2;i<ln.size();i++){
            std::string& tag=ln[i];
            if (ln.size()<=i+1){
                std::cerr << "Missing additional information for Absolute Definition: TAG: " << ln[i] << std::endl;
                return false;
            }
            
            if (tag.compare("START")==0){
                
                int tempInt;
                if (!stringToInt(ln[i+1], tempInt)){
                    std::cerr << "Value in External Definition START is not numeric: " << ln[i+1] << std::endl;
                    return false;
                }
                
                startPosition=tempInt;
                
                start=true;
                i++;
            }
            else if (tag.compare("END")==0){
                
                int tempInt;
                if (!stringToInt(ln[i+1], tempInt)){
                    std::cerr << "Value in External Definition END is not numeric: " << ln[i+1] << std::endl;
                    return false;
                }
                
                
                stopPosition=tempInt;
                stop=true;
                i++;
            }
            else if (tag.compare("TRACE")==0){
                std::string& trace=ln[i+1];
                std::vector<std::string> traces;
                split_line(traces, trace);
                
                //Check trace size;
                size_t expectedSize=stopPosition-startPosition+1;
                
                if (traces.size()!=expectedSize){
                    std::cerr << "Expected external definition trace size: " << expectedSize << "\n but got " << trace << std::endl; 
                    return false;
                }
                path.assign(traces.size(),-2);
                
                for(size_t k=0;k<traces.size();k++){
                    if (info.stateIterByName.count(traces[k])){
                        path[k]=info.stateIterByName[traces[k]];
                    }
                    else{
                        std::cerr << "External definition Trace state name: " << traces[k] << " doesn't exist in the model" << std::endl;
                        return false;
                    }
                }
                state=true;
                i++;
            }
            else{
                std::cerr << "Invalid tag found in sequence external definition: " << ln[i] <<std::endl;
            }
        }
        
        //If everything is defined correctly then return true
        if (start && stop && state){
            for(size_t i=startPosition-1;i<stopPosition;i++){
                size_t state=path[i-(startPosition-1)];
                if (defs.defined(i)){
                    if (!defs[i]->absolute){
                        std::cerr << "Absolute overlaps weighted state definition" << std::endl;
                        
                    }
                    else if (state!=defs[i]->getState()){
                        std::cerr << "Two absolute paths defined for " << i+1 << " position in sequence" <<std::endl;
                    }
                }
                else{
					defs[i] = new (std::nothrow) ExDef;
					defs[i]->setState(state);
                }
            }
            return true;
        }
        
        return false;
    }
    
    //! Parse the weightDef from stringList give the stateInfo
    bool ExDefSequence::_parseWeightDef(stringList& ln, stateInfo& info){
        bool start(false);
        bool stop(false);
        bool state(false);
        bool value(false);
        bool valType(false);
        
        size_t startPosition(0);
        size_t stopPosition(0);
        
        std::set<size_t> definedStates;
        std::set<size_t>::iterator setIterator;
        double val(0);
        
        
        for (size_t line_iter=2; line_iter < ln.size();line_iter++){
            std::string& tag=ln[line_iter];
            
            if (ln.size()<=line_iter+1){
                std::cerr << "External Definition for Weighted Def is missing values: " << ln[line_iter] << std::endl;
                return false;
            }
            
            if (tag.compare("START")==0){
                
                
                int tempInt;
                if (!stringToInt(ln[line_iter+1], tempInt)){
                    std::cerr << "Value in External Definition START is not numeric: " << ln[line_iter+1] << std::endl;
                    return false;
                }
                
                
                startPosition=tempInt;
                start=true;
                line_iter++;
            }
            else if (tag.compare("END")==0){
                
                size_t tempInt;
                if (!stringToInt(ln[line_iter+1], tempInt)){
                    std::cerr << "Value in External Definition END is not numeric: " << ln[line_iter+1] << std::endl;
                    return false;
                }
                
                
                stopPosition=tempInt;
                stop=true;
                line_iter++;
            }
            else if (tag.compare("STATE_NAME")==0){
                definedStates.insert(info.stateIterByName[ln[line_iter+1]]);
                state=true;
                line_iter++;
            }
            else if (tag.compare("STATE_LABEL")==0){
                std::vector<size_t>& temp=info.stateIterByLabel[ln[line_iter+1]];
                for(size_t temp_iter=0; temp_iter < temp.size();temp_iter++){
                    definedStates.insert(temp[temp_iter]);
                }
                state=true;
                line_iter++;
            }
            else if (tag.compare("STATE_GFF")==0){
                std::vector<size_t>& temp = info.stateIterByGff[ln[line_iter+1]];
                for(size_t temp_iter=0; temp_iter < temp.size(); temp_iter++){
                    definedStates.insert(temp[temp_iter]);
                }
                state=true;
                line_iter++;
            }
            else if (tag.compare("VALUE")==0){
                
                double tempValue;
                if(!stringToDouble(ln[line_iter+1], tempValue)){
                    std::cerr << "VALUE couldn't be converted to numerical value: " << ln[line_iter+1] << std::endl;
                }
                
                val=tempValue;
                value=true;
                line_iter++;
            }
            else if (tag.compare("VALUE_TYPE")==0){
                std::string &type=ln[line_iter+1];
                if (type.compare("P(X)")==0){
                    val=log(val);
                    valType=true;
                }
                else if (type.compare("LOG")==0){
                    valType=true;
                }
                else{
                    valType=false;
                }
                line_iter++;
            }
            else{
                std::cerr << "Invalid tag found in Sequence external definition: " << ln[line_iter] << std::endl;
            }
        }
        
        //Check to see that there are valid states defined in the set
        if (definedStates.size()>0){
            state=true;
        }
        else{
            std::cerr << "No valid states defined by External definition" <<std::endl;
            state=false;
        }
        
        //If everything is define correctly then return true
        if (start && stop && state && value && valType){
            for (size_t position=startPosition-1; position < stopPosition; position++){
                //Already have a defined external def at position
                if (defs.defined(position)){
                    if (defs[position]->absolute){
                        std::cerr << "Can't add weight to absolute external definition" << std::endl;
                        return false;
                    }
                }
                else{
                    defs[position] = new(std::nothrow) weightDef(info.stateIterByName.size());
                    
                    if (defs[position]==NULL){
                        std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                        exit(1);
                    }
                }
                
                // Add states and values to definition
                for (setIterator=definedStates.begin();setIterator!=definedStates.end();setIterator++){
                    defs[position]->assignWeight(*setIterator, val);
                }
            }
            return true;
        }
        return false;
    }

}