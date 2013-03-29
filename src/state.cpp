//state.cpp
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

#include "state.h"
namespace StochHMM{
	
    //!Create a state object
    state::state():endi(NULL), stateIterator(SIZE_MAX){
        transi = new (std::nothrow) std::vector<transition*>;
    }
	
    //!Create a state from string
    //! Parses the string to create a state
    //! \param txt String representation of a state
    //! \param names Names of all the states
    //! \param trcks Tracks defined for model
    //! \param wts Pointer to all weight defined for the model
    //! \param funcs State functions defined for the model
    state::state(std::string& txt, stringList& names,tracks& trcks, weights* wts, StateFuncs* funcs):endi(NULL), stateIterator(SIZE_MAX){
        
        //endi=new transition(STANDARD);
        transi = new std::vector<transition*>;
        parse(txt,names,trcks, wts,funcs);
    }
    
    state::~state(){
        delete transi;
        transi=NULL;
    }
    
    
	
	//state::state(int size){
	//	name="";
	//	gff="";
	//	label="";
	//	endi=new transitions(STANDARD);
	//	for (int i=0;i<size;i++){
	//		//trans.push_back(0);
	//		//log_trans.push_back(MIN);
	//		transi.push_back(initial);
	//	}
	//}
    
    
    //!Parses state from string
    //! \param txt String of state definition
    //! \param names StringList with names of other states. Used to identify position of transitions in transition
    //! In future, may want to organize transitions in non-linear fashion
    bool state::parse(std::string& txt, stringList& names,tracks& trks, weights* wts, StateFuncs* funcs){
        size_t stateHeaderInfo = txt.find("STATE:");
        size_t transitionsInfo = txt.find("TRANSITION:");
        size_t emissionInfo    = txt.find("EMISSION:");
        stringList lines;
        
        if (stateHeaderInfo==std::string::npos){
            std::cerr << "Couldn't identify state Header Information. Please check the formatting.  The supplied text was : " << txt << std::endl;
            return false;
        }
        
        if (transitionsInfo == std::string::npos){
            std::cerr << "Couldn't identify state Transitions Information.  Please check the formatting.  The supplied text was : " << txt << std::endl;
            return false;
        }
        
        //Extract and Parse Header Information
        std::string header = txt.substr(stateHeaderInfo,transitionsInfo-stateHeaderInfo);
        if (!_parseHeader(header)){
            return false;
        }
		
        
        //Extract and Parse Transition Information
        std::string trans = (emissionInfo==std::string::npos) ? txt.substr(transitionsInfo) : txt.substr(transitionsInfo, emissionInfo - transitionsInfo);
        //std::cout << trans << std::endl;
        if (!_parseTransition(trans, names, trks, wts, funcs)){
            return false;
        }
        
        
        //Check emissions existence  (only INIT state can have no emission)
        if (emissionInfo != std::string::npos){
            std::string emmis = txt.substr(emissionInfo);
            if (!_parseEmission(emmis, names, trks, wts, funcs)){
                std::cerr << "Couldn't parse the emissions for state: " << name << std::endl;
                return false;
            }
        }
        else{
            //Check name is INIT.   Only INIT state can not have an emission
            if (name.compare("INIT")!=0){
                std::cerr << "No emission defined for State: " << name << std::endl;
                return false;
            }
        }
        
        return true;
    }
    
    //!Parse the state header information and assign to class variables
    //!Required header information includes then NAME and PATH_LABEL
    //!Optonal header information includes GFF_DESC
    //! \param txt String representation of state header
    bool state::_parseHeader(std::string& txt){
        
        clear_whitespace(txt,"\t ");
        stringList lst = splitString(txt,"\n:");
        
        size_t idx;
        
        //parse NAME  (REQUIRED)
        if (lst.contains("NAME")){
            idx= lst.indexOf("NAME");
            if (idx+1 < lst.size()){
                idx++;
                name=lst[idx];
            }
            else{
                std::cerr << "The states NAME couldn't be parsed from " << txt << std::endl;
                return false;
            }
            
            
            if (name.compare("INIT")==0){
                return true;
            }
        }
        else{
            std::cerr << "No NAME found in state header information. Please check formatting. Following is header information recieved: " << txt << std::endl;
            return false;
        }
        
        //parse PATH_LABEL  (REQUIRED)
        if (lst.contains("PATH")){
            idx= lst.indexOf("PATH");
            if (idx+1 < lst.size()){
                idx++;
                label=lst[idx];
            }
            else{
                std::cerr <<  "The states PATH_LABEL couldn't be parsed from " << txt << std::endl;
                return false;
            }
            
        }
        else{
            std::cerr << "No PATH_LABEL found in state header information. Please check formatting. Following is header information recieved: " << txt << std::endl;
            return false;
        }
        
        //Parse GFF_DESC  (Optional)
        if (lst.contains("GFF")){
            idx= lst.indexOf("GFF");
            if (idx+1 < lst.size()){
                idx++;
                gff=lst[idx];
            }
            else{
                std::cerr <<  "The states GFF_DESC couldn't be parsed from " << txt << std::endl;
                return false;
            }
        }
		
        return true;
    }
    
    
    //!Parses the transitions of the state from a text string
    //! \param txt Text representation of the transitions
    //! \param name StringList of all state names defined in the model
    //! \param trks Reference to tracks for the model
    //! \param wts Weight defined in the model
    //! \param funcs State functions defined for the model
    bool state::_parseTransition(std::string& txt, stringList& names, tracks& trks, weights* wts, StateFuncs* funcs){
        //SPLIT UP TRANSITIONS AND APPLY SEPARATELY
        stringList lst;
        lst.splitND(txt,"TRANSITION:");
        
        
        for(size_t iter=0;iter<lst.size();iter++){
            
            stringList line=splitString(lst[iter],"\n");
            line.removeLWS("\t ");
            line.removeComments();
            //Line 1 defines transition type
            
            clear_whitespace(line[0],"\t");
            stringList head=splitString(line[0],":");
            
            //DETERMINE TYPE
            transType tp;
            if (head.contains("STANDARD")){tp=STANDARD;}
            else if (head.contains("DURATION")){tp=DURATION;}
            else if (head.contains("LEXICAL")){tp=LEXICAL;}
            else{
                std::cerr << "Unrecognized Transition type(STANDARD, DURATION, LEXICAL):" << txt << std::endl;
                return false;
                //Error not recognized transition type
            }
            
            //DETERMINE VALUE TYPE
            valueType valtyp;
            if (head.contains("P(X)")){valtyp=PROBABILITY;}
            else if (head.contains("LOG")){valtyp=LOG_PROB;}
            else if (head.contains("COUNTS")){valtyp=COUNTS;}
            else {
                std::cerr << "Unrecognized Transition value type( P(X), LOG, COUNTS): " << txt << std::endl;
                return false;
                //Error not recognized value type
            }
            
            
            
            if (tp==STANDARD){
                //Process each following from line
                for(size_t iter=1;iter<line.size();iter++){
                    transition* temp = new(std::nothrow) transition(tp);
                    
                    if (temp==NULL){
                        std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                        exit(1);
                    }
                    
                    if (!temp->parse(line[iter], names, valtyp, trks, wts, funcs)){
                        std::cerr << "Couldn't parse Transition " << std::endl;
                    }
                    
                    std::string transName = temp->getName();
                    
                    if (transName=="END"){
                        endi = temp;
                    }
                    else{
                        //size_t idx = names.indexOf(transName);
                        //idx--;
                        addTransition(temp);
                    }
                }
                
            }
			
            else if (tp==DURATION || tp==LEXICAL ){
                transition* temp = new(std::nothrow) transition(tp);
                
                if (temp==NULL){
                    std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                    exit(1);
                }
                
                if (!temp->parse(line,names,valtyp,trks,wts,funcs)){
                    std::cerr << "Couldn't parse Transition "<< std::endl;
                    return false;
                }
                std::string transName = temp->getName();
                
                
                //size_t idx = names.indexOf(transName);
                //idx--;
                addTransition(temp);
            }
        }
        return true;
    }
    
    
    //!Parses the emission for a state from a string
    //! \param txt String representation of emissions
    //! \param names stringList of all state names defined in the model
    //! \param trks Tracks defined in the model
    //! \param wts Weight defined of the model
    //! \param funcs StateFunction defined for the model
    
    bool state::_parseEmission(std::string& txt, stringList& names, tracks& trks, weights* wts, StateFuncs* funcs){
        stringList lst;
        lst.splitND(txt,"EMISSION:");
        //lst.print();
        
        for(size_t iter=0; iter<lst.size();iter++){
            emm* temp = new(std::nothrow) emm;
            
            if (temp==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
            
            if (!temp->parse(lst[iter],trks,wts,funcs)){
                return false;
            }
            emission.push_back(temp);
        }
        
        return true;
    }
    
    //Print the string representation of the state to the stdout
    void state::print(){
        std::cout<< stringify() << std::endl;
    }
    
    
    //Get a string representation for the state
    //! \return std::string
    std::string state::stringify(){
        std::string stateString;
        stateString+="STATE:\n";
        stateString+="\tNAME:\t" + name + "\n";
        
        if (name.compare("INIT")!=0){
            if (!gff.empty()){ stateString+="\tGFF_DESC:\t" + gff + "\n";}
            stateString+="\tPATH_LABEL:\t" + label + "\n";
        }
        
        std::string standardString;
        std::string distribString;
        std::string lexicalString;
        
        //Get Transitions Standard, Distribution, Lexical
        for(size_t i=0;i<transi->size();i++){
            if ((*transi)[i]==NULL){ continue;}
            
            transType tp = (*transi)[i]->getTransitionType();
            if (tp == STANDARD){
                if (standardString.empty()){
                    standardString+="TRANSITION:\tSTANDARD:\tLOG\n";
                }
                standardString+=(*transi)[i]->stringify() + "\n";
            }
            else if (tp == DURATION){
                distribString+="TRANSITION:\tDURATION:\tLOG\n";
                distribString+=(*transi)[i]->stringify();
                
            }
            else if (tp == LEXICAL){
                if ((*transi)[i]->LexFunctionDefined()){
                    lexicalString+="TRANSITION:\tLEXICAL:\tFUNCTION:\t";
                    lexicalString+=(*transi)[i]->getLexicalFunctionName();
                    lexicalString+="\n";
                }
                else{
                    lexicalString+="TRANSITION:\tLEXICAL:\tLOG\n";
                    lexicalString+=(*transi)[i]->stringify();
                }
                
            }
        }
        
        
        //Process Transition to Ending;
        if (endi!=NULL){
            if (standardString.empty()){
                standardString+="TRANSITIONS:\tSTANDARD:\tLOG\n";
            }
            standardString+=endi->stringify();
        }
        
        
        if (name.compare("INIT")==0){
            stateString+=standardString + "\n\n";
            return stateString;
        }
        else{
            if (!standardString.empty()){
                stateString+=standardString + "\n";
            }
            
            if (!distribString.empty()){
                stateString+=distribString + "\n";
            }
            
            if (!lexicalString.empty()){
                stateString+=lexicalString + "\n";
            }
        }
        
        
        
        //Print Emissions
        for(size_t i=0;i<emission.size();i++){
            stateString+=emission[i]->stringify();
        }
        
        return stateString;
    }
	
    //! Get the emission value from a state at a particular position in the sequence
    //! \param seqs Sequences the model is analyzing
    //! \param iter Position in the sequence
    double state::get_emission_prob(sequences &seqs, size_t iter){
        double value(emission[0]->get_emission(seqs,iter));
        for(size_t i=1;i<emission.size();i++){
            value+=emission[i]->get_emission(seqs,iter);  //Pass sequences type to get_emission for each emission in the state
        }
        return value;
    }
    
    
    //! Get the transition value from a state at a particular position in the sequence
    //! \param seqs Sequences the model is analyzing
    //! \param to State that transition is being calculated to
    //! \param iter Position in the sequence
    double state::get_transition_prob(sequences &seqs, size_t to, size_t iter){
        double value;
        
        if ((*transi)[to]==NULL){
            return -INFINITY;
        }
        else if ((*transi)[to]->transition_type==STANDARD){
            value = (*transi)[to]->log_trans;
        }
        else{
            std::cerr << "Need to implement this functionality" <<std::endl;
            value = (*transi)[to]->getTransition(iter,&seqs);
            
        }
        return value;
    }
    
    //! Get the log probability transitioning to end from the state
    double state::getEndTrans(){
        if (endi==NULL){
            return -INFINITY;
        }
        return endi->log_trans;
    }
    
    //TODO: complete the checkLabels function
    //! Checks the label tags for traceback and combine identifiers
    void state::checkLabels(std::set<std::string>& labels, std::set<std::string>& gff, std::set<std::string>& name){
        for(size_t i = 0;i<(*transi).size();i++){
            transitionFuncParam* func = (*transi)[i]->getExtFunction();
            if (func!=NULL){
                
                tracebackIdentifier tt = func->getTracebackType();
                std::string tn = func->getTracebackName();
                
                //Traceback Identifier
                if (tt != START_INIT && tt != DIFF_STATE){
                    if (tt == STATE_NAME){
                        if (!name.count(tn)){
                            std::cerr << "Model Function or Distribution traceback definition contains Unknown State Name:\t" << tn << std::endl;
                        }
                    }
                    else if (tt == STATE_GFF){
                        if (!gff.count(tn)){
                            std::cerr << "Model function or distribution traceback definition contains Unknown GFF Descrition:\t" << tn << std::endl;
                        }
                    }
                    else if (tt == STATE_LABEL){
                        if (!labels.count(tn)){
                            std::cerr << "Model function or distribution traceback definition contains Unknown Path Label:\t" << tn << std::endl;
                        }
                    }
                }
                
                
                combineIdentifier ct = func->getCombineType();
                std::string cn = func->getCombineName();
                
                //Combine Identifiers
                if (ct != FULL){
                    if (ct == STATENAME){
                        if (!name.count(cn)){
                            std::cerr << "Traceback Combine definition contains Unknown State Name:\t" << cn << std::endl;
                        }
                    }
                    else if (ct == STATEGFF){
                        if (!gff.count(cn)){
                            std::cerr << "Traceback Combine definition contains Unknown GFF Descrition:\t" << cn << std::endl;
                        }
                    }
                    else if (ct == STATELABEL){
                        if (!labels.count(cn)){
                            std::cerr << "Traceback Combine definition contains Unknown Path Label:\t" << cn << std::endl;
                        }
                    }
                }
				
            }
        }
        
        return;
    }
    
    
    /* On initial import of the states they are pushed on the transi vector in
     the order written in model.   However, the analysis requires that they be
     in the particular position defined by state iterator.
     
     This function puts the transitions in the proper order for analysis
     */
    void state::_finalizeTransitions(std::map<std::string,state*>& state_index){
		
		//Get size # of states, but correct by -1 because
		//initial state will be kept separate.
        size_t number_of_states = state_index.size();
        std::vector<transition*>* fixed_trans = new std::vector<transition*>(number_of_states-1,NULL);
        
        //Find the proper place for the transition and put it in the correct position
        for(size_t i = 0; i < transi->size(); i++){
            transition* temp = (*transi)[i];
            std::string name = temp->getName();
            state* st = state_index[name];
			if (st == NULL){
				std::cerr << "State: " << name << " was declared but not defined in the model." << std::endl;
				exit(2);
			}
            size_t index = st->getIterator();
            (*fixed_trans)[index]=temp;
            (*transi)[i]=NULL;
        }
        
        delete transi;  //Don't need the old transition vector anymore
        transi = fixed_trans;
        return;
    }
	
	
	bool state::hasComplexEmission(){
		for(size_t i=0; i < emission.size() ; ++i){
			if (emission[i]->isComplex()){
				return true;
			}
		}
		return false;
	}
	
	
	
	bool state::hasComplexTransition(){
		for(size_t i=0;i<(*transi).size();++i){
			
			if ((*transi)[i]==NULL){
				continue;
			}
			
			if ((*transi)[i]->isComplex()){
				return true;
			}
		}
		return false;
	}
	
	bool state::isSimple(){
		if (!hasComplexEmission() && ! hasComplexTransition()){
			return true;
		}
		return false;
	}
	
}
