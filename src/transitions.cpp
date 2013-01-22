//transitions.cpp
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

#include "transitions.h"
namespace StochHMM{

#pragma mark Initiation
    
    //!Create a transition of a certain type
    //! \param type enum transType (STANDARD, DURATION, LEXICAL)
    transition::transition(transType type){
        transition_type=type;
        traceback_identifier=DIFF_STATE;
        log_trans=-INFINITY;        
        extendedValue=-INFINITY;
		function=false;
        func	= NULL;
        lexFunc	= NULL;
        toState	= NULL;
		func	= NULL;
    }
    
    //!Create a transition of a certain type
    //! \param type Type of transType
    //! \param valType Type of Value Transition uses
    //! \param Is distributions defined as survival 
    transition::transition(transType type, valueType valtyp, bool survival){
        transition_type=type;
        traceback_identifier = DIFF_STATE;
        log_trans = -INFINITY;
        func=NULL;
        lexFunc=NULL;
        function=false;
        toState=NULL;
		func = NULL;
    }
    
    //! Parse the transition User and Standard Transition from String
    //! to create a transition for the state
    //! \param txt String representation of the transition
    //! \param names StringList of all the states
    //! \param valtyp Value type used (log, p(x)...)
    //! \param trks  Tracks of the model
    //! \param wts  Weight defined in the model
    //! \param funcs State Functions created by the user
    bool transition::parse(std::string& txt,stringList& names, valueType valtyp, tracks& trks, weights* wts , StateFuncs* funcs){
        
        if (!_processTags(txt,trks, wts, funcs)){
            std::cerr << "Couldn't properly process tag from: "<< txt << std::endl;
            return false;
        }
        
        if (transition_type == STANDARD){
            if (!_parseStandard(txt,names, valtyp)){
                std::cerr << "Couldn't parse the Standard Transition" << std::endl;
            }
        }
        
        return true;
    }
    
    //!Parse the Transition Lexical and User from String List
    //! \param txt String representation of the transition used to create the transition from model file
    //! \param names StringList of all the states
    //! \param valtyp Value type used (log, p(x)...)
    //! \param trks  Tracks of the model
    //! \param wts  Weight defined in the model
    //! \param funcs State Functions created by the user
    bool transition::parse(stringList& txt,stringList& names, valueType valtyp, tracks& trks, weights* wts , StateFuncs* funcs){
        
        //txt.print();
        
        if (!_processTags(txt[1],trks, wts, funcs)){
            std::cerr << "Couldn't properly process tag from: "<< txt[1] << std::endl;
            return false;
        }
        
        if (transition_type == DURATION ){
            if (!_parseDuration(txt, names, valtyp)){
                std::cerr << "Couldn't parse Duration Transition" <<std::endl;
                return false;
            }
        }
        else if (transition_type == LEXICAL){
            if (!_parseLexical(txt, names, valtyp, trks, funcs)){
                std::cerr << "Couldn't parse Lexical Transition" << std::endl;
                return false;
            }
        }
        
        return true;
    }
    
    
    //TODO: Check for errors within TAG information
    
    //Process the addition tags after the transition  
    //Tags may describe functions or weights to apply to the transition
    //Can also describe how a traceback is to be done and processed
    bool transition::_processTags(std::string& txt, tracks& trks, weights* wts, StateFuncs* funcs){
        stringList lst = extractTag(txt);
        
        if (lst.size() == 0){
            return true;
        }
        
        func=new(std::nothrow) transitionFuncParam();
        
        if (func==NULL){
            std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
            exit(1);
        }
        
        if (!func->parse(lst,trks, wts,funcs)){
            std::cerr << "Couldn't parse Transition Tag information" << std::endl;
            return false;
        }
        
        
        std::string trackName = func->getTrackName();
        
        //Does it contain the name
        
        return true;
    }
    
    //Parse a simple standard transition from a model file
    bool transition::_parseStandard(std::string& txt, stringList& names, valueType valtyp){
        
        //Process Transition
        stringList line;
        line.splitString(txt,"\t:");
        
        if (line.size()<2){
            std::cerr << "Line should contain 2 values (STATE  VALUE).  Couldn't parse:\n" << txt << std::endl;
            return false;
        }
        
        stateName= line[0];
        
        
        if (!names.contains(stateName) && stateName!="END"){
            std::cerr << "Tried to create a transition in the model to state named : " << stateName << " but there is no state with that name.  Please check the formatting. \n";
            return false;
        }
        
        if(valtyp==PROBABILITY){
            double tempValue;
            
            if (!stringToDouble(line[1], tempValue)){
                std::cerr << "Probability value not numeric: " << line[1] << std::endl;
                return false;
            }
            
            log_trans = log(tempValue);
        }
        else if (valtyp == LOG_PROB){
            
            double tempValue;
            
            if (!stringToDouble(line[1], tempValue)){
                std::cerr << "Log Probability value not numeric: " << line[1] << std::endl;
                return false;
            }
            
            log_trans = tempValue;
            
        }
        
        
        return true;
    }

    
    
    //Parse the user-defined distribution from the model file
    bool transition::_parseDuration(stringList& txt, stringList& names, valueType valtyp){
        stringList line;
        line.splitString(txt[1],"\t:,");
        
        if (line.size()<1){
            std::cerr << "Couldn't parse the transition STATE information: "<< txt[1] << std::endl;
        }
        
        stateName = line[0];
        
        //Determine Traceback options
        if (line.contains("TO_LABEL") || line.contains("TB->LABEL")){
            size_t idx= (line.contains("TO_LABEL")) ? line.indexOf("TO_LABEL") : line.indexOf("TB->LABEL");
            idx++;
            traceback_identifier = STATE_LABEL;
            traceback_string = line[idx];
        }
        else if (line.contains("TO_GFF") || line.contains("TB->GFF")){
            size_t idx= (line.contains("TO_GFF")) ? line.indexOf("TO_GFF") : line.indexOf("TB->GFF");
            idx++;
            traceback_identifier = STATE_GFF;
            traceback_string = line[idx];
        }
        else if (line.contains("TO_STATE") || line.contains("TB->STATE")){
            size_t idx= (line.contains("TO_STATE")) ? line.indexOf("TO_STATE") : line.indexOf("TB->STATE");
            idx++;
            traceback_identifier = STATE_NAME;
            traceback_string = line[idx];
        }
        else if (line.contains("TO_START") || line.contains("TB->START")){
            traceback_identifier = START_INIT;
        }
        else{
            traceback_identifier = DIFF_STATE;
        }
        
        //Process Distribution
        //from line 2 to end is distribution
        distribution = new(std::nothrow) std::vector<double>;
        
        if (distribution==NULL){
            std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
            exit(1);
        }
        
        for(int j=2;j<txt.size();j++){
            
            stringList ln;
            ln.splitString(txt[j],"\t");
            if (ln.size()!=2){
                std::cerr << "More than 2 values on the line: " << ln.stringify() << std::endl;
                return false;
            }
            
            int index(0);
            if (!stringToInt(ln[0], index)){
                std::cerr << "Distribution position value not numeric:\t" << ln[0] << std::endl;
                return false;
            }
            
            double value;
            if (!stringToDouble(ln[1], value)){
                std::cerr << "Distribution position value not numeric:\t" << ln[0] << std::endl;
                return false;
            }
            
            if (valtyp==PROBABILITY){
                value=log(value);
            }
            
            
            //Fill the lower range up to index with the first value
            //Extend values to lower
            while (distribution->size()<index-1){
                distribution->push_back(value);
            }
    #ifdef DEBUG_MODEL
            //cout << value <<endl;
    #endif
            distribution->push_back(value);
        }

        extendedValue = distribution->back();
        return true;

    }
    
    //Parse the lexical transition from the model file
    bool transition::_parseLexical(stringList& txt, stringList& names, valueType valtyp, tracks& trks, StateFuncs* funcs){
        
        //Process Transition
        stringList line;
        size_t idx;
        std::string functionName("");
        
        if (txt.contains("FUNCTION")){
            idx=txt.indexOf("FUNCTION");
            function=true;
            
            if (idx+1 < txt.size()){
                functionName=txt[idx+1];
            }
            else{
                std::cerr << "Couldn't parse the function name from the Lexical Transition:\n" << txt.stringify() << std::endl;
            }
            
        }
        
        
        idx=txt.indexOf("TRANSITION");
        
        if (idx+1 < txt.size()){
            idx++;
            line.splitString(txt[idx],"\t:,");
            
            if (line.size()>0){
                stateName = line[0];
            }
            else{
                std::cerr << "Couldn't parse the STATE for the Lexical transition:\n" << txt.stringify() << std::endl;
            }
        }
        else{
            std::cerr << "Couldn't parse the Lexical Transition:\n" << txt.stringify() << std::endl;
            return false;
        }
        
        
        if (!names.contains(stateName)){
            std::cerr << "Lexical transition defined transition to state named : " << stateName << " However, there doesn't appear to be any state with that name\n";
            return false;
        }
        
        //remaining tracks and Orders then set Track
        std::vector<track*> tempTracks;
        for(size_t i=1;i<line.size();i++){
            track* tk = trks.getTrack(line[i]);
            if (tk==NULL){
                std::cerr << "Lexical Transition tried to add a track named: " << line[i] << " . However, there isn't a matching track in the model.  Please check to model formatting.\n";
                return false;
            }
            else{
                tempTracks.push_back(tk);
            }
        }
        
        
        if (function){  //Lexical Function
            lexFunc=new(std::nothrow) emissionFuncParam(functionName,funcs,tempTracks[0]);
            
            if (lexFunc==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
        }
        
        else{   //Lexical Tables
            //Get Order Information
            std::vector<int> tempOrder;
            if (txt.contains("ORDER")){
                idx=txt.indexOf("ORDER");
            }
            else{
                std::cerr << "Unable to locate ORDER: for lexical transition to State Name: " << stateName << std::endl;
                return false;
            }
            line.splitString(txt[idx],"\t:,");
            
            size_t ambIdx;
            bool containsAmbig=line.contains("AMBIGUOUS");
            if (containsAmbig){
                ambIdx=line.indexOf("AMBIGUOUS");
                
                for(size_t i=1;i<ambIdx;i++){
                    
                    int tempInt;
                    if (!stringToInt(line[i], tempInt)){
                        std::cerr << "Lexical Transition Order is not numeric" << line[i] << std::endl;
                        return false;
                    }
                    
                    tempOrder.push_back(tempInt);
                }
            }
            else{
                for(size_t i=1;i<line.size();i++){
                    
                    int tempInt;
                    if (!stringToInt(line[i], tempInt)){
                        std::cerr << "Lexical Transition Order is not numeric" << line[i] << std::endl;
                        return false;
                    }
                    
                    tempOrder.push_back(tempInt);
                }
            }
            
            
            
            if (tempOrder.size() == tempTracks.size()){
                for(size_t i=0;i<tempOrder.size();i++){
                    scoreTable.addTrack(tempTracks[i], tempOrder[i]);
                }
            }
            else{
                std::cerr << "Different number of tracks and orders parsed in LEXICAL TRANSITION to state: " << stateName << " Check the formatting of the Lexical Transition" << std::endl;
                return false;
            }
            
            //Parse Ambiguous Tag Info
            if (containsAmbig){
                ambIdx++;
                if (line.size()<=ambIdx){
                    std::cerr << "No scoring type after AMBIGUOUS label\nAssuming AVG\n";
                    scoreTable.setUnkScoreType(AVERAGE_SCORE);
                }
                else if (line[ambIdx].compare("AVG")==0){scoreTable.setUnkScoreType(AVERAGE_SCORE);}
                else if (line[ambIdx].compare("MAX")==0){scoreTable.setUnkScoreType(HIGHEST_SCORE);}
                else if (line[ambIdx].compare("MIN")==0){scoreTable.setUnkScoreType(LOWEST_SCORE);}
                else if (line[ambIdx].compare("P(X)")==0)
                {
                    scoreTable.setUnkScoreType(DEFINED_SCORE);
                    ambIdx++;
                    
                    double tempValue;
                    if (!stringToDouble(line[ambIdx], tempValue)){
                        std::cerr << "Ambiguous Value couldn't be parsed: "<< line[ambIdx] << std::endl;
                        return false;
                    }
                    
                    scoreTable.setUnkScore(log(tempValue));
                }
                else if (line[ambIdx].compare("LOG")==0){
                    scoreTable.setUnkScoreType(DEFINED_SCORE);
                    ambIdx++;
                    
                    double tempValue;
                    if (!stringToDouble(line[ambIdx], tempValue)){
                        std::cerr << "Ambiguous Value couldn't be parsed: "<< line[ambIdx] << std::endl;
                        return false;
                    }
                    
                    scoreTable.setUnkScore(tempValue);
                }
            }
            
            
            
            //Get Tables
            int expectedColumns = 1;
            int expectedRows = 1;
            for(size_t i = 0; i<scoreTable.getNumberOfAlphabets(); i++){
                expectedColumns*=scoreTable.getAlphaSize(i);
                expectedRows*=POWER[scoreTable.getOrder(i)][scoreTable.getAlphaSize(i)-1];
            }
            
            std::vector<std::vector<double> >* log_prob = scoreTable.getLogProbabilityTable();
            std::vector<std::vector<double> >* prob = scoreTable.getProbabilityTable();
            std::vector<std::vector<double> >* counts = scoreTable.getCountsTable();
            
            
            idx++;
            
            for (size_t iter = idx; iter< txt.size();iter++){
                
                //If it's the first line check for a '#' indicating that the column header is present
                if (iter==idx && txt[idx][0]=='@'){
                    continue;
                }
                
                line.splitString(txt[iter],"\t ");
                
                //Check for Row header
                if (line[0][0]=='@'){
                    line.pop_ith(0);
                }
                
                //line.splitString(txt[iter],"\t ");
                std::vector<double> temp = line.toVecDouble();
                if (temp.size() != expectedColumns){
                    std::cerr << "The following line couldn't be parsed into the required number of columns.   Expected Columns: " << expectedColumns << "\n The line appears as: "  << txt[iter] << std::endl;
                    return false;
                }
                else{
                    if (valtyp == PROBABILITY){
                        prob->push_back(temp);
                        logVector(temp);
                        log_prob->push_back(temp);
                    }
                    else if (valtyp == LOG_PROB){
                        log_prob->push_back(temp);
                        expVector(temp);
                        prob->push_back(temp);
                    }
                    else if (valtyp == COUNTS){
                        counts->push_back(temp);
                        probVector(temp);
                        prob->push_back(temp);
                        logVector(temp);
                        log_prob->push_back(temp);
                    }
                    
                }
            }
            
            if (log_prob->size() != expectedRows){
                std::cerr << " The Lexical table doesn't contain enough rows.  Expected Rows: " << expectedRows << " \n Please check the Lexical Table and formatting\n";
                return false;
            }

        }
        
        return true;
    }

    
    
    #pragma mark Accessor
    
    //!Print the transition to stdout
    void transition::print(){
        std::cout << stringify() << std::endl;
    }
    
    //!Convert the transition to a string representation as shown in the model file
    //! \return std::string Representation of the transition
    std::string transition::stringify(){
        std::string transString;
        if (transition_type == STANDARD){
            transString+="\t" + stateName + ":\t" + double_to_string(log_trans);
            if (func!=NULL){
                transString+="\t" + func->stringify();
            }
        }
        else if (transition_type == DURATION){
            transString+="\t" + stateName + ":\t" ;
            transString+= (traceback_identifier==STATE_NAME) ? "TO_STATE:\t" + traceback_string :
                          (traceback_identifier==STATE_LABEL) ? "TO_LABEL:\t" + traceback_string :
                          (traceback_identifier==STATE_GFF) ? "TO_GFF:\t" + traceback_string : "DIFF_STATE" ;
            if (func!=NULL){
                transString+="\t" + func->stringify();
            }
            
            transString+="\n";
            for(size_t i=0;i<distribution->size();i++){
                transString+= "\t\t" + int_to_string((int) i+1) + "\t" + double_to_string((*distribution)[i]) + "\n";
            }
        }
        else if (transition_type == LEXICAL){
            transString+="\t" + stateName + ":\t";
            
            
            if (function){
                transString+=lexFunc->getTrackName();
                if (func!=NULL){
                    transString+="\t" + func->stringify();
                }
            }
            else{
                for(size_t i=0;i<scoreTable.trackSize();i++){
                    if (i>0){
                        transString+=",";
                    }
                    transString+=scoreTable.getTrack(i)->getName();
                }
                
                
                if (func!=NULL){
                    transString+="\t" + func->stringify();
                }
                
                transString+="\n\t\tORDER:\t";
                
                for(size_t i=0;i<scoreTable.trackSize();i++){
                    if (i>0){
                        transString+=",";
                    }
                    transString+=int_to_string(scoreTable.getOrder(i));
                }
                
                unknownCharScoringType ambtemp = scoreTable.getAmbScoringType();
                
                if (ambtemp!=NO_SCORE){
                    transString+="\tAMBIGUOUS:\t";
                    transString+=(ambtemp==HIGHEST_SCORE)? "MAX":
                    (ambtemp==LOWEST_SCORE)? "MIN":
                    (ambtemp==AVERAGE_SCORE)? "AVG": "LOG:" + double_to_string(scoreTable.getAmbDefinedScore());
                }
                
                transString+="\n";
                
                
                //Print the probability table
                //Output Column Headers
                transString+=scoreTable.stringify();

            }
        }
        
        return transString;
    }
    
    //!Calculate the transition probability for a transition at a particular position in the sequence
    //! \param pos Position in sequences that we are determining the transition
    //! \param seqs Pointer to sequences, used in determining transition
    //! \return double Score of the transition
    double transition::getTransition(size_t pos,sequences* seqs){
         if (transition_type==STANDARD){
             return log_trans;
         }
         else if (transition_type==DURATION){
             
             if (pos >= distribution->size()){   // Check that size corresponds to zero index or one index.....
                 return extendedValue;
             }
             else{
                 return (*distribution)[pos-1];
             }
         }
         else if (transition_type==LEXICAL){
             
             double value;
             
             if (function){
                 value = lexFunc->evaluate(*seqs, pos);
             }
             else{
                 value = scoreTable.getValue(*seqs,pos);
             }
            return value;
         }
         else{
            return -INFINITY; 
         }
     }
    
}
