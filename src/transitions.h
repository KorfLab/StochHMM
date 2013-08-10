///transitions.h
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
#ifndef TRANSITIONS_H
#define TRANSITIONS_H

#include <vector> 
#include <string>
#include "text.h"
#include "track.h"
#include "stochMath.h"
#include "stochTypes.h"
#include "externalFuncs.h"
#include "weight.h"
#include "sequences.h"
#include "lexicalTable.h"
#include <stdlib.h>
namespace StochHMM{
    
    //! \file transitions.h
    /*! Define transition class
     */

    class state;
    
    
    //! \class transition
    /*! \brief Transition information for each state
     
     */
    
class transition{
public:
    transition(); // Constructor
    transition(transType); // Constructor
    transition(transType,valueType,bool); // Constructor
    
    
    friend class model;
    friend class state;
    
    //MUTATORS
    
    //
    bool parse(std::string&, stringList&, valueType valtyp, tracks& ,weights*, StateFuncs*);
    bool parse(stringList&, stringList&, valueType valtyp, tracks&  ,weights*, StateFuncs*);
    
    //! Set the name of the state we are transitioning to
    //! \param txt  Name of the next state
    inline void setName(std::string& txt){stateName=txt;};
    
    //! Set the pointer to the state we are transitioning to
    //! \param st  Pointer to state
    inline void setState(state* st){toState=st;};
    
    //! Set the type of the transition (enum transType)
    //! \param type enum Transtype
    //! \sa enum transType
    inline void setTransType(transType type){transition_type=type;};
    
    //Distributions
    //! Sets teh Traceback Identifier used when doing a traceback during calculation of transition
    //! \param tbIdt enum tracebackIdentifier
    //! \sa enum tracebackIdentifier
    inline void setTB_Identifier(tracebackIdentifier tbIdt){traceback_identifier = tbIdt;};
    
    //! Set the traceback Label/GFF/StateName to use during traceback
    //! \param txt String to use for traceback
    inline void setTB_String(std::string& txt){traceback_string = txt;};
    
    //TODO: Check to see that distribution is survival function
    //! Set the length distribution of the transition
    //! Each position in the vector is a length 0=1, 99=100 length
    //! Distribution should be a survival function but could be user defined
    //! \param dst std::vector of double
    inline void setDistribution(std::vector<double>* dst){distribution=dst;}; //Check that distribution is survival function
    
    //Standard
    //! Set the standard transition probability of the transition
    //! \param value log based 2 value of transitioin probabability 
    inline void setTransProb(double value){log_trans=value;};
    
    
    //////////////////////// Lexical ////////////////////////
    
//    void addTrack(track*,int); //Add track and order
    
//    //TODO: check the table is correct size according to order and track
//    //!Set the lexical probability table from vector of vectors
//    //! \param tbl Vector of Vector of double (table of sequence frequencies (log2 value)
//    inline void setLexicalProb(std::vector<std::vector<double> >* tbl){log_prob=tbl;};
//    
//    //TODO: check ambiguous character scoring ....
//    //!Set the how to score unknown or ambiguous characters in the sequence
//    //! \param type enum unknownCharScoringType
//    inline void setUnkScoreType(unknownCharScoringType type){unknownScoreType=type;};
//    
//    
//    //!Set unknown scoring value
//    //!Value to use for unknown characters in the sequence
//    inline void setUnkScore(double val){unknownDefinedScore=val;};
    
    
    //ACCESSORS
    //! Get the name of state that transition is to
    //! \return std::string Name of the state that transition is to
    inline std::string& getName(){return stateName;};
    
    //! Get pointer to stte that transition is to
    //! \return state* Pointer to state that transitioning to
    inline state* getState(){return toState;};
    
    //inline bool isDefined(){return defined;};
    
    //! Get the transition type
    //! \return transType
    //! \sa transType
    inline transType getTransitionType(){return transition_type;};
    
    //! Get the tracebackIdentifier defined in the transitioin
    //! \return tracebackIdentifier
    //! \sa tracebackIdentifier
    inline tracebackIdentifier getTracebackIdentifier(){return traceback_identifier;};
    
    //! Get the string (GFF/Label/State Name) Traceback is to
    //! \return std::string name that traceback is to
    inline std::string& getTracebackString(){return traceback_string;};
    
    double getTransition(size_t,sequences*);  // get the transition using to and the position trellis
    double get_reduced_order(int,sequences*);
    double getTransition();
    
    //! Get pointer to externalFuncs 
    //! \return externalFuncs* Pointer to external function definition in transition
    inline transitionFuncParam* getExtFunction(){return func;};
    inline bool FunctionDefined(){if(func!=NULL){return true;} else {return false;}};
    
    inline bool LexFunctionDefined(){return function;}
    inline std::string getLexicalFunctionName(){return lexFunc->getName();}
	
	inline std::string getPDFFunctionName(){return pdfFunctionName;}
	
	inline bool isSimple(){
		if (transition_type != DURATION && func == NULL && !function){
			return true;
		}
		return false;
	}
	
	inline bool isComplex(){
		return !isSimple();
	}
    
    void print();
    std::string stringify();
    
private:    
	transType transition_type;  //0: standard  1:USER DISTRIBUTION  2:INTERNAL DISTRIBUTION	3:LEXICAL
    
    //Transition to State
    std::string stateName;  //What state we are transitioning to (Fill out when parsing)
    state* toState;    //pointer to the state (Filled during HMM finalization)
	
    /*--------------- STANDARD TRANSITIONS ----------------*/
	double log_trans; //Log of standard transition probability
    
    /*--------------- EXPLICIT DURATION DISTRIBUTION TABLES ----------------*/
    std::vector<double>* distribution;  //! Transition Length Distribution
    double extendedValue;
    
    tracebackIdentifier traceback_identifier;   //0:until different state	1:STATE_NAME	2:STATE_LABEL	3:STATE_GFF_TAG   4:START(INIT)  
	//if not defined it should traceback until same state ends.... default to zero
    
    std::string traceback_string;	//name,label,or gff tag to traceback to.
    
	
	/*--------------- Lexical Table ----------------*/
    bool function;
    emissionFuncParam* lexFunc;
    lexicalTable scoreTable;
    
	/*--------------- PDF Table ----------------*/
	pdfFunc* pdfFunction;
	std::string pdfFunctionName;
	track* pdfTrack; //Track to pass the PDFFunction
    size_t track_number;
	
    /*--------------- External Functions ----------------*/
    //Extern function Identifies function to be called at tagged transition;
    //Function has to be initialized before importing the model
    //Function parameters are <int> position, <string> sequence, <string> TracebackSequence
    //Function return value should be double in log space;
    
    transitionFuncParam* func;
    //track* externTrack;
    
    
    //Private Methods
    bool _parseStandard(std::string&,stringList&, valueType);
    bool _parseDuration(stringList&, stringList&, valueType);              
    bool _parseLexical(stringList&, stringList&, valueType, tracks&, StateFuncs*);
	bool _parsePDF(stringList&,stringList&,valueType,tracks&, StateFuncs*);
    bool _processTags(std::string&, tracks& , weights*, StateFuncs*);
};


}
#endif //TRANSITIONS_H//