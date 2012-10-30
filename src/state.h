//state.h
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


#ifndef STATE_H
#define STATE_H

#include <string>
#include <vector>
#include "text.h"
#include "emm.h"
#include "transitions.h"
#include <stdint.h>
#include <stdlib.h>

#ifndef SIZE_MAX
#define SIZE_MAX ((size_t)-1)
#endif

namespace StochHMM{

    class transition;
    
class state{
public:
    state();
	//state(int);
    state(std::string&,stringList&,tracks&,weights*, StateFuncs*); //!Create state from string
    ~state();
    
    friend class model;
    
    //ACCESSOR
    
    //!Get the states integer index value
    //! \return integer
    inline size_t getIterator(){return stateIterator;};
    
    //!Get the name of the state
    //! \return std::string Name of state
    inline std::string& getName(){return name;};
    
    //!Get the GFF tag to use for the state
    //! \results std::string GFF tag
    inline std::string& getGFF(){return gff;};
    
    //!Get the Label used for the state
    //! \result std::string 
    inline std::string& getLabel(){return label;};
    
    //!Get all transition defined for the state
    //!\return std::vector<transitions*>*  Pointer to all transitions in state
    inline std::vector<transition*>* getTransitions(){return transi;};
    
    //!Get all states that this state has transitions to
    //! \return std::vector<state*>* Pointer to all state that are transitioned to
    inline std::vector<state*>* getTo(){return &to;};
    inline std::vector<state*>::iterator getToBegin(){return to.begin();};
    inline std::vector<state*>::iterator getToEnd(){return to.end();};
    
    
    
    //!Get all states that transition to this state
    //!\return std::vector<state*>* Pointer to all state that transition to this state
    inline std::vector<state*>* getFrom(){return &from;};
    inline std::vector<state*>::iterator getFromBegin(){return from.begin();};
    inline std::vector<state*>::iterator getFromEnd(){return from.end();};
    
    
    //TODO: Check that undefined return values are NULL
    //!Get transition at index
    //! \param iter Index to get transition for
    //! \return transition* pointer to the transition
    inline transition* getTrans(size_t iter){return (*transi)[iter];};
    
    //!Get the ending transition 
    //! \return transition* Pointer to ending transition
    //! \If not defined return value should be NULL
    inline transition* getEnding(){return endi;};
    
    //!Get the emission defined at index
    //!\param iter Index of emm to get
    //!\results emm* pointer to emission
    inline emm* getEmission(size_t iter){return emission[iter];};
    
    double get_emission(sequences&, size_t);  //get emission for given (position)
	double get_trans(sequences&,int,int);  //get transition  (position,from or too)
    double getEndTrans();

    void print();
    std::string stringify();
    
    //MUTATORS
    bool parse(std::string&,stringList&,tracks&,weights*,StateFuncs*);
    
    //!Add the transition to the state
    //!\param trans Pointer to transition to add to the state
    inline void addTransition(transition* trans){transi->push_back(trans);};
    
    //!Set the ending transition to a given transition
    //!\param trans Pointer to transition to be used as ending transition
    inline void setEndingTransition(transition* trans){endi=trans;};
    
    //!Add emission to the state
    //!\param em Pointer to the emission to be added
    inline void addEmission(emm* em){emission.push_back(em);};
    
    //!Set the name of the state
    //!\param txt Name of the state
    inline void setName(std::string& txt){name=txt;};
    
    //!Set the GFF Tag for the state
    inline void setGFF(std::string& txt){gff=txt;};
    
    //!Set the Label for the state
    inline void setLabel(std::string& txt){label=txt;};
    
    //!Add state that this state transitions to
    inline void addToState(state* st){to.push_back(st);};
    
    //!Add state that transitions to this state
    inline void addFromState(state* st){from.push_back(st);};
    
    //!Set the index value to be used for the state 
    inline void setIter(size_t val){stateIterator=val;};
    
    void checkLabels(std::set<std::string>& ,std::set<std::string>& ,std::set<std::string>& );
    
    void _finalizeTransitions(std::map<std::string,state*>& state_index);
	
private:
    std::string name;	/* State name */
    std::string gff ;	/* State features description */
    std::string label;	/* State feature path label */
    
    std::vector<transition*>* transi;
	transition* endi;
    
    // Track Emissions //
    std::vector<emm*> emission;
	
    
	//Linking State Information (These are assigned at model finalization)
	size_t stateIterator;  //index of state in HMM
    std::vector<state*> to;
    std::vector<state*> from;
    
    bool _parseHeader(std::string&);
    bool _parseTransition(std::string&,stringList&, tracks&, weights* , StateFuncs*);
    bool _parseEmission(std::string&,stringList&, tracks&, weights*, StateFuncs*);
};

}
#endif /*STATE_H*/