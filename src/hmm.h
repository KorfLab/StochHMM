//hmm.h
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
// TODO:  add survival function distribution option
// TODO:  Test internal distribution output options

#ifndef HMM_H
#define HMM_H


#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <map>
#include <list>
#include <set>
#include <stdlib.h>
#include <bitset>

#include "state.h"
#include "track.h"
#include "text.h"
#include "stochTypes.h"
#include "weight.h"
//#include "transInfoParse.h"   // Transitions Information Parsing
#include "modelTemplate.h"
namespace StochHMM{
	
	
	
	/*! Hidden Markov Model Class
	 Model class combines the States, and model information together in a single unit
	 Model is used by trellis classes to evaluates sequences
	 
	 
	 */
	class model {
	public:
		model();
		
		//    model(std::string&,StateFuncs*); //! Construct model from model file
		//    model(std::string&,std::string&,StateFuncs*); //!Construct model from model file and template file
		//    model(std::string&,StateFuncs*,templates*,weights*);
		
		//ACCESSOR FUNCTIONS
		
		//Model Information
		
		//!Get the Name of the Model
		inline std::string& getName(){return name;} // Get name of model
		
		//! Get Description of model
		inline std::string& getDescription(){return desc;}
		
		//! Get Date of model
		inline std::string& getDate(){return date;}
		
		//!Get Creation Command of model
		inline std::string& getCommand(){return command;}
		
		//! Get Author of model
		inline std::string& getAuthor(){return author;}
		//inline float getLowerRange(){return range[0];};
		//inline float getUpperRange(){return range[1];};
		
		//State Information
		
		//!Get the number of states that are defined in the model
		inline size_t state_size(){return states.size();}
		
		//!Get the name of the state at index
		//! \param iter Index of state
		inline std::string& getStateName(size_t iter){return states[iter]->getName();};
		
		//!Get the Label of the state at index
		inline std::string& getStateLabel(size_t iter){return states[iter]->getLabel();}
		
		//!Get the GFF Tag of the state at index
		inline std::string& getStateGFF(size_t iter) {return states[iter]->getGFF();}
		
		//!Get pointer to the state at index
		inline state*  getState(size_t iter){return states[iter];}
		
		state* getState(const std::string&);
		
		inline state* operator[](size_t iter){return states[iter];}
		
		//!Get vector of states that state at index transitions to
		inline std::bitset<STATE_MAX>* getStateXTo(size_t iter){return &(states[iter]->to);}
		
		//!Get vector of states that the initial state transitions to
		inline std::bitset<STATE_MAX>* getInitialTo(){return &(initial->to);}
		
		//!Get vector of states that transfer to the state at index
		inline std::bitset<STATE_MAX>* getStateXFrom(size_t iter){return &(states[iter]->from);}
		
		//!Get list of states that transition to the ending state
		inline std::bitset<STATE_MAX>* getEndingFrom(){return &(ending->from);}
		
		//q0 transitions
		
		//!Get pointer to the initial state
		inline state*  getInitial(){return initial;}
		
		//!Get pointer to the ending state
		inline state*  getEnding(){return ending;}
		
		
		//Scaling Factors
		weight* getScalingFactor(std::string&);
		
		//Attrib
		double getDistanceToAttrib(double);
		
		
		//Track Information
		
		//!Get the number of tracks defined in the model
		inline size_t track_size(){return trcks.size();}
		
		//!Get pointer to track at the index
		//!\param iter Index of track to get
		inline track* getTrack(size_t iter){return trcks[iter];}
		
		track* getTrack(const std::string&);
		
		//!Get index iterator of the track with a particular name
		//!\param txt Name of track to get index for
		inline size_t getTrackIter(const std::string& txt){return trcks.indexOf(txt);}
		
		//!Get pointer to the tracks of the model
		inline tracks* getTracks(){return &trcks;}
		
		inline bool isBasic(){return basicModel;}
		
		//!  Print model to stdout
		void print();
		
		
		//void writeGraphViz(std::string);
		void writeGraphViz(std::string,bool);
		
		//! Get text representation of the model
		std::string stringify();
		
		
		//MUTATORS
		bool import(std::string&,StateFuncs*); //! Parse Model from File
		bool import(std::string&);
		bool import(std::string&, StateFuncs*, templates*, weights*);
		
		bool parse(const std::string&, StateFuncs*, templates*, weights*);
		bool parse(std::string&,std::string&);
		
		///@{
		//!Mutator for building model internally also called by import
		
		//!Set the name of the model
		inline void setName(std::string& txt){name=txt;};
		
		//! Set the model desc
		inline void setDesc(std::string& txt){desc=txt;};
		
		//! Set the model creation date
		inline void setDate(std::string& txt){date=txt;};
		
		//! Set the model creation command
		inline void setCommand(std::string& txt){command=txt;};
		
		//! Set the model author
		inline void setAuthor(std::string& txt){author=txt;};
		
		//! Set the attribute of the model
		//! If model is going to be chosen from different values, then you can assign value
		inline void setNumericalAttrib(float value){range[0]=value;attribTwo=false;};
		
		//! Set the upper range of the attribute
		inline void setUpperRange(float& value){range[1]=value;attribTwo=true;};
		
		//! Set the lower range of the attribute
		inline void setLowerRange(float& value){range[0]=value;attribTwo=true;};
		
		//! Add track to the model
		inline void addTrack(track* trk){trcks.push_back(trk);};
		
		void addState(state*);
		
		//! Set the initial state pointer
		//! \param st pointer to initial state
		inline void setInit(state* st){initial=st;};
		
		//! Set teh ending state pointer
		inline void setEnd(state* st){ending=st;};
		//inline void addWeight(std::string& txt,weight* wt){scaling[txt]=wt;};
		
		//!Finalizes model references from and to states
		//!Each model must be finalized before being used to decode
		//!Check the Functions and Labels of the States
		void finalize();
		
		//!Check model topology
		//!Iterates through all states to check to see if there are any:
		//! 1. Orphaned States
		//! 2. Dead end States
		//! 3. Uncompleted States
		bool checkTopology();
		
		//Get a vector<bool> of states that are explicit duration states
		inline std::vector<bool>* get_explicit(){return explicit_duration_states;}
		
		bool hasComplexEmission(){
			if (complex_emission_states){
				return true;
			}
			return false;
		}
		
	private:
		bool finalized;
		bool basicModel;
		
		std::string name;	//! Model Name
		std::string desc;	//! Model Description
		std::string date;   //! Model Creation Date
		std::string command; //!Model Creation Command
		std::string author;  //!Model Author
		float range[2];      //!Model Attrib Values
		bool attribTwo;      //!Two attrib Values
		
		tracks trcks; //tracks...
		
		std::map<std::string,state*> statesByName;
		
		std::vector<state*> states;
		
		state* initial;
		state* ending;
		
		weights* scaling;  //Change to weights
		
		//std::map<std::string,weight*> scaling;
		templates* templatedStates;
		
		std::vector<bool>* explicit_duration_states;
		
		std::vector<bool>* complex_transition_states;
		std::vector<bool>* complex_emission_states;
		
		
		bool _parseHeader(std::string&);
		bool _parseTracks(std::string&);
		bool _parseAmbiguous(std::string&);
		bool _parseScaling(std::string&);
		bool _parseTemplates(std::string&);
		
		bool _parseStates(std::string&,StateFuncs*);
		bool _splitStates(std::string&,stringList&);
		bool _getOrderedStateNames(stringList&,stringList&);
		bool _processTemplateState(std::string&, stringList&);
		
		std::string _stringifyHeader();
		std::string _stringifyTracks();
		std::string _stringifyAmbig();
		std::string _stringifyScaling();
		std::string _stringifyStates();
		void _addStateToFromTransition(state*);
		
		void checkBasicModel();
		void checkExplicitDurationStates();
		void _checkTopology(state* st, std::vector<uint16_t>& visited);
		
	};
	
	
	//----------------------------------------------------------------------------//
	// Description:   multimodel class
	// Stores multiple HMM models and contains the get functions for specific models
	//
	//
	//----------------------------------------------------------------------------//
	class models{
	public:
		//CONSTRUCTOR
		
		
		//ACCESSOR
		
		//TODO: Fix so pointer is null if out of bound
		//! Get model located at index
		//! \param iter Index iterator for model
		//! \return pointer to model
		inline model* operator[](size_t iter){return hmms[iter];};
		
		//!Get the number of model
		//! \return size_t
		inline size_t size(){return hmms.size();};
		
		//FIXME: implement getModelByAttrib
		//model* getModelByAttrib(float);
		
		model* getModel(size_t);
		
		//MUTATOR
		void importModels(std::string&,StateFuncs*);
		void addModel(model*);
		
	private:
		std::vector<model*> hmms;
		weights* scaling;
		templates* modelTemplates;
		//int numberModels;
		//model* getGCModel(float);
	};
	
	
	
	void print_vec (std::vector<std::vector<double> >&);
	
	
}
#endif /*HMM_H*/
