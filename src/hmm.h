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
#include "stateInfo.h"
namespace StochHMM{
	
	
	
	/*! Hidden Markov Model Class
	 \class model class combines the States, and model information together in a single unit.
	 This includes the states(emissions, transitions), initial and ending states, track 
	 information(alphabet and ambiguous character definitions).
	 
	 Provides functions to import the model from a text file
	 
	 Model is used by trellis class to evaluates sequences.  
	 
	 
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
		
		
		
		//-----------State Information-------------/
		
		//!Get the number of states that are defined in the model
		inline size_t state_size(){return states.size();}
		
		//!Get the name of the state at index
		//! \param iter Index of state
		inline std::string& getStateName(size_t iter){
			if (iter>= states.size()){
				std::cerr << "Attempting to access State Name which is out of range\n";
				exit(2);
			}
			return states[iter]->getName();
		};
		
		//!Get the Label of the state at index
		//! \param iter Index of state
		inline std::string& getStateLabel(size_t iter){
			if (iter>= states.size()){
				std::cerr << "Attempting to access State Label which is out of range\n";
				exit(2);
			}
			return states[iter]->getLabel();
		}
		
		//!Get the GFF Tag of the state at index
		//! \param iter Index of state
		inline std::string& getStateGFF(size_t iter) {
			if (iter>= states.size()){
				std::cerr << "Attempting to access State GFF which is out of range\n";
				exit(2);
			}
			return states[iter]->getGFF();
		}
		
		//!Get pointer to the state at index
		//! \param iter Index of state
		//! \return ptr_state Pointer to state
		inline state*  getState(size_t iter){
			if (iter>= states.size()){
				return NULL;
			}
			return states[iter];
		}
		
		//!Get pointer to state by the name
		//! \param const std::string  Name of state
		//! \return ptr_state Pointer to state
		state* getState(const std::string&);
		
		
		//!Get state by using iterator value
		inline state* operator[](size_t iter){
			if (iter>= states.size()){
				return NULL;
			}
			return states[iter];
		}
		
		//!Get vector of states that state at index transitions to
		inline std::bitset<STATE_MAX>* getStateXTo(size_t iter){
			if (iter>= states.size()){
				return NULL;
			}
			return &(states[iter]->to);
		}
		
		//!Get vector of states that the initial state transitions to
		inline std::bitset<STATE_MAX>* getInitialTo(){return &(initial->to);}
		
		//!Get vector of states that transfer to the state at index
		inline std::bitset<STATE_MAX>* getStateXFrom(size_t iter){
			if (iter>= states.size()){
				return NULL;
			}
			return &(states[iter]->from);
		}
		
		//!Get list of states that transition to the ending state
		inline std::bitset<STATE_MAX>* getEndingFrom(){return &(ending->from);}
		
		inline stateInfo* getStateInfo(){return &info;}
		
		//--------- Initial and Ending States
		
		//!Get pointer to the initial state
		inline state*  getInitial(){return initial;}
		
		//!Get pointer to the ending state
		inline state*  getEnding(){return ending;}
		
		
		//--------- Scaling Factors
		//!Get Scaling factor defined in model by name
		//! \param std::string Name of Scaling or Weight defined in model
		//! \return ptr to weight 
		weight* getScalingFactor(std::string&);
		
		
		//--------- Attrib
		//!Get distance to value
		//!This is used when from selecting multiple models
		//!User can set an attribute value.   Then when evaluating
		//!an attribute of sequence they can see which model is closest
		//!and choose that model
		double getDistanceToAttrib(double);
		
		
		//---------- Track Information
		
		//!Get the number of tracks defined in the model
		inline size_t track_size(){return trcks.size();}
		
		//!Get pointer to track at the index
		//!\param iter Index of track to get
		//! \return if iter is within range then return pointer to track
		//! \return else return NULL
		inline track* getTrack(size_t iter){
			if (iter >= trcks.size()){
				return NULL;
			}
			return trcks[iter];
		}
		
		//!Get pointer to track based on Name associated with the track
		//! \param const std::string Name associated with Track
		//! \return if name is found returns pointer to name.
		//! \return if name is not found return NULL
		track* getTrack(const std::string&);
		
		//!Get index iterator of the track with a particular name
		//!\param txt Name of track to get index for
		//!\return size_t index of track with name
		inline size_t getTrackIter(const std::string& txt){return trcks.indexOf(txt);}
		
		//!Get pointer to the tracks of the model
		//!\return pointer to tracks defined in model
		inline tracks* getTracks(){return &trcks;}
		
		//!Check to see if model is a basic HMM
		//!\return false if model contains explicit duration transition, or user-defined functions for emission/transitions
		inline bool isBasic(){return basicModel;}
		
		
		
		
		//----------------  Printing or getting String representation of Model
		
		//! Print model by std::cout
		void print();
		
		//! Get text representation of the model
		//! \return std::string of model.
		std::string stringify();
		
		//void writeGraphViz(std::string);
		//! Write a simple GraphViz graph
		//! Formatting is very basic and function may disappear.
		//void writeGraphViz(std::string,bool);

		
		
		//MUTATORS
		//!Import and Parse the model from text file
		//!\param std::string Filename
		//!\param StateFuncs ptr  Pointer to StateFuncts, if no State Functions
		//! are (Univariate, Multivariate, Emission Functs, Transition Functions)
		//! described then you can use NULL
		//! \return true if import was successful
		bool import(std::string&,StateFuncs*);
		
		//!Import and Parse the model from text file
		//! \param std::string Filename
		//! \return true if import was successful
		bool import(std::string&);
		
		//!Import and Parse the model from text file
		//! \param std::string Filename
		//! \param StateFuncs ptr  Pointer to StateFunctions
		//! \param templates ptr Pointer to Templated State template
		//! \param weight ptr  Pointer to weighting factors
		//! \return true if import was successful
		bool import(std::string&, StateFuncs*, templates*, weights*);
		
		
		//!Import and parse the model from std::string
		//! \sa import(std::string&)
		bool importFromString(std::string&);
		
		//!Import and parse the model from std::string
		bool importFromString(std::string&,StateFuncs*);
		
		//!Import and parse the model from std::string
		bool importFromString(std::string&, StateFuncs*, templates*, weights*);
		
		//!Parse the model from std::string
		//!This is used by import functions to parse the model
		bool parse(const std::string&, StateFuncs*, templates*, weights*);
		
		//!Parse the model from std::string
		bool parse(std::string&,std::string&);
		
		
		//--------------  Set Model Data
		
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
		
		
		//----------------- Finalize and Check Final Model
		
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
		//!Flag set to tell whether the transitions bitsets have been set foreach
		//!state.  Model is also checked for correct order of states
		bool finalized;   
		
		
		//!Flag for whether model contains anything other than simple transitions and emissions
		//!If False then the model either contains additional function or emissions
		bool basicModel;
		
		std::string name;	//! Model Name
		std::string desc;	//! Model Description
		std::string date;   //! Model Creation Date
		std::string command; //! Model Creation Command
		std::string author;  //! Model Author
		float range[2];      //! Model Attrib Values
		bool attribTwo;      //! Two attrib Values
		
		tracks trcks; //! Tracks defined by model (Contains alphabet and ambiguous character definitions
		
		std::vector<state*> states; //!  All the states contained in the model

		std::map<std::string,state*> stateByName; //Ptr to state stored by State name;
		stateInfo info;
		
		
		state* initial; //!Initial state q0
		state* ending;	//!Ending state
		
		weights* scaling;  //! Weights or scaling fractors associated with the model
		
		//std::map<std::string,weight*> scaling;
		templates* templatedStates; //!Templated states
		
		std::vector<bool>* explicit_duration_states;	//! States that are explicit duration states
		
		std::vector<bool>* complex_transition_states;	//! States that have functions associated with transitions
		std::vector<bool>* complex_emission_states;		//! States that have functions associated with emissions
		
		bool _parseHeader(std::string&);	//! Function to parse header of the model from text file
		bool _parseTracks(std::string&);	//! Parse Tracks definitions from text file
		bool _parseAmbiguous(std::string&);	//! Parse Ambiguous definitions from text file
		bool _parseScaling(std::string&);	//! Parse Scaling definitions from text file
		bool _parseTemplates(std::string&);	//! Parse Templated States definitions from text file
		
		bool _parseStates(std::string&,StateFuncs*); //!Parse state from text file
		bool _splitStates(std::string&,stringList&); //!Split the state definitions into individual states from text file
		bool _getOrderedStateNames(stringList&,stringList&); //! Gets list of states names from model
		bool _processTemplateState(std::string&, stringList&); //! Adds templated states to using template
		
		std::string _stringifyHeader();	//!Converts Header information from model to string representation found in text file
		std::string _stringifyTracks(); //!Converts Tracks information from model to string representation found in text file
		std::string _stringifyAmbig();  //!Converts Ambiguous Character information from model to text string
		std::string _stringifyScaling();//!Converts Scaling definitions from model to text string
		std::string _stringifyStates(); //!Converts States definitions from model to text string
		
		
		void _addStateToFromTransition(state*); //!Processes each statea and defines definitions to state and from state for use
												//!in banding the trellis decoding functions
		
		void checkBasicModel();	//!Checks to see if the model has basic transitions and emissions(no addtl functions)
		void checkExplicitDurationStates();  //!Checks to see which states are explicit duration states
		void _checkTopology(state* st, std::vector<uint16_t>& visited); //!Checks to see that all states are connected and there
			
		
	};
	
	
	//----------------------------------------------------------------------------//
	//! models is a class to store multiple models.  This allows StochMM the ability, to
	//! load multiple models, then select the model that appropriate for the sequence
	//! based on a used-defined attribute.
	
	//! Stores multiple HMM models and contains the get functions for specific models
	//!
	//!
	//----------------------------------------------------------------------------//
	class models{
	public:
		//CONSTRUCTOR
		
		
		//ACCESSOR
		
		//! Get model located at index
		//! \param iter Index iterator for model
		//! \return pointer to model
		inline model* operator[](size_t iter){
			if (iter>hmms.size()-1){
				return NULL;
			}
			
			return hmms[iter];
		};
		
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
	};
	
	
	//!Print 2D vector to std::cout
	void print_vec (std::vector<std::vector<double> >&);
	
	//    void markov_length_distribution(model*);
	//
	//    void markov_generate_sequence(model*);
}
#endif /*HMM_H*/
