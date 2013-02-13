//
//  Lexical.h
//  StochHMM
//
//  Created by Paul Lott on 4/2/12.
//  Copyright (c) 2012 University of California, Davis. All rights reserved.
//

#ifndef StochHMM_Lexical_h
#define StochHMM_Lexical_h
#include <string>
#include <vector>
#include <math.h>
#include <ctype.h>
#include <algorithm>
#include <stdint.h>
#include <stdlib.h>
#include "track.h"
#include "index.h"
#include "externalFuncs.h"
#include "weight.h"
#include "sequences.h"
#include "stochTypes.h"
//#include "simpleTable.h"

namespace StochHMM{
	
    //! \class lexicalTable
    //! \brief Lexical table stores the log2 probabilities for both emissions and lexical transitions
    //!
    class lexicalTable{
	private:
		
		
    public:
        
        lexicalTable();
        
        ~lexicalTable();
        
        double getValue(sequences&,size_t);
		
		//!Initialize the final emission table with ambiguous characters
		//Creates the log_emission simpleTable
		void initialize_emission_table();
		double get_reduced_order(sequences& seq, size_t position);
                
        std::vector<std::vector<double> >* getCountsTable();
        std::vector<std::vector<double> >* getProbabilityTable();
        std::vector<std::vector<double> >* getLogProbabilityTable();
        
        void createTable(int rows, int columns, int pseudocount, valueType typ);
        
        void addTrack(track*,int);
        void assignTable(std::vector<std::vector<double> >*, valueType);
        
        //!Set how the emission will deal with unknown alphabet
        //! \param type enum UnknownCharScoringType
        inline void setUnkScoreType(unknownCharScoringType type){unknownScoreType=type;};
        
        //!Set a given score to be returned for unknownCharScoringType
        inline void setUnkScore(double val){unknownDefinedScore=val;};
        
        //!Get pointer to track at index position of emission
        //!\param iter Index iterator of position
        //!\return track* Track in emission
        inline track* getTrack(size_t iter){return trcks[iter];};
        
        //!Get the number of tracks defined in emission
        //!\return size_t
        inline size_t trackSize(){return trcks.size();};
        
        //!Get Orders of lexical emission will use for all tracks
        //!\return std::vector<int>
        inline std::vector<uint8_t>& getOrder(){return order;};
        inline uint8_t getOrder(size_t i){return order[i];}
        
        //! Get Log(prob) emission table
        //! \return std::vector<std::vector<double> >
        inline std::vector<std::vector<double> >& getLogEmm(){return *logProb;}
        
        //! Get the alphabet sizes for all tracks used in emission
        //! \return std::vector<int>
        inline std::vector<uint8_t>& getAlphaSize(){return alphabets;}
        inline uint8_t getAlphaSize(size_t i){return alphabets[i];}
        inline size_t getNumberOfAlphabets(){return alphabets.size();}
        
        inline unknownCharScoringType getAmbScoringType(){return unknownScoreType;}
        inline double getAmbDefinedScore(){return unknownDefinedScore;}
        
        //!Increment counts
        inline void incrementCounts(size_t word_index, size_t char_index) { if (counts != NULL) (*counts)[word_index][char_index]++; }
        
        //!Increment counts by double
        inline void incrementCountsDouble(size_t word_index, size_t char_index, double val) { if (counts != NULL) (*counts)[word_index][char_index]+= val; }
        
        std::string stringify();
        
        void print();
        
    private:
		unknownCharScoringType unknownScoreType;  //! What type of score to use with unknown
        double unknownDefinedScore;  //!Undefined character score
		
		size_t number_of_tracks;
		std::vector<track*> trcks;  //Pointer to tracks of interest
        std::vector<uint8_t> alphabets;  //alphabet sizes for each emission
		std::vector<uint8_t> max_unambiguous;
        std::vector<uint8_t> order;  //Orders for each emission
        uint8_t max_order;
		
        size_t y_dim;
		size_t* x_subarray;
		size_t* y_subarray;
        std::vector<std::vector<double> >* prob;     //p(x)
        std::vector<std::vector<double> >* counts;   //counts
        std::vector<std::vector<double> >* logProb;  //log2(P(x))
		
		
		size_t array_size;
		size_t dimensions;
		std::vector<size_t> subarray_value;   //Values used to decompose index into sequenece AAA(A)B(B)
		std::vector<size_t> subarray_sequence;
		std::vector<size_t> subarray_position;

		std::vector<size_t> decompose_values; //Values used to compose index from sequences AAAB(AB)
		std::vector<size_t> decompose_sequence;
		
		std::vector<double>* log_emission;		
		std::vector<std::vector<double>* > low_order_emissions;
		std::vector<std::vector<std::pair<size_t,size_t>* > >low_order_info;
		
		void init_table_dimension_values();
		void init_array_dimension_values();
		size_t convertIndex(size_t,size_t);
		
		void decompose(size_t row, size_t column, std::vector<uint8_t>& letters);
		void decompose(size_t index, std::vector<uint8_t>& letters);
		
		void transferValues(std::vector<bool>& transferred);
		size_t calculateArrayIndex(std::vector<uint8_t>& kmer);
		void expand_ambiguous(std::vector<uint8_t>& letters, std::vector<double>& expanded);
		std::vector<std::vector<uint8_t> >* expand_ambiguous(std::vector<std::vector<uint8_t> >* words, size_t letter);
		size_t calculateIndexFromDecomposed(std::vector<uint8_t>& word);
		double getAmbiguousScore(std::vector<uint8_t>& letters);
    };
    
	
	    
    
}


#endif
