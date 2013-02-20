//
//  stoch_table.h
//  StochHMM
//
//  Created by Paul Lott on 2/4/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#ifndef __StochHMM__stoch_table__
#define __StochHMM__stoch_table__

#include <iostream>
#include <sstream>
#include <string>
#include "stochMath.h"
#include "traceback_path.h"

namespace StochHMM {
	
	/*! \class stochTable
	 * Stochastic table stores stochastic values in a vector. Data structure is 
	 * implemented to reduce the memory that would be needed if we allocated a 2D
	 * table with all the transitions.   This lookup table is 1D and only stores
	 * those values which we add.   Each the value are stores
	 * (Current state, Previous State, Score as P(X)).  In addition, the indices
	 * of each state and position are recorded so that they can be referenced
	 * quickly.  So the memory size necessary to store the values for a sparsely
	 * connected model are greatly reduce.  However, the time to recreate the table
	 * is increased because it doesn't allocate the memory up front, but rather 
	 * pushes values as they are recorded.
	*/
	class stochTable{
	public:
		
		/*! \struct stoch_value;
		 * Structure used to store the traceback values (Current state, 
		 * Previous State, Traceback probability). Also includes the iterator of
		 * the previous state (which iterator is the previous state in relation
		 * to all states in the last position
		 */
		struct stoch_value{
			stoch_value(uint16_t id, uint16_t prev, float p): state_id(id), state_prev(prev), prob(p){}
			uint16_t state_id;
			uint16_t state_prev;
			uint16_t prev_cell;
			float prob;
		};
		
		
		stochTable(size_t);
		~stochTable();
		void push(size_t pos, size_t st, size_t st_to, float val);
		std::string stringify();
		void print();
		void finalize();
		size_t get_state_position(size_t pos,uint16_t);
		
		void traceback(traceback_path& path);
		
	private:
		size_t last_position;
		std::vector<stoch_value>* state_val;
		std::vector<size_t>* position;
	};

}

#endif /* defined(__StochHMM__stoch_table__) */
