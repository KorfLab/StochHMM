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
#include "stochMath.h"
#include "traceback_path.h"

namespace StochHMM {

	
	class stochTable{
	public:
		
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
		void print();
		void finalize();
		uint16_t get_state_position(size_t pos,uint16_t);
		
		void traceback(traceback_path& path);
		
	private:
		size_t last_position;
		std::vector<stoch_value>* state_val;
		std::vector<size_t>* position;
	};

}

#endif /* defined(__StochHMM__stoch_table__) */
