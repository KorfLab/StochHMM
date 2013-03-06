//
//  stateInfo.h
//  StochHMM
//
//  Created by Paul Lott on 2/26/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#ifndef __StochHMM__stateInfo__
#define __StochHMM__stateInfo__

#include <iostream>
#include <map>
#include <vector>

namespace StochHMM{
	
	
	//!\struct stateInfo
	//!\brief Contains state information for quick reference between states information
	struct stateInfo{
		std::map<std::string,size_t> stateIterByName; //Map string iterators by name
		std::map<std::string,std::vector<size_t> > stateIterByGff; //Map states iterators by GFF tags
		std::map<std::string,std::vector<size_t> > stateIterByLabel; //Map states iterators by Label tags
	};
	
	
}


#endif /* defined(__StochHMM__stateInfo__) */
