//table.cpp
 
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

#include "table.h"
namespace StochHMM{


//constructor inline
Table::Table(int seq, int states){
	header=std::vector<std::string>(states);
    std::vector<double> initial(states,0);
    std::vector<std::vector<double> > next(seq,initial);
	data=next;
	x=states;
	y=seq;
}


void Table::print(){
	
	//print header;
    std::cout << "Position";
	for(int i=0;i<x;i++){
        std::cout << "\t" << header[i];
	}
    std::cout << std::endl;
	
	//For every position print the row of states
	for(int i=0;i<y;i++){
        std::cout << i+1;
		for(int j=0;j<x;j++){
            std::cout <<"\t" << data[i][j];
		}
        std::cout << std::endl;
	}
}


}