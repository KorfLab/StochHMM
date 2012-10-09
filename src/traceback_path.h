///traceback_path.h

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

#ifndef TRACEBACK_PATH_H
#define TRACEBACK_PATH_H

#include <vector>
#include <string>
#include <iostream>
#include <math.h>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <stdint.h>
#include "options.h"
#include "hmm.h"
#include "stochMath.h"
namespace StochHMM{

    //! \struct gff_feature
    //!Each gff_feature represents a single GFF line
    struct gff_feature{
        std::string seqname;    //Column 1
        std::string source;     //Column 2
        std::string feature;    //Column 3
        size_t start;              //Column 4
        size_t end;                //Column 5
        char score;             //Column 6
        char strand;            //Column 7
        char frame;             //Column 8
        std::string attribute;  //Column 9
    };


    //! \class traceback_path
    //! Stores one traceback path information for a sequence
    class traceback_path{
    public:
        traceback_path(model*);
        
        void push_back(int);
        void clear();
        size_t size() const;
        int val(int);
        inline model* getModel() const {return hmm;};
        
        void fprint_path(std::ofstream&);
        
        void label(std::vector<std::string>&);
        void label(std::string&);
        
        void gff(std::vector<gff_feature>&,std::string&);
        void name(std::vector<std::string>&);
        void path(std::vector<int>&);
        
        void print_path();
        void print_label();
        void print_gff(std::string,double,int,int,double);
        void print_gff(std::string);
        
        inline double getScore(){
            return score;
        }
        
        inline void setScore(double scr){
            score = scr;
            return;
        };
        
        //double path_prob (const HMM&, sequence&);  //Need to rewrite function
        inline int operator[](int val) const {return trace_path[val];};
        bool operator== (const traceback_path&) const;
        bool operator<  (const traceback_path&) const;
        bool operator>  (const traceback_path&) const;
        bool operator<= (const traceback_path&) const;
        bool operator>= (const traceback_path&) const;
    private:
        model* hmm;
        std::vector<int> trace_path;
        double score;
    };

    //Rows are the Sequence Position
    //Columns are the States
    //! \def std::vector<std::vector<int> > heatTable;
    //! Table for heat map data generated from multiple tracebacks
    typedef std::vector<std::vector<int> > heatTable;

    //! \class multiTraceback
    //! Contains multiple tracebacks.  Will store them in sorted unique list (sorted in order of number of occurances);
    class multiTraceback{
    public:
        multiTraceback();
        ~multiTraceback();
        
        //Iteratorate through map;
        void begin();
        void end();
        void operator++();
        void operator--();
        void operator=(size_t);
        
        
        //Access the data at a point
        traceback_path path();
        int counts();
        traceback_path operator[](size_t);
        
        
        //Assign
        void assign(traceback_path&);
        
        //Finalize multiTraceback (Sort and setup Iterators);
        void finalize();
            
        inline void clear(){paths.clear();return;};
        inline size_t size(){return paths.size();};
        
        heatTable* heat();

    private:
        size_t vectorIterator;
        size_t maxSize;
        std::vector<std::map<traceback_path,int>::iterator> pathAccess;
        std::map<traceback_path,int> paths;
        heatTable* table;
    };
     
    bool sortTBVec(std::map<traceback_path,int>::iterator, std::map<traceback_path,int>::iterator);
}
#endif /*TRACEBACK_PATH_H*/



