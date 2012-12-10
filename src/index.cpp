//
//  index.cpp

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

#include "index.h"
namespace StochHMM{

    //!Create an index
    Index::Index(){
        index=0;
        amb=NULL;
        ambiguous=false;
    }
    
    //!Create an Index by copy
    Index::Index(const Index& a){
        index=a.index;
        ambiguous=a.ambiguous;
        
        if (amb!=NULL){
            amb=new(std::nothrow) std::vector<size_t>(a.amb->size());
            
            if (amb==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
            
            copy( a.amb->begin(), a.amb->end(), amb->begin());
        }
    }
    
    //!Destroy Index
    Index::~Index(){
        delete amb;
        amb=NULL;
    }



    //!Assign Index by copy
    //!\param rhs
    Index& Index::operator= (const Index& rhs){
        index=rhs.index;
        ambiguous=rhs.ambiguous;
        
        if (amb!=NULL){
            amb=new(std::nothrow) std::vector<size_t>(rhs.amb->size());
            
            if (amb==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
            
            copy( rhs.amb->begin(), rhs.amb->end(), amb->begin());
        }
        return *this;
    }
    
    //! Add indexes together
    //! \param lhs
    //! \return Index with values added
    Index Index::operator+ (const Index& lhs){
        Index temp;
        temp.index=lhs.index + this->index;
        if (ambiguous){
            if (lhs.ambiguous){
                temp.ambiguous=true;
                temp.amb=new(std::nothrow) std::vector<size_t>;
                
                if (temp.amb==NULL){
                    std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                    exit(1);
                }
                
                addVectorCombinatorial(*temp.amb, *this->amb, *lhs.amb);
                
            }
            else{
                temp.ambiguous=true;
                temp.amb=new(std::nothrow) std::vector<size_t>(*this->amb);
                
                if (temp.amb==NULL){
                    std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                    exit(1);
                }
                
                addValueToVector(*temp.amb, lhs.index);
            }
        }
        else{
            if (lhs.ambiguous){
                temp.ambiguous=true;
                temp.amb=new(std::nothrow) std::vector<size_t>(*lhs.amb);
                
                if (temp.amb==NULL){
                    std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                    exit(1);
                }
                
                addValueToVector(*temp.amb, this->index);
            }
        }
        return temp;
    }
    
    //!Add and Assign
    
    void Index::operator+= (const Index& lhs){
        if (!ambiguous && !lhs.ambiguous){
            index+=lhs.index;
            return;
        }
        else if (ambiguous){
            if (lhs.ambiguous){
                std::vector<size_t>* temp = new(std::nothrow) std::vector<size_t>;
                
                if (temp==NULL){
                    std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                    exit(1);
                }
                
                addVectorCombinatorial(*temp, *(this->amb), *(lhs.amb));
                delete amb;
                amb=temp;
            }
            else{
                addValueToVector(*this->amb, lhs.index);
            }
        }
        else{
            if (lhs.ambiguous){
                amb=new(std::nothrow) std::vector<size_t>(*lhs.amb);
				ambiguous = true;
                
                if (amb==NULL){
                    std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                    exit(1);
                }
                
                addValueToVector(*amb, index);
            }
        }
        return;
    }

    //!Add integer index to Index
    //! \param lhs Integer index to add to Index
    void Index::operator+= (const size_t lhs){
        if (ambiguous){
            addValueToVector(*amb, lhs);
        }
        else{
            index+=lhs;
        }
        return;
    }

    //!Multiply Indices by integer
    //! \param lhs Integer value
    void Index::operator*= (const size_t lhs){
        if (ambiguous){
            multiplyValueToVector(*amb, lhs);
        }
        else{
            index*=lhs;
        }
        return;
    }


    //!Multiply Indices by integer
    Index Index::operator* (const size_t lhs){
        Index temp;
        if (ambiguous){
            multiplyValueToVector(*amb, lhs);
        }
        else{
            index*=lhs;
        }
        return temp;
    }
    
    //!Get the number of indices in Index
    size_t Index::size(){
        if (ambiguous){
            return amb->size();
        }
        else{
            return 1;
        }
    }
    
    
    //! Get the index at position in Index
    //!\param iter Integer iterator of position
    //!\return int Integer value of index
    size_t Index::operator[](const size_t iter){
        if (!ambiguous && iter==0){
            return index;
        }
        else if (ambiguous){
            return (*amb)[iter];
        }
        else{
            std::cerr << "iterator: " << iter << " called on unambiguous value.   Can only access with Iterator = 0" << std::endl;
            exit(5);
        }
    }

        
    /*! \fn void setAmbiguous(const std::vector<int>& vec)
     \brief Sets the index class to ambiguous
     \param vec vector of ints that correspond to multiple indices created by ambiguous characters
    */
    void Index::setAmbiguous(const std::vector<size_t>& vec){
        amb=new(std::nothrow) std::vector<size_t>(vec);
        
        if (amb==NULL){
            std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
            exit(1);
        }
        
        ambiguous = true;
        return;
    }

}
