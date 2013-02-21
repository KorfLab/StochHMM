//
//  sequence.cpp
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

#include "sequence.h"

namespace StochHMM{
    
    
    //!Create a sequence datatype
    sequence::sequence(){
        realSeq=false;
		real = NULL;
        seq=NULL;
        mask=NULL;
        seqtrk  = NULL;
        length  = 0;
        max_mask=-1;
        external= NULL;
        attrib  = -INFINITY;
    }
    
    //!Create a sequence data type
    //! \param realTrack  True if the sequence is a list of real numbers
    sequence::sequence(bool realTrack):mask(NULL){
        if (realTrack){
            real=new(std::nothrow) std::vector<double>;
            
            if (real==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
            
            seq=NULL;
        }
        else{
            real=NULL;
            seq=new(std::nothrow) std::vector<uint8_t>;
            
            if (seq==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
        }
        realSeq = realTrack;
        seqtrk  = NULL;
        length  = 0;
        max_mask=-1;
        external= NULL;
        attrib  = -INFINITY;
    }
    
    //! \brief Create a sequence typ
    //! \param vec Vector of doubles to used as the real numbers for the sequence
    //! \param tr  Pointer to the track to be used
    sequence::sequence(std::vector<double>* vec ,track* tr):mask(NULL){
        realSeq = true;
        seqtrk  = tr;
        real    = vec;
        external= NULL;
        attrib  = -INFINITY;
        length  = vec->size();
        max_mask=-1;
    }
    
    //! Create a sequence type
    //! \param sq Character string that represent sequence
    //! \param tr Track to be used to digitize sequence
    sequence::sequence(char* sq, track* tr):mask(NULL){
        length  = 0;
        max_mask=-1;
        real    = NULL;
        seq     = new(std::nothrow) std::vector<uint8_t>;
        
        if (seq==NULL){
            std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
            exit(1);
        }
        
        external= NULL;
        attrib  = -INFINITY;
        realSeq = false;
        seqtrk  = tr;
        undigitized = sq;
        _digitize();
        length=seq->size();
    }
    
    
    //! Create a sequence type
    //! \param sq std::string string that represent sequence
    //! \param tr Track to be used to digitize sequence
    sequence::sequence(std::string& sq, track* tr):mask(NULL){
        length  = 0;
        max_mask=-1;
        real    = NULL;
        seq     = new(std::nothrow) std::vector<uint8_t>;
        
        if (seq==NULL){
            std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
            exit(1);
        }
        
        external= NULL;
        attrib  = -INFINITY;
        realSeq = false;
        seqtrk  = tr;
        undigitized = sq;
        _digitize();
        length=seq->size();
    }
    
    //!Destroy sequence type
    sequence::~sequence(){
        delete seq;
        delete real;
        delete mask;
        
        seq     = NULL;
        real    = NULL;
        mask    = NULL;
        seqtrk  = NULL;
        external= NULL;
    }
    
    
    
    //! Copy constructor for sequence 
    //!
    sequence::sequence(const sequence& rhs){
        realSeq = rhs.realSeq;
        header  = rhs.header;
        attrib  = rhs.attrib;
        length  = rhs.length;
        seqtrk  = rhs.seqtrk;
        external= rhs.external;  //Need copy constructor for this
        max_mask= rhs.max_mask;
        undigitized=rhs.undigitized;
        
        if (rhs.seq!=NULL){
            seq = new(std::nothrow) std::vector<uint8_t>(*rhs.seq);
            if (seq==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
        }
        else{
            seq=NULL;
        }
        
        if (rhs.real!=NULL){
            real = new(std::nothrow) std::vector<double>(*rhs.real);
            if (real==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
        }
        else{
            real=NULL;
        }
        
        if (rhs.mask!=NULL){
            mask = new(std::nothrow) std::vector<int>(*rhs.mask);
            if (mask==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
        }
        else{
            mask = NULL;
        }
    }
    
    sequence& sequence::operator= (const sequence& rhs){
		
		//Clean up if necessary 
		if (external != NULL){
			delete external;
			external = NULL;
		}
		
		if (seq!= NULL){
			delete seq;
			seq = NULL;
		}
		
		if (real!=NULL){
			delete real;
			seq = NULL;
		}
		
		if (mask!= NULL){
			delete mask;
			mask = NULL;
		}
		
		//Copy rhs over to this
        realSeq = rhs.realSeq;
        header  = rhs.header;
        attrib  = rhs.attrib;
        length  = rhs.length;
        seqtrk  = rhs.seqtrk;
        external= rhs.external;
        max_mask= rhs.max_mask;
        undigitized=rhs.undigitized;
        
        if (rhs.seq!=NULL){
            seq = new(std::nothrow) std::vector<uint8_t>(*rhs.seq);
            if (seq==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
        }
        else{
            seq=NULL;
        }
        
        if (rhs.real!=NULL){
            real = new(std::nothrow) std::vector<double>(*rhs.real);
            if (real==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
        }
        else{
            real=NULL;
        }
        
        if (rhs.mask!=NULL){
            mask = new(std::nothrow) std::vector<int>(*rhs.mask);
            if (mask==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
        }
        else{
            mask = NULL;
        }
        
        return *this;
    }
    
    
    //!Get digitized sequence value at a position
    //! \param position Position in the sequence to get the value for
    //! \return short integer value of symbol at positiion
    uint8_t sequence::seqValue (size_t position){
        if (!realSeq){
            if (seq!=NULL){
                return (*seq)[position];
            }
            else{
                std::cerr << "sequence has not been digitized. \n";
                exit(1);
            }
        }
        
        std::cerr << "Invalid sequence type queried in sequence class" <<std::endl;
        exit(1);
        
    };

    //!Get real value of sequence at a position
    //! \param position Position in the sequence to get the value
    //! \return positioin 
    double sequence::realValue(size_t position){
        if (realSeq){
            if (real!=NULL){
                return (*real)[position];
            }
            else{
                std::cerr << "Values have not been imported.\n";
                exit(1);
            }
                
        }
        std::cerr << "Invalid sequence type queried in sequence class" <<std::endl;
        exit(1);
    };
    
    //!Get std::string representation of the string
    //!If the string is a real track, then it will return a string of doubles
    //!If the string is a non-real track, then it will return a string of shorts, where the shorts are the digitized value of the sequence according to the track
    //! \return std::string String representation of the sequence
    std::string sequence::stringify(){
        std::string output;
		if (!header.empty()){
			 output+= header +"\n";
		}
        
        if (!seq && !realSeq){
            output+=undigitized;
        }
        
        if (realSeq){
            for(int i=0;i<length;i++){
                output+= double_to_string((*real)[i]) + " ";
            }
        }
        else{
            for(int i=0;i<length;i++){
                output+= int_to_string((int)(*seq)[i]) + " ";
            }
        }
        
        
        if (mask){
            output += "\n";
            for(int i=0;i<length;i++){
                output+= int_to_string((int)(*mask)[i]) + " ";
            }
        }
                
        output+="\n";
        return output;
    }
    
    //!Get the undigitized value of the string
    //!If the string is a real-track then it will return the same as stringify()
    //!If the string is a non-real track, it will return undigitized sequence
    std::string sequence::undigitize(){
        
        if (realSeq){
            return stringify();
        }
        
        if (!seq){  //If the sequence is not digitized yet.  Return the undigitized version 
            return undigitized;
        }
		
        std::string output;
		if (!header.empty()){
			output+= header +"\n";
		}
		
        if (seqtrk!=NULL){
            size_t alphaMax = seqtrk->getAlphaMax();
            
            for (int i=0;i<length;i++){
                output+=seqtrk->getAlpha((*seq)[i]);
                if (alphaMax!=1){
                    output+=" ";
                }
            }
        }
        else {
            std::cerr << "track pointer is not defined.   Can't undigitize sequence without valid track\n";
        }
        
        return output;
    }
    
    
    //TODO: Fix return value to return false if no sequence was returned or parsed
    //!Extract sequence from a fasta file
    //! \param file File stream 
    //! \param trk Track to use for digitizing sequence
    //! \return true if function was able to get a sequence from the file
    bool sequence::getFasta(std::ifstream& file, track* trk){
        
        seqtrk=trk;
        
        if (!file.good()){
            return false;
        }
        
        //Find next header mark 
        while(file.peek() != '>'){
            std::string temp;
            getline(file,temp,'\n');
            
            if (!file.good()){
                std::cerr << "Sequence doesn't contain a header \">\" "<< std::endl;
                return false;
            }
        }
        
        getline(file,header,'\n');
        bool success;
        //get sequence
        std::string line;
        while(getline(file,line,'\n')){
            undigitized+=line;
            
            char nl_peek=file.peek();  // see if we have new sequence header on the next line
            if (nl_peek=='>'){
                success = _digitize();
                break;
            }
            else if (nl_peek=='['){
                success = _digitize();
                //external=getExDef(file,nextSeq->size());
            }
            else if (nl_peek==EOF){
                success = _digitize();
                break;
            }
            else{
                continue;
            }
        }
        
        length=seq->size();
        return success;
    }
    
    bool sequence::getMaskedFasta(std::ifstream& file, track* trk){
        seqtrk=trk;
        
        if (!file.good()){
            return false;
        }
        
        //Find next header mark 
        while(file.peek() != '>'){
            std::string temp;
            getline(file,temp,'\n');
            
            if (!file.good()){
                std::cerr << "Sequence doesn't contain a header \">\" "<< std::endl;
                return false;
            }
        }
        
        bool success;
        std::string line, mask_string;
        
        getline(file,header,'\n');
        
        if (file.good()) {
            getline(file,undigitized,'\n');
            success = _digitize();
            length=seq->size();
        }
        if (file.good()) {
            getline(file,mask_string,'\n');
            
            stringList lst;
            //split string on space, comma, tab
            clear_whitespace(mask_string,"\n\r");
            lst.splitString(mask_string, " ,\t");
            //alloc mask vector
            
            if (lst.size() == length) {
                
                if (mask!=NULL){
                    delete mask;
                }
                
                mask=new(std::nothrow) std::vector<int>;
                if (mask==NULL){
                    std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                    exit(1);
                }
            }
            else {
                std::cerr << "Mask length not equal to Sequence length." << std::endl;
                return false;
            }
            //pass vector to lst.toVecInt
            lst.toVecInt(*mask);
            max_mask=*max_element(mask->begin(),mask->end());
        }
        else {
            success = false;
        }

        return success;
        
    }
    
    //!Digitize the sequence 
    bool sequence::_digitize(){
        
        if (seqtrk==NULL){
            std::cerr << "Can't digitize sequence without a valid track defined\n";
            return false;
        }
		
		
        
        stringList lst;
        clear_whitespace(undigitized,"\n");
        if (seqtrk->getAlphaMax()>1){
            lst.splitString(undigitized, " ,\t");
        }
        else{
            lst.fromAlpha(undigitized, 1);
        }
		
		if (seq == NULL){
			seq = new std::vector<uint8_t>(lst.size());
		}
		else{
			seq->assign(lst.size(),0);
		}
        
        for (size_t i=0;i<lst.size();i++){
            short symbl = seqtrk->symbolIndex(lst[i]);
            
            //Check ambiguous here
            if (!seqtrk->isAmbiguousSet() && symbl<0){
                std::cerr << "Can't digitize ambiguous character without ambiguous characters being allowed in the model." << std::endl;
                return false;
            }
            
            (*seq)[i] = symbl;
        }
        
        undigitized.clear();  //Once sequence is digitized we don't need the old seqeunce string
        
        return true;
    }
    

    //! Import one fastq entry from the file
    //!FastQ format:
    //!Line 1: Start with @
    //!Line 2: Sequence  , Can be multiple lines
    //!Line 3: Start with +
    //!Line 4: Quality Score , Can be multiple lines
    //! \param file File stream to file
    //! \param trk Track to used to digitize
    //! \return true if the sequence was successfully imported
    bool sequence::getFastq(std::ifstream& file, track* trk){
        seqtrk=trk;
            
        if (!file.good()){
            return false;
        }
        
        //Find first line that starts with "@" and Get header
        
        //Move down until the next line has a "@"
        while(file.peek() != '@' && file.good()){
            std::string temp;
            getline(file,temp,'\n');
        }
        
        //Get Header (One line)
        if (file.good()){
            getline(file,header,'\n');
        }
        else{
            return false;
        }
        
        
        std::string sequence="";
        //Get sequence (Multiple Lines)
        if (file.good()){
            while(getline(file,sequence,'\n')){
                undigitized+=sequence;

                char nl_peek=file.peek();  // see if we have new sequence header on the next line
                if (nl_peek=='+'){  //If next line begins with + then we are at a 
                    break;
                }
                else if (nl_peek==EOF){
                    break;
                }
                else{
                    continue;
                }
            }
        }
        else{
            return false;
        }
        
        _digitize();
        
        std::string quality_string;
        
        //Get Quality String (Multiple Lines)
        if (file.good()){
            while(getline(file,sequence,'\n')){
                quality_string+=sequence;
                
                char nl_peek=file.peek();  // see if we have new sequence header on the next line
                if (nl_peek=='@'){  //If next line begins with + then we are at a
                    if (quality_string.size() < undigitized.size()){
                        continue;
                    }
                    
                    break;
                }
                else if (nl_peek==EOF){
                    break;
                }
                else{
                    continue;
                }
            }
        }
        else{
            return false;
        }
        
        length=seq->size();
        return true;
    }
    
    
    //TODO: fix the return, so it returns false if sequence wasn't imported from file
    //! Import one Real number sequence from the file
    //! \param file File stream to file
    //! \param trk Track to used to digitize
    //! \return true if the sequence was successfully imported
    bool sequence::getReal(std::ifstream& file, track* trk){
        seqtrk=trk;
                
        //get header
        while(file.peek() != '@' && !file.eof()){
            std::string temp;
            getline(file,temp,'\n');
        }
		
		if (file.eof()){
			std::cerr << "No header found for sequence.  Header should start line with \"@\".\n";
			exit(2);
		}
        
        std::string sequence="";
        getline(file,sequence,'\n');
        header=sequence;
        
        //get sequence
        std::string line;
        stringList lst;
        while(getline(file,line,'\n')){
            lst.splitString(line,",\t ");
            std::vector<double> temp (lst.toVecDouble());
            for(size_t i=0;i<temp.size();i++){
                real->push_back(temp[i]);
            }
            
            char nl_peek=file.peek();  // see if we have new sequence header on the next line
            if (nl_peek=='>'){
                _digitize();
                break;
            }
            else if (nl_peek=='['){
                _digitize();
                //external=getExDef(file,nextSeq->size());
            }
            else if (nl_peek==EOF){
                _digitize();
                break;
            }
            else{
                continue;
            }
        }
        
        length=real->size();
        return true;
    }
    
    //! Return the mask at sequence position
    int sequence::getMask(size_t position) {
        
        if (mask!=NULL){
            if (position < getLength()) {
                return (*mask)[position];
            }
            else {
                std::cerr << "Position exceeds sequence length.\n";
                exit(1);
            }
        }
        else{
            std::cerr << "No Mask information.\n";
            exit(1);
        }
    }
    
    //!Set the sequence from a std::string
    //! \param sq Sequence to be used as sequence
    //! \param tr Track to be used to digitize sequence
    void sequence::setSeq (std::string& sq, track* tr){
        realSeq= false;
        undigitized = sq;
		
		if (tr!= NULL){
			seqtrk = tr;
			_digitize();
			length = seq->size();
		}
		else{
			length = sq.size();
		}
		
		return;
    }
    
    //!Set the sequence from a vector of doubles
    //! \param rl Vector of doubles to be used as real number sequence
    //! \param tr Track to be used to digitize sequence
    void sequence::setRealSeq(std::vector<double>* rl, track* tr){
        seqtrk = tr;
        realSeq = true;
        real=rl;
		
		length = rl->size();
		
		return;
    }
    
    
    //! Get the symbol (alphabet character or word) for a a given position of a alphanumerical sequence
    //! \param pos  Position within sequence
    //! \return string String of the symbol
    //! \sa track::getAlpha(short)
    std::string sequence::getSymbol(size_t pos) const{
        if (realSeq){
            std::cerr << "Can't get symbol of real values\n";
            return "";
        }
        
        if (seqtrk==NULL){
            std::cerr << "track is undefined for sequence\n";
            return "";
        }
        
        return seqtrk->getAlpha((*seq)[pos]);
    }
    
    
//    //! Get the index value of word at a given position and with a order of dependence)
//    //! Calulates the index for a character/word and returns an Index type.   If the word is non-ambiguous then it will return 1 value for row and column.  First pair is column index,   Second pair is row index.  If symbol is ambiguous, index contains all possible words that match.
//    //! \param[in] position Position within the sequence
//    //! \param[in] order Order of dependence
//    //! \param[out] word index pair of index values
//    //! \sa Index  
//    
//    void sequence::get_index(size_t position, int order, std::pair<Index,Index>& word_index){
//        int letter=seqValue(position);
//        size_t alphabet=seqtrk->getAlphaSize();
//                
//        Index& letter_index=word_index.first;
//        
//        if (letter<0){
//            letter_index.setAmbiguous(seqtrk->getAmbiguousSet(letter));
//        }
//        else{
//            letter_index+=letter;
//        }
//        
//        Index& y_subtotal = word_index.second;
//        
//        if (order!=0){
//            
//            for(int k=order;k>=1;k--){
//                int prev_letter=seqValue(position-k);
//                Index word;
//                if (prev_letter<0){
//                    word.setAmbiguous(seqtrk->getAmbiguousSet(prev_letter));
//                    y_subtotal+=word * POWER[k-1][alphabet-1];
//                }
//                else{
//                    y_subtotal+=prev_letter * POWER[k-1][alphabet-1];
//                }
//            }
//        }
//        
//        return;
//    }
    
    //! Reverse the sequence; If mask is defined, the mask will also be reversed
    //! \return true if reverse was successfully performed on sequence
    bool sequence::reverse(){
        if (realSeq){
            if (real!=NULL){
                std::reverse(real->begin(), real->end());
                return true;
            }
        }
        else if (seq!=NULL){
            std::reverse(seq->begin(), seq->end());
            
            if (mask!=NULL){
                std::reverse(mask->begin(), mask->end());
            }
            return true;
        }
        
        //Handle non-digitized sequence
        if (seqtrk!=NULL){
            size_t max_size = seqtrk->getAlphaMax();
            if (max_size ==1){
                if (undigitized.size()>0){
                    std::reverse(undigitized.begin(),undigitized.end());
                    return true;
                }
                else{
                    std::cerr << "No sequence is defined to reverse\n";
                }
            }
            else{
                std::cerr << "Reverse on undigitized sequence isn't defined for track alphabets that are more than one character\n";
            }
        }
        
        return false;
    }
    
    
    //! Complements the sequence.  
    //! \return true if complementation was successful
    bool sequence::complement() {
        if (realSeq){
            std::cerr<< "sequence::complement isn't defined for real valued sequences\n";
        }
        else if (seqtrk==NULL){
            std::cerr << "track is not defined.  Can't complement without defined complement in track\n";
        }
        else if (seq!=NULL){
            
            for (int i = 0; i < seq->size(); i++) {
                (*seq)[i] = seqtrk->getComplementIndex((*seq)[i]);
            }
            
            return true;
        }
        else {
            size_t undigitized_size=undigitized.size();
            if (undigitized_size>0){
                size_t max_size = seqtrk->getAlphaMax();
                if (max_size ==1){
                    for(size_t i=0;i<undigitized_size;i++){
                        std::string character = seqtrk->getComplementSymbol(undigitized[i]);
                        undigitized[i] = character[0];
                    }
                    return true;
                    
                }
                else{
                    std::cerr << "Complement on undigitized sequence isn't defined for track alphabets that are more than one character\n";
                }
                
            }
            else{
                std::cerr << "No sequence defined\n";
            }
            
        }
        return false;
    }
    
    
    //! Reverses and complements the sequence
    //! \return true if both reverse and complement were successful
    bool sequence::reverseComplement() {
        if (!reverse()){
            std::cerr << "Unable to perform reverseComplement on sequence because reverse failed\n";
            return false;
        };
        if (!complement()){
             std::cerr << "Unable to perform reverseComplement on sequence because complement failed\n";
            return false;
        };
        return true;
    }
    
    //! Digitize a sequences
    bool sequence::digitize(){
        if (realSeq){
            return true;
        }
        else if (undigitized.size() > 0){
            _digitize();
            return true;
        }
        else if (seq!=NULL){
            std::cerr << "Digitized sequence already exists\n";
            return true;
        }
        
        std::cerr << "No undigitized sequence exists to convert\n";
        return false;
    }
	
	
	void sequence::shuffle(){
		
		if (realSeq){
			std::random_shuffle(real->begin(), real->end());
		}
		else if (seq!=NULL){
			std::random_shuffle(seq->begin(), seq->end());
		}
		else{
			std::random_shuffle(undigitized.begin(), undigitized.end());
		}
		return;
	}
	
	
	
	sequence random_sequence(std::vector<double>& freq, size_t length, track* tr){
        sequence random_seq;
        
        if (tr==NULL){
            std::cerr << "Track is not defined" << std::endl;
            return random_seq;
        }
        
        size_t alphaSize=tr->getAlphaSize();
        size_t freqSize=freq.size();
		
        if (alphaSize!=freqSize){
            std::cerr << "Frequency distribution size and Alphabet size must be the same." << std::endl;
            return random_seq;
        }
        
        //Create CDF of frequency distribution
        std::vector<std::pair<double,std::string> > cdf;
        double sum = 0.0;
        for(size_t i=0;i<freqSize;++i){
            sum+=freq[i];
            std::pair<double,std::string> val (sum, tr->getAlpha(i));
            cdf.push_back(val);
        }
        
        //Generate random sequence
        std::string random_string;
        for(size_t j=0;j<length;++j){
            double val = ((double)rand()/((double)(RAND_MAX)+(double)(1))); //Generate random
            for (size_t m=0;m<freqSize;++m){ //Check to see which value is
                if (cdf[m].first>=val){
                    random_string+=cdf[m].second;
                    break;
                }
            }
        }
        
        random_seq.setSeq(random_string, tr);
        return random_seq;
    }
    
}
