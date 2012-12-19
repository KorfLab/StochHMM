//traceback_path.cpp
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

#include "traceback_path.h"

namespace StochHMM{


    int WID=80;
    //!Create a traceback_path
    //!\param modl Pointer to model file 
    traceback_path::traceback_path(model* modl){
        hmm=modl;
    }

    //!Pushes a state index onto the end of the path
    //!\param state Index of state to add
    void traceback_path::push_back(int state){
        trace_path.push_back(state);
    }
    
    //!Returns the size (ie. length) of the traceback_path
    size_t traceback_path::size()const {
        return trace_path.size();
    }
    
    //!Clears all traceback path information
    void traceback_path::clear(){
        trace_path.clear();
    }

    //!Returns the state index at a given position (it) within the sequence
    int traceback_path::val(int it){
        return trace_path[it];
    }
    
    
    //TODO: change assignment to the lhs
    //! Get the path to std::vector<int>
    //! \param [out] pth std::vector<int> that represents path
    void traceback_path::path(std::vector<int>& pth){
        pth=trace_path;
        return ;
    }

    //TODO: change assignment to lhs
    //! Get the label of the traceback_path and assigns to vector<string> ref
    void traceback_path::label(std::vector<std::string>& pth){
        
        for(size_t k=trace_path.size()-1; k!=SIZE_MAX; k--){
            state* st = hmm->getState(trace_path[k]);
            pth.push_back(st->getLabel());
        }
        return;
    }
    
    //!Get string of path label traceback
    //! \param[out] pth std::string
    void traceback_path::label(std::string& pth){
        
        if (pth.size()>0){
            pth.clear();
        }
        
        
        for(size_t k=trace_path.size()-1; k!=SIZE_MAX; k--){
            state* st = hmm->getState(trace_path[k]);
            pth+=st->getLabel();
        }
        return;
    }
    
    //!Get names of traceback path
    //!\param [out] pth vector<string>
    void traceback_path::name(std::vector<std::string>& pth){
        for(size_t k = trace_path.size()-1; k != SIZE_MAX; k--){
            state* st = hmm->getState(trace_path[k]);
            pth.push_back(st->getName());
        }
        return;
    }
    
    
    //!Get GFF output of traceback path
    //!\param [out] pth
    void traceback_path::gff(std::vector<gff_feature>& pth,std::string& sequenceName){
        std::string current_label="";
        long long start=0;
        size_t path_size=size();
        
        for(size_t k = path_size-1;k != SIZE_MAX; k--){
            state* st = hmm->getState(trace_path[k]);
            std::string new_label=st->getGFF();
            if (new_label.compare("")==0){
                if (start>0){
                    gff_feature ln;
                    ln.seqname=sequenceName;
                    ln.source=hmm->getName();
                    ln.feature=current_label;
                    ln.start=start;
                    ln.end=path_size-(k+1);
                    ln.score='.';
                    ln.strand='+';
                    ln.frame='.';
                    ln.attribute="";
                    
                    pth.push_back(ln);
                    
                    start=0;
                    current_label=new_label;
                }
                else {
                    continue;
                }
            }
            else {
                if(k==0){
                    gff_feature ln;
                    ln.seqname=sequenceName;
                    ln.source=hmm->getName();
                    ln.feature=current_label;
                    ln.start=start;
                    ln.end=path_size-(k+1);
                    ln.score='.';
                    ln.strand='+';
                    ln.frame='.';
                    ln.attribute="";
                    
                    pth.push_back(ln);
                }
                else if (start==0){
                    start=path_size-k;
                    current_label=new_label;
                }
                else if (new_label.compare(current_label)==0){
                    continue;
                }
                else {
                    gff_feature ln;
                    ln.seqname=sequenceName;
                    ln.source=hmm->getName();
                    ln.feature=current_label;
                    ln.start=start;
                    ln.end=path_size-(k+1);
                    ln.score='.';
                    ln.strand='+';
                    ln.frame='.';
                    ln.attribute="";
                    
                    pth.push_back(ln);
                    
                    start=path_size-k;
                    current_label=new_label;
                }
            }
            
        }
        
        return;
    }
    
    
    //!Print the path to stdout
    void traceback_path::print_path() const{
        int line=0;
        for(size_t k = this->size()-1; k != SIZE_MAX; k--){
            std::cout << trace_path[k]<< " ";
            line++;
        }
        std::cout << std::endl << std::endl;
    }

    //!Print the path to file stream
    void traceback_path::fprint_path(std::ofstream &file){
        int line=0;
        for(size_t k=this->size()-1;k != SIZE_MAX; k--){
            file << trace_path[k]<< " ";
            line++;
        }
        file << std::endl;
    }

    //!Check to see if paths are the same
    bool traceback_path::operator== (const traceback_path &rhs) const{
        if (rhs.trace_path==trace_path){
            return true;
        }
        else {
            return false;
        }
    }

    //!Comparison operators for path
    bool traceback_path::operator<  (const traceback_path &rhs ) const{
        if (trace_path<rhs.trace_path){
            return true;
        }
        else {
            return false;
        }
    }

    //!Comparison operators for path
    bool traceback_path::operator>  (const traceback_path &rhs) const{
        if (trace_path>rhs.trace_path){
            return true;
        }
        else {
            return false;
        }
    }

    //!Comparison operators for path
    bool traceback_path::operator<=  (const traceback_path &rhs) const{
        if (trace_path<=rhs.trace_path){
            return true;
        }
        else {
            return false;
        }
    }

    //!Comparison operators for path
    bool traceback_path::operator>=  (const traceback_path &rhs) const{
        if (trace_path>=rhs.trace_path){
            return true;
        }
        else {
            return false;
        }
    }


    //!Print traceback_path labels to stdout
    void traceback_path::print_label() const {
        int line=0;
        for(size_t k = trace_path.size()-1;k != SIZE_MAX;k--){
//            if(line==WID && WID>0){
//                std::cout<< std::endl;
//                line=0;
//            }
            state* st = hmm->getState(trace_path[k]);
            std::cout << st->getLabel() << " ";
            line++;
        }
        std::cout << std::endl << std::endl;
        
    }

    //!Outputs the gff formatted output for the traceback to stdout
    void traceback_path::print_gff(std::string sequence_name, double score, int ranking, int times, double posterior) const {
        std::string current_label="";
        long long start=0;
        size_t path_size=this->size();
        
        for(size_t k=path_size-1;k != SIZE_MAX;k--){
            state* st = hmm->getState(trace_path[k]);
            std::string new_label=st->getGFF();
            if (new_label.compare("")==0){
                if (start>0){
                    std::cout << sequence_name << "\tStochHMM\t" << current_label <<"\t"<< start << "\t" << path_size-(k+1) << "\t" << score << "\t+\t.\tRank:"<<ranking<<",Counts:" << times<<",Posterior:"<<posterior<<std::endl;
                    start=0;
                    current_label=new_label;
                }
                else {
                    continue;
                }
            }
            else {
                if(k==0){
                    std::cout << sequence_name << "\tStochHMM\t" << current_label <<"\t"<< start << "\t" << path_size << "\t" << score << "\t+\t.\tRank:"<<ranking<<",Counts:" << times<<",Posterior:"<<posterior<<std::endl;
                }
                else if (start==0){
                    start=path_size-k;
                    current_label=new_label;
                }
                else if (new_label.compare(current_label)==0){
                    continue;
                }
                else {
                    std::cout << sequence_name << "\tStochHMM\t" << current_label <<"\t"<< start << "\t" << path_size-(k+1) << "\t" << score << "\t+\t.\tRank:"<<ranking<<",Counts:" << times<<",Posterior:"<<posterior<<std::endl;
                    start=path_size-k;
                    current_label=new_label;
                }
            }
            
        }
    }	


    //!outputs the gff formatted output for the traceback
    void traceback_path::print_gff(std::string sequence_name) const {
        std::string current_label="";
        long long start=0;
        size_t path_size=size();
        //cout << trace_path.size() << endl;
        
        for(size_t k = path_size-1; k != SIZE_MAX; k--){
            state* st = hmm->getState(trace_path[k]);
            std::string new_label=st->getGFF();
            if (new_label.compare("")==0){
                if (start>0){
                    std::cout << sequence_name << "\tStochHMM\t" << current_label <<"\t"<< start << "\t" << path_size-(k+1) << "\t.\t+"<<std::endl;
                    start=0;
                    current_label=new_label;
                }
                else {
                    continue;
                }
            }
            else {
                if(k==0){
                    std::cout << sequence_name << "\tStochHMM\t" << current_label <<"\t"<< start << "\t" << path_size << "\t.\t+"<<std::endl;
                }
                else if (start==0){
                    start=path_size-k;
                    current_label=new_label;
                }
                else if (new_label.compare(current_label)==0){
                    continue;
                }
                else {
                    std::cout << sequence_name << "\tStochHMM\t" << current_label <<"\t"<< start << "\t" << path_size-(k+1) << "\t.\t+"<<std::endl;
                    start=path_size-k;
                    current_label=new_label;
                }
            }
            
        }
        
        std::cout << std::endl<<std::endl;
    }


    /* Need to re-write to handle new HMM Types
    //fix to handle higher order model 
    double traceback_path::path_prob (const HMM &model){
        int size= trace_path.size();
        
        int seq_size= model.seq_size();
        int alpha_size=model.alpha.size();
        
        
        if (size!=seq_size){
            cerr << "Sequence and Path different sizes\n";
            exit(1);
        }
        
        vector<vector<double> > log_emm=convert_order(model, trace_path[size-1], 0);
        double prob=model.initial.get_trans(trace_path[size-1]) + log_emm[0][seq.val(0)];
        
        
        for (unsigned int i=1;i<size;i++){
            int index=0;
            
            if (model.states[trace_path[size-i-1]].order>i){
                vector<vector<double> > emmiss=convert_order(model, trace_path[size-i-1], i);
                for(int n=i;n>=1;n--){
                    index+=POWER[n-1][alpha_size-1]*seq.val(i-n);
                }
                
                //cout <<i<<"\t"<<prob <<endl;
                prob+=model.states[trace_path[size-i]].get_trans(trace_path[size-i-1])+emmiss[index][seq.val(i)];
            }
            
            
            else{
                
                for(int n=model.states[trace_path[size-i-1]].order;n>=1;n--){
                    index+=POWER[n-1][alpha_size-1]*seq.val(i-n);
                    //cout <<"Index\t"<< index<<"\t"<< k << endl;
                }
                //cout <<i<<"\t"<<prob <<endl;
                prob+=model.states[trace_path[size-i]].get_trans(trace_path[size-i-1])+model.states[trace_path[size-i-1]].log_emm[index][seq.val(i)];
            }
        }
        
        return prob;	
    }
     */


    //!Create multiTraceback()
    multiTraceback::multiTraceback(){
        maxSize=0;
        vectorIterator=0;
        table=NULL;
    }

    //!Destroy multiTraceback
    multiTraceback::~multiTraceback(){
        if(table!=NULL){
            delete table;
        }
    }

    //!Set position in multiTraceback to the beginning
    void multiTraceback::begin(){
        vectorIterator=0;
        return;
    }
    
    //!Set position in multiTraceback to the ending
    void multiTraceback::end(){
        vectorIterator=paths.size();
        return;
    }
    
    //!Increment the iterator to next position
    void multiTraceback::operator++(){
        if (vectorIterator<maxSize){
            vectorIterator++;
        }
        return;
    }

    //!Decrement the iterator to previous position
    void multiTraceback::operator--(){
        if (vectorIterator>0){
            vectorIterator--;
        }
        return;
    }
    
    //!Set iterator to index val
    //!\param val Index value to set
    void multiTraceback::operator=(size_t val){
        if (val<=maxSize){
            vectorIterator=val;
        }
        return;
    }
    
    
    //!Get traceback_path at index position
    //! \param val Index position
    traceback_path multiTraceback::operator[](size_t val){
        return (*pathAccess[val]).first;
    }
    
    
    //!Get traceback_path at currently set index in multiTraceback
    traceback_path multiTraceback::path(){
        return (*pathAccess[vectorIterator]).first;
    }

    //!Get the number times that traceback_path was recorded in multiple traceback
    int multiTraceback::counts(){
        return (*pathAccess[vectorIterator]).second;
    }

    
    //!Add traceback_path to multiTraceback
    //!\param path Traceback path to add
    void multiTraceback::assign(traceback_path& path){
        paths[path]++;
        return;
    }


    //!Sorts the multiTraceback by the number of time a particular tracback path occurred
    void multiTraceback::finalize(){
        maxSize=paths.size();
        vectorIterator=0;
        
        std::map<traceback_path,int>::iterator pathsIterator;
        for(pathsIterator=paths.begin();pathsIterator!=paths.end();pathsIterator++){
            pathAccess.push_back(pathsIterator);
        }
        
        sort(pathAccess.begin(),pathAccess.end(),sortTBVec);
        return;
    }
    
    //!Generate a hit table from a multiple traceback paths
    //! Hit table is 2D table describing how many times a state was called at a particular position in the sequence
    heatTable* multiTraceback::get_hit_table(){
        if (table!=NULL){
            delete table;
        }
        
        //Over the lenght of the sequence
        model* hmm = ((*pathAccess[0]).first).getModel();
        size_t sequenceSize=((*pathAccess[0]).first).size();
        size_t stateSize=hmm->state_size();
        
        
        std::vector<int> states(stateSize,0);
        table = new heatTable(sequenceSize,states);
        
        std::map<traceback_path,int>::iterator it;
        
        for( it =paths.begin(); it!=paths.end();it++){
            int count = (*it).second;
            for(size_t position=0;position<sequenceSize;position++){
                int tbState=(*it).first[position];
                (*table)[position][tbState]+=count;
            }
        }
        return table;
    }
    
    
    void multiTraceback::print_hits(){
        if (table==NULL){
            get_hit_table();
        }
        
        std::string header_row = "Position";
        model* hmm = ((*pathAccess[0]).first).getModel();
        for (size_t state_iter =0; state_iter<hmm->state_size(); state_iter++){
            header_row+=",";
            header_row+=hmm->getStateName(state_iter);
        }

        std::cout << header_row << std::endl;
        
        for(size_t position = 0; position < table->size(); position++){
            std::string line = join((*table)[position], ',');
            std::cout << position << "," << line << std::endl;
        }
        
        return;
    }
    
    
    void multiTraceback::print_path(){
        this->finalize();
        for(size_t iter=0; iter<this->size(); iter++){
            std::cout << "Traceback occurred:\t " << (*pathAccess[iter]).second << std::endl;
            (*pathAccess[iter]).first.print_path();
            std::cout << std::endl;
        }
        return;
    }
    
    void multiTraceback::print_label(){
        this->finalize();
        for(size_t iter=0; iter<this->size(); iter++){
            std::cout << "Traceback occurred:\t " << (*pathAccess[iter]).second << std::endl;
            (*pathAccess[iter]).first.print_label();
            std::cout << std::endl;
        }
        return;
    }
    
    void multiTraceback::print_gff(std::string& header){
        this->finalize();
        for(size_t iter=0; iter<this->size(); iter++){
            std::cout << "Traceback occurred:\t " << (*pathAccess[iter]).second << std::endl;
            (*pathAccess[iter]).first.print_gff(header);
            std::cout << std::endl;
        }
        return;
    }
    
    
    bool sortTBVec(std::map<traceback_path,int>::iterator lhs, std::map<traceback_path,int>::iterator rhs)
    {
        return ((*lhs).second < (*rhs).second);
    }

}