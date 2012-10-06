//
//  trainingSeq.cpp
//  StochHMM
//
//  Created by Paul Lott on 3/5/12.
//  Copyright (c) 2012 University of California, Davis. All rights reserved.
//

#include "trainingSeq.h"

sequence::sequence(string& fileName){
    file=new std::ifstream;
    file->open(fileName.c_str());
    if (!file->is_open()){
        cerr << "Couldn't find file: " << fileName <<endl;
        exit (1);
    }
    masked=false;
    
    return;
}

bool sequence::next(){
    header="";
    DNA="";
    seq.clear();
    
    if (file==NULL || file->eof()){
        return false;
    }
    
    //Find first Header line then import header
    while(file->peek() != '>'){
        std::string temp;
        
        if (file->eof()){
            return false;
        }
        
        getline(*file,temp,'\n');
    }
    getline(*file,header,'\n');
    
    if (file->eof()){
        return false;
    }
    
    std::string temp;
    //get sequence
    while(getline(*file,temp,'\n')){
        
    	for(int j=0;j<temp.size();j++){
            temp[j]=toupper(temp[j]);
            switch (temp[j]){
                case 'A':
                    seq.push_back(0);
                    break;
                case 'C':
                    seq.push_back(1);
                    break;
                case 'G':
                    seq.push_back(2);
                    break;
                case 'T':
                    seq.push_back(3);
                    break;
                default:
                    temp[j]='N';
                    seq.push_back(-1);
            }
    	}
        
        DNA+=temp;
        
    	char nl_peek=file->peek();  // see if we have new sequence header on the next line
    	if (nl_peek=='>'){
    		break;
    	}
    	else if (nl_peek==EOF){
    		break;
    	}
    	else{
    		continue;
    	}
    }
    
    if(masked){
        importMask();
        if (mask.size() != seq.size()){
            cerr << "Mask not the same size as sequence\n";
        }
    }
    
    return true;
}

bool sequence::importMask(){
    mask.clear();
    
    if (file==NULL || file->eof()){
        return false;
    }
    
    std::string temp;
    //Find first Header line then import header
    while(file->peek() != '>'){
        
        if (file->eof()){
            return false;
        }
        
        getline(*file,temp,'\n');
    }
    getline(*file,temp,'\n');
    
    if (file->eof()){
        return false;
    }
    
    //get sequence
    while(getline(*file,temp,'\n')){
        
        char * str=new char[temp.size()+1];
        str[temp.size()]=0;
        memcpy(str, temp.c_str(), temp.size());
        
        char * pch;
        pch = strtok (str," ,.");
        while (pch != NULL)
        {
            mask.push_back(atoi(pch));
            pch = strtok (NULL, " ,.");
        }
        
    	char nl_peek=file->peek();  // see if we have new sequence header on the next line
    	if (nl_peek=='>'){
    		break;
    	}
    	else if (nl_peek==EOF){
    		break;
    	}
    	else{
    		continue;
    	}
    }
    
    return true;
}

//Creates the reverse compliment of a sequence
void sequence::reverseComplement () {
	size_t length = seq.size();
	for (int i=0;i<length;i++){
        switch (seq[i]){
            case 0:
                seq[i]=3; break;
            case 1:
                seq[i]=2; break;
            case 2:
                seq[i]=1; break;
            case 3:
                seq[i]=0; break;
            default:
                seq[i]=-1;
        }
	}
    reverse(seq.begin(),seq.end());
	return;
}


bool trainingSeqs::openFiles(std::string& filename){
    
}


bool trainingSeqs::importMask(){
    
}


//Read upto million characters and determine the alphabet of that sequence
//Return a vector<std::string>
void trainingSeqs::determineAlphabet(size_t iter){
    const std::ifstream* file(seqFile[iter]);
    
    //Import a buffer full of the sequence
    std::string temp;
    std::string undigitized;
    
    while(getline(*file,temp,'\n')){
        if (temp[0]== '>'){
            continue;
        }
        else{
            undigitized+=temp;
        }
        
        if (undigitized.size() >= MAX_BUFFER){
            break;
        }
    }
    
    std::set<std::string> alphabet;
    pair<std::set<int>::iterator,bool> ret;
    
    for(size_t i=0;i<undigitized.size();i++){
        ret = alphabet.insert(undigitized[i]);
        if (ret.second()){
            indexAlpha[iter].push_back(undigitized[i]);
        }
    }
    
    for(size_t i=0;i<indexAlpha[iter].size();i++){
        string &temp= indexAlpha[iter][i];
        alphaIndex[iter][temp]=i;
    }
    
    file.seekg(0,std::ios::beg);
    file.clear();
    
    return;
}

void trainingSeqs::setAlphabet(size_t iter, stringList& lst){
    for(size_t i=0;i<lst.size();i++){
        string& temp = lst[i];
        indexAlpha[iter].push_back(temp);
        alphaIndex[iter][temp]=i;
    }
    return;
}
          

//          bool sequence::getFasta(std::ifstream& file, track* trk){
//              
//              seqtrk=trk;
//              
//              if (!file.good()){
//                  return false;
//              }
//              
//              //Find next header mark 
//              while(file.peek() != '>'){
//                  std::string temp;
//                  getline(file,temp,'\n');
//              }
//              
//              getline(file,header,'\n');
//              
//              //get sequence
//              std::string line;
//              while(getline(file,line,'\n')){
//                  undigitized+=line;
//                  
//                  char nl_peek=file.peek();  // see if we have new sequence header on the next line
//                  if (nl_peek=='>'){
//                      _digitize();
//                      break;
//                  }
//                  else if (nl_peek=='['){
//                      _digitize();
//                      //external=getExDef(file,nextSeq->size());
//                  }
//                  else if (nl_peek==EOF){
//                      _digitize();
//                      break;
//                  }
//                  else{
//                      continue;
//                  }
//              }
//              
//              length=seq->size();
//              return true;
//          }
//          
//          void sequence::_digitize(){
//              stringList lst;
//              clear_whitespace(undigitized,"\n");
//              if (seqtrk->getAlphaMax()>1){
//                  lst.splitString(undigitized, " ,\t");
//              }
//              else{
//                  lst.fromAlpha(undigitized, 1);
//              }
//              
//              for (size_t i=0;i<lst.size();i++){
//                  short symbl = seqtrk->symbolIndex(lst[i]);
//                  seq->push_back(symbl);
//              }
//              
//              return;
//          }
