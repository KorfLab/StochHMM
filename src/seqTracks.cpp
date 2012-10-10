//
//  seqTracks.cpp
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



#include "seqTracks.h"
namespace StochHMM{

    //!Create an empty initialized seqTrack
    seqTracks::seqTracks(){
        _init();
        return;
    }
    
    //Single Model, Single Seq File
    
    //!Create, initialize, and start loading sequence file
    //!\param mdl  Single model
    //!\param filename  Sequence File filename
    //!\param format  Sequence file Format (FASTA or FASTQ)
    seqTracks::seqTracks(model& mdl , std::string& filename, SeqFileFormat format){
        _init();
        
        loadSeqs(mdl, filename, format, NULL);
        
        return;
    }
    
    //!Create, initialize, and start loading sequence file
    //!\param mdl Single model
    //!\param filename  Sequence File filename
    //!\param format  Sequence file Format (FASTA or FASTQ)
    //!\param trFunc TrackFunc* functions to create tracks based on imported sequences
    seqTracks::seqTracks(model& mdl , std::string& filename , SeqFileFormat format, TrackFuncs* trFunc){
        _init();
        
        loadSeqs(mdl, filename,format, trFunc);
        
        return;
    }
    
    
    
    //Single Model, Multiple Seq File
    
    //!Create, initialize, and start loading sequence file
    //!\param mdl  Single model
    //!\param filenames  List of Sequence File filenames
    //!\param format  Sequence file Format (FASTA or FASTQ)
    //!\param type Single track per file (SINGLE) or multiple tracks per file (MULTI)
    seqTracks::seqTracks(model& mdl , std::vector<std::string>& filenames, SeqFileFormat format, SeqFilesType type){
        _init();
        
        loadSeqs(mdl, filenames, format, type, NULL);
        
        return;
    }
    
    
    //!Create, initialize, and start loading sequence file
    //!\param mdl  Single Model
    //!\param filenames  List of Sequence File filenames
    //!\param format  Sequence file Format (FASTA or FASTQ)
    //!\param type Single track per file (SINGLE) or multiple tracks per file (MULTI)
    //!\param trFunc TrackFunc* functions to create tracks based on imported sequences
    seqTracks::seqTracks(model& mdl , std::vector<std::string>& filenames, SeqFileFormat format, SeqFilesType type, TrackFuncs* trFunc){
        _init();
        
        
        loadSeqs(mdl, filenames, format, type, trFunc);
        
        return;
    }

    
    //Multiple Models, Single Seq File

    
    //!Create, initialize, and start loading sequence file
    //!\param mdls Multiple models 
    //!\param filename Sequence File filename
    //!\param format  Sequence file Format (FASTA or FASTQ)
    //!\param attribFunc  pt2Attrib function to choose model based on sequence attributes
    seqTracks::seqTracks(models& mdls , std::string& filename , SeqFileFormat format, pt2Attrib* attribFunc){
        _init();
                
        loadSeqs(mdls,filename,format,attribFunc, NULL);
            
        return; 
    }
    
    
    //!Create, initialize, and start loading sequence file
    //!\param mdls  Multiple models
    //!\param filename Sequence File filename
    //!\param format  Sequence file Format (FASTA or FASTQ)
    //!\param attribFunc  pt2Attrib function to choose model based on sequence attributes
    //!\param trFunc TrackFunc* functions to create tracks based on imported sequences
    seqTracks::seqTracks(models& mdls , std::string& filename , SeqFileFormat format, pt2Attrib* attribFunc, TrackFuncs* trFunc){
        _init();
       
        loadSeqs(mdls,filename,format,attribFunc,trFunc);
        
        return;
        
    }
    
    //Multiple Models, Multiple Seq Files    
    
    //!Create, initialize, and start loading sequence file
    //!\param mdls  Multiple models
    //!\param filenames  List of Sequence File filenames
    //!\param format  Sequence file Format (FASTA or FASTQ)
    //!\param type Single track per file (SINGLE) or multiple tracks per file (MULTI)
    //!\param attribFunc  pt2Attrib function to choose model based on sequence attributes
    seqTracks::seqTracks(models& mdls , std::vector<std::string>& filenames, SeqFileFormat format, SeqFilesType type, pt2Attrib* attribFunc){
        _init();
                
        loadSeqs(mdls, filenames, format, type, attribFunc, NULL);
        
        return;
    }
    
    
    //!Create, initialize, and start loading sequence file
    //!\param mdls  Multiple models
    //!\param filenames  List of Sequence File filenames
    //!\param format  Sequence file Format (FASTA or FASTQ)
    //!\param type Single track per file (SINGLE) or multiple tracks per file (MULTI)
    //!\param attribFunc  pt2Attrib function to choose model based on sequence attributes
    //!\param trFunc TrackFunc* functions to create tracks based on imported sequences
    seqTracks::seqTracks(models& mdls , std::vector<std::string>& filenames, SeqFileFormat format, SeqFilesType type, pt2Attrib* attribFunc, TrackFuncs* trFunc){
        _init();
        
        loadSeqs(mdls, filenames, format, type, attribFunc, trFunc);
        
        return;
    }
    
    
    
    //!Destroy seqTracks
    seqTracks::~seqTracks(){
        for(int i=0;i<filehandles.size();i++){
            if (filehandles[i]!=NULL){
                filehandles[i]->close();
            }
        }
        
        hmms=NULL;
        trackFunctions = NULL;
        attribModelFunc = NULL;
        
        while(!jobQueue.empty()){
            seqJob* element=jobQueue.front();
            jobQueue.pop();
            delete element;
            element = NULL;
        }
        
        return;
    }
    
    
    void seqTracks::_init(){
        numImportJobs=1;
        jobs=0;
        
        hmms    = NULL;
        hmm     = NULL;
        trackFunctions   = NULL;
        attribModelFunc  = NULL;
        
        
        seqFormat = FASTA;  //Set default file format to fasta
        fileType = SINGLE_TRACK;
        
        good=false;
        
        //Make seqTrack thread-safe
#ifdef THREADS
        pthread_mutex_init(&exit_thread_flag_mutex, NULL);
        pthread_cond_init(&exit_thread_flag_cv, NULL);
#endif
        exit_thread=1;
    }
    
    //Reset the Queue, Files, and Import Tracks
    //Does not reset track Functions or Attribute model selection functions 
    void seqTracks::_reset(){
        for(size_t i=0;i<filehandles.size();i++){
            delete filehandles[i];
        }
        
        filehandles.clear();
        seqFilenames.clear();
        
        hmms    = NULL;
        hmm     = NULL;
        seqFormat = FASTA;
        fileType  = SINGLE_TRACK;
        
        good = false;
        
        exit_thread = 1;
        
        importTracks.clear();
    }


    
    
    /////////////////////////////////   Importing Fasta files //////////////////////////////////
    
    
    ///////////////////////////////////////////////////////////////////////////////
    ////////////////  Single Model , Single Sequence File  ////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    
    
    //! Load the fasta sequence file
    //! \param mod  Model to be used
    //! \param seqFile  Fasta sequence filename
    bool seqTracks::loadSeqs(model& mod, std::string& seqFile, SeqFileFormat format){     
        return loadSeqs(mod,seqFile,format,NULL);
    }
    
    
    //! Load the fasta sequence file
    //! \param mod  Model to be used
    //! \param seqFile  Sequence filename
    bool seqTracks::loadSeqs(model& mod, std::string& seqFile, SeqFileFormat format, TrackFuncs* trFuncs){
        
        if (filehandles.size()>0){
            _reset();
        }
        
        hmm = &mod;
        seqFormat = format;
        seqFilenames.push_back(seqFile);
        
        //Assign valid Track Functions
        if (trFuncs!=NULL){
            trackFunctions = trFuncs;
        }
        
        //Get State Information and Determine # of tracks to import
        _initStateInfo();
        fileType = (importTracks.size()>1) ? MULTI_TRACK : SINGLE_TRACK;
        
        _open();
        
        //Fill Job Queue
        importJobs();

        return true;
    }
    
    
    
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////  Single Model , Multiple Sequence File  ///////////////////
    ///////////////////////////////////////////////////////////////////////////////
    
    //! Load the fasta sequence files
    //! \param mod  Model to be used
    //! \param seqFile  Fasta sequence filenames
    //! \param type  SINGLE track per file or MULTI tracks per file
    bool seqTracks::loadSeqs(model& mod, std::vector<std::string>& seqFiles, SeqFileFormat format, SeqFilesType type){
        return loadSeqs(mod, seqFiles, format,type,NULL);
    }
    
    
    //! Load the fasta sequence files
    //! \param mod  Model to be used
    //! \param seqFile  Fasta sequence filenames
    //! \param type SINGLE track per file or MULTI tracks per file
    //! \param trFuncs  Track Functions to create tracks using imported seqs
    bool seqTracks::loadSeqs(model& mod, std::vector<std::string>& seqFiles, SeqFileFormat format, SeqFilesType type, TrackFuncs* trFuncs){
        if (filehandles.size()>0 || importTracks.size()>0){
            _reset();
        }
        
        hmm = &mod;
        seqFormat = format;
        fileType = type;
        seqFilenames = seqFiles;

        
        //Assign valid Track Functions
        if (trFuncs!=NULL){
            trackFunctions = trFuncs;
        }
        
        //Get State Information and Determine # of tracks to import
        _initStateInfo();
        
        
        size_t tracksToImport = importTracks.size();
        if (fileType == SINGLE_TRACK && tracksToImport>1){
            size_t sequenceFiles = seqFiles.size();

            if (tracksToImport!=sequenceFiles){
                std::cerr << "Number of tracks to import and sequenced don't match.  # Files == # Tracks to import " << std::endl;
                _reset();
                return false;
            }
        }
        
        //Open File
        _open();
        
//        for(size_t i = 0; i<seqFiles.size();i++){
//            std::ifstream *SEQ= new std::ifstream;
//            filehandles.push_back(SEQ);
//            
//            filehandles[i]->open(seqFiles[i].c_str());
//            
//            if (!filehandles[i]->is_open()){
//                std::cerr << "Can't open sequence file: "  << seqFiles[i] << std::endl;
//                return false;
//            }
//            
//            if (filehandles[i]->good()){
//                good = true;
//            }
//            else{
//                std::cerr << "Can't read from file: " << seqFiles[i] << std::endl;
//                
//                if (fileType == SINGLE_TRACK && importTracks.size()>1){
//                    std::cerr << "Failed import of " << seqFiles[i] << " causes there to be a missing track in sequence data." << std::endl;
//                    return false;
//                }
//                else{
//                    std::cerr << "Skipped processing of " << seqFiles[i] << "." << std::endl;
//                }
//            }
//        }
                
        //Fill Job Queue
        importJobs();
        
        return true;
    }
    
    
    ///////////////////////////////////////////////////////////////////////////////
    //////////////// Multiple Models , Single Sequence File  //////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    
    
    //! Load the fasta sequence file
    //! \param mModels  Models to be used
    //! \param seqFile  Sequence filename
    //! \param format   Sequence File Format
    //! \param attribFunc   Pointer to Attribute calculation function
    //! \param trFuncs      Pointer to Track Functions
    bool seqTracks::loadSeqs(models &mModels, std::string &seqFile, SeqFileFormat format, pt2Attrib* attribFunc, TrackFuncs* trFuncs){
        
        if (filehandles.size()>0){
            _reset();
        }
        
        hmms = &mModels;
        seqFormat = format;
        seqFilenames.push_back(seqFile);
        
        //Assign valid Attribute Function
        if (attribFunc!=NULL){
            attribModelFunc = attribFunc;
        }
        
        if (attribModelFunc==NULL){
            std::cerr << "No valid Attribute calculating function" << std::endl;
            return false;
        }
        
        //Assign valid Track Functions
        if (trFuncs!=NULL){
            trackFunctions = trFuncs;
        }
        
        //Get State Information and Determine # of tracks to import
        _initStateInfo();
        fileType = (importTracks.size()>1) ? MULTI_TRACK : SINGLE_TRACK;
        
        
        //Open File 
        _open();
//        std::ifstream *SEQ= new std::ifstream;
//        filehandles.push_back(SEQ);
//        
//        filehandles[0]->open(seqFile.c_str());
//        
//        if (!filehandles[0]->is_open()){
//            std::cerr << "Can't open sequence file: "  << seqFile << std::endl;
//            return false;
//        }
//        
//        if (filehandles[0]->good()){
//            good = true;
//        }
//        else{
//            std::cerr << "Can't read from file: " << seqFile << std::endl;
//            return false;
//        }
        
        //Fill Job Queue
        importJobs();
        
        
        return true;
    }
    
    
    //! Load the fasta sequence file
    //! \param mModels  Models to be used
    //! \param seqFile  Fasta sequence filename
    bool seqTracks::loadSeqs(models &mModels, std::string &seqFile, SeqFileFormat format){
        return loadSeqs(mModels, seqFile, format, NULL,NULL);
    }
    
    
    ///////////////////////////////////////////////////////////////////////////////
    ////////////////// Multiple Models , Multiple Sequence File  //////////////////
    ///////////////////////////////////////////////////////////////////////////////
    
    
    //! Load the fasta sequence files
    //! \param mModels  Models to be used
    //! \param seqFiles  Fasta sequence filenames
    bool seqTracks::loadSeqs(models& mModels, std::vector<std::string>& seqFiles, SeqFileFormat format, SeqFilesType type, pt2Attrib* attribFunc, TrackFuncs* trFuncs){
        
        if (filehandles.size()>0){
            _reset();
        }
        
        hmms = &mModels;
        seqFormat = format;
        fileType = type;
        seqFilenames = seqFiles;

        
        //Assign valid Track Functions
        if (trFuncs!=NULL){
            trackFunctions = trFuncs;
        }
        
        //Assign valid Attribute Function
        if (attribFunc!=NULL){
            attribModelFunc = attribFunc;
        }
        
        if (attribModelFunc==NULL){
            std::cerr << "No valid Attribute calculating function" << std::endl;
            return false;
        }
        
        //Get State Information and Determine # of tracks to import
        _initStateInfo();
        
        size_t tracksToImport = importTracks.size();
        if (fileType == SINGLE_TRACK && tracksToImport>1){
            size_t sequenceFiles = seqFiles.size();
            
            if (tracksToImport!=sequenceFiles){
                std::cerr << "Number of tracks to import and sequenced don't match.  # Files == # Tracks to import " << std::endl;
                return false;
            }
        }
        
        //Open File
        _open();
//        for(size_t i = 0; i<seqFiles.size();i++){
//            std::ifstream *SEQ= new std::ifstream;
//            filehandles.push_back(SEQ);
//            
//            filehandles[i]->open(seqFiles[i].c_str());
//            
//            if (!filehandles[i]->is_open()){
//                std::cerr << "Can't open sequence file: "  << seqFiles[i] << std::endl;
//                return false;
//            }
//            
//            if (filehandles[i]->good()){
//                good = true;
//            }
//            else{
//                std::cerr << "Can't read from file: " << seqFiles[i] << std::endl;
//                
//                if (fileType == SINGLE_TRACK && importTracks.size()>1){
//                    std::cerr << "Failed import of " << seqFiles[i] << " causes there to be a missing track in sequence data." << std::endl;
//                    return false;
//                }
//                else{
//                    std::cerr << "Skipped processing of " << seqFiles[i] << "." << std::endl;
//                }
//            }
//        }
        
        //Fill Job Queue
        importJobs();
        
        return true;
    }

    
    
    //! Load the fasta sequence files
    //! \param mModels  Models to be used
    //! \param seqFiles  Fasta sequence filenames
    bool seqTracks::loadSeqs(models& mModels, std::vector<std::string>& seqFiles , SeqFileFormat format, SeqFilesType type){
        return loadSeqs(mModels, seqFiles, format, type, NULL, NULL);
    }
    
    
    //! Load the fasta sequence files
    //! \param mModels  Models to be used
    //! \param seqFiles  Fasta sequence filenames
    bool seqTracks::loadSeqs(models& mModels, std::vector<std::string>& seqFiles, SeqFileFormat format, SeqFilesType type, pt2Attrib* attribFunc){
        return loadSeqs(mModels, seqFiles, format, type,attribFunc,NULL);
    }

    
    //!Get the next sequence(s) and model from the job queue
    //!If the number of jobs falls below MIN_JOBS/2 then refill the queue
    //!Else give pointer to the next job in the queue
    seqJob* seqTracks::getJob(){
        seqJob *jb= NULL;
        
        //Check job queue and if necessary fill the job queue
        if ((jobs==0 || (jobs < (numImportJobs/2))) && good){
            while(good && jobs<numImportJobs){
                getNext();
            }
        }
        
        if (!good){
            if (fileType==SINGLE_TRACK){
                if (seqFilenames.size()>importTracks.size()){
                    _close();
                    _open();
                }
                else{
                    _close();
                }
            }
            else {
                if (seqFilenames.size()>0){
                    _close();
                    _open();
                }
                else{
                    _close();
                }
            }
        }
        
        if (jobs>0){
            jb=jobQueue.front();
            jobQueue.pop();
            jobs--;
        }
        
        return jb;
    }
    
      
    
    
    // TODO: fix PT2TRACKFUNC function assignment
    
    bool seqTracks::_initStateInfo(){
        model* temp = NULL;
        
        if (hmms!=NULL){
            temp = (*hmms)[0];
            
            if (temp == NULL){
                std::cerr << "seqTracks initialization error:  Model not defined at index zero of models datatype.  Can't initialize seqTrack with necessary model information" << std::endl;
                return false;
            }
            
        }
        else if (hmm!=NULL){
            temp = hmm;
        }
        else{
            std::cerr << "seqTracks initialization error:  Model is not defined. Therefore, can't initiate seqTrack with necessary model inforamation" << std::endl;
            return false;
        }
        
        modelTracks = temp->getTracks();
        
        //Get State Names, Labels and Paths for ExDefs
        for(size_t i=0;i<temp->state_size();i++){
            
            info.gff[temp->getStateGFF(i)].push_back(i);
            info.label[temp->getStateLabel(i)].push_back(i);
            info.names[temp->getStateName(i)]=i;
        }
        
        
        //Determine which tracks to import and which to get by using track functions
        track* tempTrack;
        trackCount= temp->track_size();
        for(size_t i=0;i<trackCount;i++){
            tempTrack = temp->getTrack(i);
            if (tempTrack->isTrackFuncDefined()){
                ppTrack tmp;
                tmp.func = NULL;
                tmp.trackNumber=i;
                tmp.trackToUse=temp->getTrackIter(tempTrack->getTrackToUse());
                std::string functionTouse=tempTrack->getTrackFunction();
                
                //tmp.func = funcs->getFunction(functionTouse);
                postprocessTracks.push_back(tmp);
            }
            else{
                importTracks.push_back(std::make_pair(i,tempTrack->getAlphaType()));
            }
        }
        return true;
    }

    //!Print the seqTracks to stdout
    void seqTracks::print(){
        
        for (int i=0;i<jobQueue.size();i++){
            //jobQueue[i]->print_seq();
        }
        return;
    }


    //!Get a job and Add to the queue
    bool seqTracks::importJobs(){
        while(good && jobs<numImportJobs){
            getNext();
        }
        
        if (!good && seqFilenames.size()>importTracks.size()){
            _close();
            _open();
        }
        
        return true;
    }


    //For each track it will get the required sequence type and put it in the sequences 
    
    //TODO: Handle tracks that are created using functions
    
    //!Get the next sequence job
    bool seqTracks::getNext(){
        seqJob* temp_job=new(std::nothrow) seqJob(trackCount);
        
        if (temp_job==NULL){
            std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
            exit(1);
        }
        
        bool valid=true;
        
        sequence* sq;
        
        //std::cout << importTracks.size() << std::endl;
        
        for(int i=0;i<importTracks.size();i++){
            bool success;

            if (importTracks[i].second == REAL){
                sq=new(std::nothrow) sequence(true);
                
                if (sq==NULL){
                    std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                    exit(1);
                }
                
                if (fileType == SINGLE_TRACK){
                    success = sq->getReal(*filehandles[i], (*modelTracks)[importTracks[i].first]);
                }
                else{
                    success = sq->getReal(*filehandles[0], (*modelTracks)[importTracks[i].first]);
                }
                
            }
            else if (seqFormat == FASTA){ // AlphaNum and Fasta
                sq=new(std::nothrow) sequence(false);
                
                if (sq==NULL){
                    std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                    exit(1);
                }
                
                if (fileType == SINGLE_TRACK){
                    success = sq->getFasta(*filehandles[i], (*modelTracks)[importTracks[i].first]);
                }
                else{
                    success = sq->getFasta(*filehandles[0], (*modelTracks)[importTracks[i].first]);
                }             
            }
            else{
                sq=new(std::nothrow) sequence(false);
                if (sq==NULL){
                    std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                    exit(1);
                }
                if (fileType == SINGLE_TRACK){
                    success = sq->getFastq(*filehandles[i], (*modelTracks)[importTracks[i].first]);
                }
                else{
                    success = sq->getFastq(*filehandles[0], (*modelTracks)[importTracks[i].first]);
                }
                
            }
            
            
            if (fileType == SINGLE_TRACK){
                if (!filehandles[i]->good()){
                    good=false;
                }
            }
            else{
                if (!filehandles[0]->good()){
                    good=false;
                }
            }
            
            if (!success){
                std::cerr << "Failed to import data track from " << seqFilenames[i] << std::endl;
                delete sq;
                sq = NULL;
            }
            
            if (sq==NULL){ // If sequence is bad break 
                valid=false;
                break;
            }
            else{
                //If exDef is defined in sequence put it in sequences
                if (sq->exDefDefined()){
                    temp_job->set->setExDef(sq->getExDef());
                }
                temp_job->set->addSeq(sq,importTracks[i].first);
                
                if (fileType==SINGLE_TRACK){
                    temp_job->setSeqFilename(seqFilenames[i]);
                }
                else{
                    temp_job->setSeqFilename(seqFilenames[0]);
                }
            }
        }
        
        
        if (valid){
            
            //Get sequences defined by sequence external function that is user-defined
           for (size_t i =0;i<postprocessTracks.size();i++) {
               std::vector<double>* rl = NULL;
               if (postprocessTracks[i].func != NULL ){
                   rl = (*postprocessTracks[i].func)(temp_job->set->getUndigitized(postprocessTracks[i].trackToUse));
               }
               else{
                   std::cerr << "Sequence external function not defined for track number: " << postprocessTracks[i].trackNumber << std::endl;
                   std::cerr << "Using Sequences from track: " << postprocessTracks[i].trackToUse << std::endl;
                   rl = new(std::nothrow) std::vector<double>;
                   if (rl==NULL){
                       std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                       exit(1);
                   }
               }
                
               sequence* sq = new(std::nothrow) sequence(rl , (*modelTracks)[postprocessTracks[i].trackNumber]);
               
               if (sq==NULL){
                   std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                   exit(1);
               }
               
               temp_job->set->addSeq(sq,postprocessTracks[i].trackNumber);
            }
            
            //TODO:  Fix the selection of models based upon the attribute model function
            //Select Model based on models distance from attributes
            if (hmms){
                
                //If we have a attribute calculation function we'll access which model to use
                if (attribModelFunc){
                    //Calculate attribute
                    double attb = (*attribModelFunc)(temp_job->set->getUndigitized(TrackToUseForAttrib));
                    
                    //assign first by default
                    double min = (*hmms)[0]->getDistanceToAttrib(attb);
                    temp_job->hmm=(*hmms)[0];  
                    
                    //Check other models for one that is closer to attb
                    for(size_t i=1;i<hmms->size();i++){
                        double newVal=(*hmms)[i]->getDistanceToAttrib(attb);
                        
                        //If it's closer assign it to the job
                        if (newVal<min){
                            temp_job->hmm=(*hmms)[i];
                        }
                    }

                }
                else{
                    temp_job->hmm=(*hmms)[0];  //Default to first HMM
                }
            }
            
            else{
                temp_job->hmm=hmm;  //assign single model to job
            }
            
            
            //Check that all sequences are same length.
            size_t lengthOfAll=SIZE_MAX;
            for (int i=0;i<trackCount;i++){
                
                size_t length=temp_job->set->getLength(i);
                if (lengthOfAll==SIZE_MAX){
                    lengthOfAll=length;
                }
                else if (lengthOfAll!=length){
                    std::cerr << "Sequence Lengths not the same" <<std::endl;
                    delete temp_job;
                    return false;
                }
                else {
                    continue;
                }

            }
            temp_job->set->setLength(lengthOfAll);
            jobQueue.push(temp_job);
            jobs++;
        }
        else{
            delete temp_job;
        }
        
        return true;
    }
    
    
    bool seqTracks::_open(){
        
        if (seqFilenames.size()==0){
            return false;
        }
        
        if (fileType == SINGLE_TRACK ){
            if (seqFilenames.size()<importTracks.size()){
            std::cerr << "Number of sequences provided doesn't match the number of tracks. Given that the sequences files contain a single track per file. " << std::endl;
            }
            
            for(size_t i = 0; i<importTracks.size();i++){
                std::ifstream *SEQ= new(std::nothrow) std::ifstream;
                
                if (SEQ==NULL){
                    std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                    exit(1);
                }
                
                filehandles.push_back(SEQ);
                
                if (seqFilenames.size()<i+1){
                    return false;
                }
                
                filehandles[i]->open(seqFilenames[i].c_str());
                
                if (!filehandles[i]->is_open()){
                    std::cerr << "Can't open sequence file: "  << seqFilenames[i] << std::endl;
                    return false;
                }
                
                if (filehandles[i]->good()){
                    good = true;
                }
                else{
                    std::cerr << "Can't read from file: " << seqFilenames[i] << std::endl;
                    
                    if (fileType == SINGLE_TRACK && importTracks.size()>1){
                        std::cerr << "Failed import of " << seqFilenames[i] << " causes there to be a missing track in sequence data." << std::endl;
                        return false;
                    }
                    else{
                        std::cerr << "Skipped processing of " << seqFilenames[i] << "." << std::endl;
                    }
                }
            }
        }
        else {
            
            std::ifstream *SEQ= new(std::nothrow) std::ifstream;
            
            if (SEQ==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
            
            filehandles.push_back(SEQ);
            
            filehandles[0]->open(seqFilenames[0].c_str());
            
            if (!filehandles[0]->is_open()){
                std::cerr << "Can't open sequence file: "  << seqFilenames[0] << std::endl;
                return false;
            }
            
            if (filehandles[0]->good()){
                good = true;
            }
            else{
                std::cerr << "Can't read from file: " << seqFilenames[0] << std::endl;
                return false;
            }
        }
        
        
                
        return true;
    }
                   
                   
    bool seqTracks::_close(){
        if (fileType == SINGLE_TRACK){
            for(size_t i=0;i<importTracks.size();i++){
                //std::cout << importTracks.size() << "\t" << filehandles.size() << std::endl;

                if (filehandles.size()>=importTracks.size()){
                    filehandles[i]->close();
                    delete filehandles[i];
                    filehandles[i]=NULL;
                    filehandles.erase(filehandles.begin());
                    seqFilenames.erase(seqFilenames.begin());
                }
                
            }
        }
        else{
            filehandles[0]->clear();
            delete filehandles[0];
            filehandles.erase(filehandles.begin());
            seqFilenames.erase(seqFilenames.begin());
        }
        
        return true;
    }
                    
}
