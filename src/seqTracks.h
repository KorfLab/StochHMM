//
//  seqTracks.h
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

#ifndef SEQTRACK_H
#define SEQTRACK_H

#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <queue>
#include "text.h"
#include "traceback_path.h"
#include "hmm.h"
#include "sequences.h"
#include <stdlib.h>




namespace StochHMM{
    
        
#ifdef THREADS
    extern pthread_cond_t exit_thread_flag_cv;
    extern pthread_mutex_t exit_thread_flag_mutex;
#endif 
    
    //!\file seqTracks.h
    //! Contains functions to import FASTA/FASTQ sequences from files and select the applicable model to deal with that sequence.
    //! It was set up to generate a seqJob for each sequence, then select a model, and then allow the programmer to thread the decoding algorithm.
    
    //!\enum enum SeqFileFormat {FASTA,  FASTQ};
    //!File format of the sequences
    enum SeqFileFormat {FASTA, FASTQ, CSV};
    
    //!\enum enum SeqFilesType {SINGLE, MULTI};
    //!Sequence files have single track or multiple track sequences per file
    enum SeqFilesType {SINGLE_TRACK, MULTI_TRACK};
    
    
    class seqJob;
    

    
    //!\struct ppTrack
    //!Stores what track is determined by a Track Function
    //!trackNumber is the index reference the derived track
    //!trackToUse is the track to use to generate the derived track
    //!Example: Function would be SIDD, the track we would use to get the sidd track
    //!would be a DNA sequence track.
    struct ppTrack{
        int trackNumber;
        int trackToUse;
        StochHMM::pt2TrackFunc* func;
    };

    
    /*! \class SeqTracks
        \brief SeqTracks are used to integrate the information provided by the model with the sequences that are being imported
    */
    class seqTracks{
    public:
        
        //Constructor
        seqTracks();
        
        //Single Model, Single Seq File
        seqTracks(model&, std::string&, SeqFileFormat);
        seqTracks(model&, std::string&, SeqFileFormat, TrackFuncs*);
        
        //Single Model, Multiple Seq File)
        seqTracks(model&, std::vector<std::string>&, SeqFileFormat, SeqFilesType);
        seqTracks(model&, std::vector<std::string>&, SeqFileFormat, SeqFilesType, TrackFuncs*);
        
        //Multiple Models, Single Seq File
        seqTracks(models&,std::string&, SeqFileFormat, pt2Attrib*);
        seqTracks(models&,std::string&, SeqFileFormat, pt2Attrib*, TrackFuncs*);
        
        
        //Multiple Models, Multiple Seq Files
        seqTracks(models&,std::vector<std::string>&, SeqFileFormat, SeqFilesType, pt2Attrib*);
        seqTracks(models&,std::vector<std::string>&, SeqFileFormat, SeqFilesType, pt2Attrib*, TrackFuncs*);
        
        
        //Destructor
        ~seqTracks();
        
        //MUTATORS
                
        ////////////////  Single Model , Single Sequence File  ////////////////////////
        bool loadSeqs(model&, std::string&, SeqFileFormat); // <-Main Function
        bool loadSeqs(model&, std::string&, SeqFileFormat, TrackFuncs*); // <-Main Function
        
        
        ////////////////  Single Model , Multiple Sequence File  //////////////////////
        bool loadSeqs(model&, std::vector<std::string>&, SeqFileFormat, SeqFilesType);
        bool loadSeqs(model&, std::vector<std::string>&, SeqFileFormat, SeqFilesType, TrackFuncs*);  // <- Main Function

        
        ////////////////  Multiple Models , Single Sequence File  //////////////////////
        
        bool loadSeqs(models&, std::string&, SeqFileFormat, pt2Attrib*, TrackFuncs*); // <-Main Function
        
        bool loadSeqs(models&, std::string&, SeqFileFormat); //only allow if pt2Attrib is set else error
        
        
        ////////////////  Multiple Models , Multiple Sequence File  ////////////////////
        
               
        //// Main Function ////
        bool loadSeqs(models&, std::vector<std::string>&, SeqFileFormat, SeqFilesType, pt2Attrib*, TrackFuncs*); // <-Main Function
        //// Facade Functions ////
        bool loadSeqs(models&, std::vector<std::string>&, SeqFileFormat, SeqFilesType); //only allow if pt2Attrib is set else error
        bool loadSeqs(models&, std::vector<std::string>&, SeqFileFormat, SeqFilesType, pt2Attrib*);
        

        
        inline void setAttribFunc(pt2Attrib* func){attribModelFunc=func;}
        inline void setTrackFunc(TrackFuncs* func){trackFunctions=func;}
        inline void setNumImportJobs(size_t value){numImportJobs=value;}
        
        void setTrackFilename(std::string&,std::string&);

        //ACCESSORS
        bool getNext();
        bool importJobs();
        
        seqJob* getJob();

        sequence* getFasta(int);
        sequence* getFastq(int);
        sequence* getReal(int);
        
        size_t size(void){return jobQueue.size();}
        inline int getTrackCount(){return trackCount;}
        
        void print();
        void getInformation();
        float remainingSeqs(){
            float val = (float) seqFilenames.size() / (float) importTracks.size();
            //std::cout << seqFilenames.size() <<"\t" <<  importTracks.size() <<"\t" <<  val << std::endl;
            return val;}
        
    private:
        
        std::vector<std::ifstream*> filehandles; //input file stream
        std::vector<std::string> seqFilenames; //input filenames
        size_t numImportJobs;
        bool good;
        
        TrackFuncs* trackFunctions;
        pt2Attrib* attribModelFunc;
        
        
        int TrackToUseForAttrib;
        
        stateInfo info;
        
        std::vector<std::pair<int,trackType> > importTracks;
        std::vector<ppTrack> postprocessTracks;
        
        
        models* hmms; //Models
        model* hmm; //Models
        tracks* modelTracks; //Tracks defined by models
        
        SeqFileFormat seqFormat;  //File format (Fasta or FastQ);
        SeqFilesType fileType;    //File Type (Single File or Multiple Files);
        
        
        //bool fastq;  
        //bool multiFile;
        
        size_t trackCount;
        stringList input;
        
        
        //Threading Variables
        std::queue<seqJob*> jobQueue; //used to be trcks
        int jobs;  //Counts of # of jobs waiting
        int exit_thread;  //set to 0 if file stream is EOF
        
        
        //External Definition import function for Sequence
        ExDefSequence* getExDef(int,int);

        bool _loadFASTA(std::string&, SeqFilesType);
        bool _loadFASTA(std::vector<std::string>&, SeqFilesType);
        
        bool _loadFASTQ(std::string&, SeqFilesType);
        bool _loadFASTQ(std::vector<std::string>&, SeqFilesType);
        
        
        bool _initStateInfo(void);
        void _reset();
        void _init();
        bool _open();
        bool _close();
    };
    
    
    

    
        
    //!\class seqJob
    //!Stores the model and sequence information for each job
    class seqJob{   //Could make a derivative of sequences
    public:
        //Constructor
        seqJob(size_t);
        
        //Destructor
        ~seqJob();
        
        friend class seqTracks;
        
        
        //MUTATORS
        void evaluateFunctions();
        
        
        //ACCESSORS
        inline size_t size(){return set->getLength();};
        
        inline model* getModel(){return hmm;};
        inline sequences* getSeqs(){return set;};
        
        inline std::string getHeader(){return set->getHeader();};
        
        inline bool evaluated(){return funcEvaluated;};
        
        
        inline void printModel(){if(hmm) hmm->print();};
        
        inline void printSeq(){set->print();};
        
        inline traceback_path* getPath(){if (decodingPerformed) return path;else return NULL;};
        
        double getSeqAttrib(){return attrib;};
        
        inline std::string getSeqFilename(size_t iter){return seqFilename[iter];};
        inline void setSeqFilename(std::string& filename){seqFilename.push_back(filename); return;};
        inline void printFilenames(){for(size_t i=0;i<seqFilename.size();i++){ std::cout << seqFilename[i] << std::endl;}};
        
    private:
        model* hmm;
        sequences* set;
        
        std::vector< std::string>  seqFilename;
        
        double attrib;
        bool funcEvaluated;
        
        TrackFuncs* functions;
        
        bool decodingPerformed;
        traceback_path* path;
        
    };


}
#endif /*SEQTRACK_H*/