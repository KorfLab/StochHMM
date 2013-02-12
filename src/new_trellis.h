//
//  new_trellis.h
//  StochHMM
//
//  Created by Paul Lott on 11/13/12.
//  Copyright (c) 2012 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#ifndef __StochHMM__new_trellis__
#define __StochHMM__new_trellis__

#include <iostream>
#include <vector>
#include <stdint.h>
#include <vector>
#include <string>
#include <fstream>
#include <bitset>
#include <map>
#include "stochTypes.h"
#include "sequences.h"
#include "hmm.h"
#include "traceback_path.h"
#include "stochMath.h"
#include <stdint.h>
#include <iomanip>
#include "stochTable.h"

namespace StochHMM{
	
	
	typedef std::vector<std::vector<int16_t> > int_2D;
	typedef std::vector<std::vector<float> > float_2D;
	typedef std::vector<std::vector<double> > double_2D;
	typedef std::vector<std::vector<long double> > long_double_2D;

	typedef std::vector<std::vector<std::vector<uint16_t> > > int_3D;
	typedef std::vector<std::vector<std::vector<float> > > float_3D;
	typedef std::vector<std::vector<std::vector<double> > > double_3D;
	typedef std::vector<std::vector<std::vector<long double> > > long_double_3D;
//	typedef std::vector<std::vector<std::vector<std::pair<int16_t,int16_t> >
	
	class nthScore{
	public:
		int16_t st_tb;
		int16_t score_tb;
		double score;
		nthScore():st_tb(0),score_tb(0),score(-INFINITY){};
		nthScore(int16_t st, int16_t tb, double sc):st_tb(st),score_tb(tb),score(sc){};
		
	};
	
//	struct nthTrace{
//		int16_t st_tb;
//		int16_t score_tb;
//		nthTrace():st_tb(0),score_tb(0){};
//		nthTrace(int16_t st, int16_t tb):st_tb(st),score_tb(tb){};
//	};
	
	class nthTrace{
	public:
		std::map<int32_t,int32_t> tb;
		
		inline void assign(int16_t st, int16_t n_score, int16_t tb_ptr, int16_t tb_sc){
			tb[(((int32_t) st) << 16 | n_score)] = (((int32_t)tb_ptr) << 16) | tb_sc;
		}
		
		inline void get(int16_t& st, int16_t& n_score){
			int32_t key_val = (((int32_t) st) << 16 | n_score);
			if (tb.count(key_val)){
				st = -1;
			}
			int32_t val = tb[key_val];
			st	= (val >> 16) & 0xFFFF;
			n_score	= val & 0xFFFF;
		}
	};

	
	/*! \class Trellis
	 *	\brief Implements the HMM scoring trellis and algorithms
	 *
	 *
	 *
	 */
	class trellis{
	public:
		trellis();
		trellis(model* h , sequences* sqs);
		~trellis();
		void reset();
		inline model* getModel(){return hmm;}
		inline sequences* getSeq(){return seqs;}
		
		/*-----------   Decoding Algorithms ------------*/
		
		//TODO: Fix these functions so that they evaluate the model and choose a
		// a default algorithm based on the whether the model defines explicit duration
		// states
		
		void viterbi();
		void viterbi(model* h, sequences* sqs);
		
		void forward();
		void forward(model* h, sequences* sqs);
		
		void forward_viterbi();
		void forward_viterbi(model* h, sequences* sqs);
		
		void backward();
		void backward(model* h, sequences* sqs);
		
		void posterior();
		void posterior(model* h, sequences* sqs);
		
		void stochastic_viterbi();
		void stochastic_viterbi(model* h, sequences* sqs);
		
		void stochastic_forward();
		void stochastic_forward(model* h, sequences* sqs);
		
		void nth_viterbi();
		void nth_viterbi(model* h, sequences* sqs);
		
		void baum_welch();
		
		
		/*-----------	Naive Model Decoding Algorithms -------------*/
		/* These algorithms are the simply coded algorithms with little to no
		 optimizations.  
		 */
		
		void naive_forward();
		void naive_forward(model* h, sequences* sqs);
		
		void naive_backward();
		void naive_backward(model* h, sequences* sqs);
		
		void naive_viterbi();
		void naive_viterbi(model* h, sequences* sqs);
		
		void naive_baum_welch();
		void naive_baum_welch(model* h, sequences* sqs);
		
		void naive_stochastic_viterbi();
		void naive_stochastic_viterbi(model* h, sequences* sqs);
		
		void naive_stochastic_forward();
		void naive_stochastic_forward(model* h, sequences* sqs);
		
		void naive_nth_viterbi(size_t n);
		void naive_nth_viterbi(model* h, sequences* sqs, size_t n);
		
		/*-----------   Simple Model Decoding Algorithms ------------*/
		/* These algorithms are for use with models that do not define explicit
			duration states or external functions.   Therefore, probabilities can
			be easily computed or referenced.
		 
		 */
		
		void simple_viterbi();
		void simple_viterbi(model* h, sequences* sqs);
				
		void simple_forward_viterbi();
		void simple_forward_viterbi(model* h, sequences* sqs);
		
		void simple_forward();
		void simple_forward(model* h, sequences* sqs);
		
		void simple_backward();
		void simple_backward(model* h, sequences* sqs);
		
		void simple_posterior();
		void simple_posterior(model* h, sequences* sqs);
				
		void simple_stochastic_viterbi();
		void simple_stochastic_viterbi(model* h, sequences* sqs);
		
		void simple_stochastic_forward();
		void simple_stochastic_forward(model* h, sequences* sqs);
		
		void simple_baum_welch();
		void simple_baum_welch(model* h, sequences* sqs);
		
		void simple_nth_viterbi(size_t n);
		void simple_nth_viterbi(model* h, sequences* sqs, size_t n);

		
		/*-----------   Fast Complex Model Decoding Algorithms  ----------*/
		/*	These algorithms are for use with models that define external functions
			or explicit duration states.
		 
			These algorithms store the associated emissions and transitions
			probabilities that were calculated during viterbi algorith.  For later
			use in the forward/backward algorithms.
		 
			The values are stored in dense arrays, so the memory is accessed faster
			at the cost of amount of memory required.
		 */
		
		void fast_complex_viterbi();
		void fast_complex_viterbi(model* h, sequences* sqs);
		
		void fast_complex_forward_viterbi();
		void fast_complex_forward_viterbi(model* h, sequences* sqs);
		
		void fast_complex_backward();
		void fast_complex_backward(model* h, sequences* sqs);
		
		void fast_complex_stochastic_viterbi();
		void fast_complex_stochastic_viterbi(model* h, sequences* sqs);
		
		void fast_complex_stochastic_forward();
		void fast_complex_stochastic_forward(model* h, sequences* sqs);
		
		void fast_complex_baum_welch();
		void fast_complex_baum_welch(model* h, sequences* sqs);
		
		
		/*-----------   Sparse Complex Model Decoding Algorithms  --------*/
		/*	These algorithms store the associated emissions and transitions
			probabilities that were calculated during viterbi algorithm.  For later
			use in the forward/backward algorithms.
		
			The values are stored in sparse tables (hash map), so access to values 
			is slower but with lower amount of memory required.
		 */
		
		void sparse_complex_viterbi();
		void sparse_complex_viterbi(model* h, sequences* sqs);
		
		void sparse_complex_foward_viterbi();
		void sparse_complex_foward_viterbi(model* h, sequences* sqs);
		
		void sparse_complex_backward();
		void sparse_complex_backward(model* h, sequences* sqs);
		
		void sparse_complex_stochastic_forward();
		void sparse_complex_stochastic_forward(model* h, sequences* sqs);
		
		void sparse_complex_baum_welch();
		void sparse_complex_baum_welch(model* h, sequences* sqs);
		
		
		/*---------   Standard Traceback Algorithms --------*/
				
		void traceback(traceback_path& path);
        void traceback(traceback_path& path ,size_t position, size_t state);
		void traceback_nth(traceback_path& path, size_t n);
		void traceback_posterior(traceback_path& path);
		
		void traceback_stoch_posterior(traceback_path& path);
		void traceback_stoch_posterior(multiTraceback&, size_t reps);
		void traceback_stoch_forward(multiTraceback&,size_t);
		void traceback_stoch_viterbi(multiTraceback&,size_t);
		void traceback_nth_viterbi(multiTraceback&);
		
		void stochastic_traceback(traceback_path& path);

		
		
		
		inline bool store(){return store_values;}
		inline void store(bool val){store_values=val; return;}

		void print();
		std::string stringify();
		void export_trellis(std::ifstream&);
		void export_trellis(std::string& file);
		inline model* get_model(){return hmm;}
		
		//Model Baum-Welch Updating
		void update_transitions();
		void update_emissions();
		
		
		/*----------- Accessors ---------------*/
		inline double_2D* get_naive_forward_scores(){return dbl_forward_score;}
		inline double_2D* get_naive_backward_scores(){return dbl_backward_score;}
		inline double_2D* get_naive_viterbi_scores(){return dbl_viterbi_score;}
		inline double_2D* get_dbl_posterior(){return dbl_posterior_score;}
		
		
		inline float_2D* getForwardTable(){return forward_score;}
		inline float_2D* getBackwardTable(){return backward_score;}
		inline float_2D* getPosteriorTable(){return posterior_score;}
		
		inline double getForwardProbability(){return ending_forward_prob;}
		inline double getBackwardProbability(){return ending_backward_prob;}

		
	private:
		double getEndingTransition(size_t);
        double getTransition(state* st, size_t trans_to_state, size_t sequencePosition);
        size_t get_explicit_duration_length(transition* trans, size_t sequencePosition,size_t state_iter, size_t to_state);
        double exFuncTraceback(transitionFuncParam*);
		
		
		model* hmm;		//HMM model
        sequences* seqs; //Digitized Sequences
		size_t nth_size; //Size of N to calculate;
		
		
		size_t state_size;	//Number of States
		size_t seq_size;	//Length of Sequence
		
		trellisType type;
		
		bool store_values;
		bool exDef_defined;
		
		//Traceback Tables
		int_2D*		traceback_table;	//Simple traceback table
//		int_3D*		nth_traceback_table;//Nth-Viterbi traceback table
		stochTable* stochastic_table;
		
		//Score Tables
		float_2D*	viterbi_score;      //Storing Viterbi scores
		float_2D*	forward_score;      //Storing Forward scores
		float_2D*	backward_score;     //Storing Backward scores
		float_2D*	posterior_score;	//Storing Posterior scores
		
		double_2D*  dbl_forward_score;
		double_2D*	dbl_viterbi_score;
		double_2D*  dbl_backward_score;
		double_2D*	dbl_posterior_score;
		double_3D*  dbl_baum_welch_score;
		std::vector<std::vector<std::vector<nthScore >* > >* naive_nth_scores;
		
		
		std::vector<nthTrace>*  nth_traceback_table;
		
		//Ending Cells
		double	ending_viterbi_score;
		int16_t	ending_viterbi_tb;
		double	ending_forward_prob;
		double	ending_backward_prob;
		std::vector<nthScore>* ending_nth_viterbi;

//		std::vector<tb_score>*    ending_stoch_tb;
//		std::vector<tb_score>*    ending_nth_viterbi;
		
		//Cells used for scoring each cell
		std::vector<double>* scoring_current;
		std::vector<double>* scoring_previous;
		std::vector<double>* alt_scoring_current;
		std::vector<double>* alt_scoring_previous;
		
		std::vector<size_t>* explicit_duration_current;
		std::vector<size_t>* explicit_duration_previous;
		std::vector<size_t>* swap_ptr_duration;
		
		std::vector<double>* swap_ptr;
		
		std::vector<std::vector<double>* > complex_emissions;
		std::vector<std::vector<std::map<uint16_t,double>* >* >* complex_transitions;
		
		std::vector<std::vector<nthScore> >* nth_scoring_current;
		std::vector<std::vector<nthScore> >* nth_scoring_previous;
		std::vector<std::vector<nthScore> >* nth_swap_ptr;
	};
	
	void sort_scores(std::vector<nthScore>& nth_scores);
	bool _vec_sort(const nthScore& i, const nthScore& j);
	
}

#endif /* defined(__StochHMM__new_trellis__) */
