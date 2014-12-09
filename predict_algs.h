/* MTD program. Copyright (C) 2012. Created by Maria Chatzou and Pantelis Bagos.
   This program is free software: you can redistribute it and/or modify 
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.For more details see
   the GNU General Public License <http://www.gnu.org/licenses/>. */


#ifndef _predict_algs_h
#define _predict_algs_h


void score_pred(char **testseqArray, int alphabet, int dimension, int order, int seq_num, double Lambda[order], double a[order][dimension][dimension], double Score[seq_num]);


void null_Score(char **testseqArray, int alphabet, int dimension, int seq_num, double null_Score[seq_num]);


void null_model_score(char **testseqArray, int alphabet, int dimension, int order, int seq_num, double Lambda[order],double a[order][dimension][dimension],\
				 double null_Lambda[order],double null_a[order][dimension][dimension], double Score[seq_num], double null_Score[seq_num]);


void predict_fastaSet(char *file_name, char **testseqArray, char *res_name, int option, int alphabet, int dimension, int order, int seq_num, \
				double lambda[order], double a[order][dimension][dimension], double null_lambda[order], double null_a[order][dimension][dimension] );


#endif /*_test_algs_h */
