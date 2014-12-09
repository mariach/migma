/* MTD program. Copyright (C) 2012. Created by Maria Chatzou and Pantelis Bagos.
   This program is free software: you can redistribute it and/or modify 
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.For more details see
   the GNU General Public License <http://www.gnu.org/licenses/>. */


#ifndef _train_algs_h
#define _train_algs_h


void initialize(char **trainseqArray, int seq_num, int alphabet, int diastasi, int order, int lamda_option, \
			 double a[order][diastasi][diastasi], double Lambda[order]);

void emAlg(char **trainseqArray, int seq_num, int alphabet, int diastasi, int order, double a[order][diastasi][diastasi], double lambda[order], \
                         double estim_lambda[order],double estim_a[order][diastasi][diastasi] );
			 
void viterbiAlg(char **trainseqArray, int seq_num, int alphabet, int diastasi, int order, double a[order][diastasi][diastasi], double lambda[order], \
                         double estim_lambda[order],double estim_a[order][diastasi][diastasi] );

void gradientAlg(char **trainseqArray, int seq_num, int alphabet, int diastasi, int order, double a[order][diastasi][diastasi], double lambda[order], \
                         double estim_lambda[order], double estim_a[order][diastasi][diastasi] );

void score(char **trainseqArray, int seq_num, int alphabet, int diastasi, int order, double Lambda[order],double a[order][diastasi][diastasi], double *Score);

void trainMTD(char **trainseqArray, int seq_num, char *model, int alphabet, int dimension, int order, int alg_option);



#endif /*_train_algs_h */

