/* MTD program. Copyright (C) 2012. Created by Maria Chatzou and Pantelis Bagos.
   This program is free software: you can redistribute it and/or modify 
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.For more details see
   the GNU General Public License <http://www.gnu.org/licenses/>. */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>	
#include <math.h>
#include <ctype.h>
#include <unistd.h>

#include "file_handle.c"
#include "predict_algs.c"
#include "train_algs.c"

#define PROT_LENGTH 20
#define DNA_LENGTH 4

#include "negative_models/uniprot_sprot.mtd_1.h"
#include "negative_models/uniprot_sprot.mtd_2.h"
#include "negative_models/uniprot_sprot.mtd_3.h"
#include "negative_models/uniprot_sprot.mtd_4.h"
#include "negative_models/uniprot_sprot.mtd_5.h"
#include "negative_models/uniprot_sprot.mtd_6.h"
#include "negative_models/uniprot_sprot.mtd_7.h"
#include "negative_models/uniprot_sprot.mtd_8.h"
#include "negative_models/uniprot_sprot.mtd_9.h"
#include "negative_models/uniprot_sprot.mtd_10.h"
#include "negative_models/uniprot_sprot.mtd_11.h"
#include "negative_models/uniprot_sprot.mtd_12.h"
#include "negative_models/uniprot_sprot.mtd_13.h"
#include "negative_models/uniprot_sprot.mtd_14.h"
#include "negative_models/uniprot_sprot.mtd_15.h"
#include "negative_models/uniprot_sprot.mtd_16.h"
#include "negative_models/uniprot_sprot.mtd_17.h"
#include "negative_models/uniprot_sprot.mtd_18.h"
#include "negative_models/uniprot_sprot.mtd_19.h"
#include "negative_models/uniprot_sprot.mtd_20.h"



#define PROT_LENGTH 20
#define DNA_LENGTH 4


int main ( int argc, char *argv[] )
{  
		int alphabet, order, dimension, alg=0;		
		int optFlag=0, runflag=0, i, I, J, option=0, seq_num=0;	
	        char *model_name, *file_name, *train_file_name, *out_file, *out_mtd_file;
	        char **testseqArray,  **trainseqArray;


		switch (argc)
    		{
    			 /**
    			  * Short reminder 
    			  */
    			case 1: 
      			//system("clear");
			fprintf(stderr,"MTD program. Copyright (C) 2012. \nCreated by Maria Chatzou and Pantelis Bagos.\n");
    			fprintf(stderr,"This program is free software: you can redistribute it and/or modify\n");
    			fprintf(stderr,"it under the terms of the GNU General Public License as published by\n");
    			fprintf(stderr,"the Free Software Foundation, either version 3 of the License, or\n");
    			fprintf(stderr,"(at your option) any later version.For more details see\n");
			fprintf(stderr,"the GNU General Public License <http://www.gnu.org/licenses/>.\n\n");

   
      			fprintf(stderr,"Usage:\n");
   			fprintf(stderr,"  (predict mode) %s <mtd_file> <input_file> <output_file> \n",argv[0]);
   			fprintf(stderr,"  (train mode) %s <input_training_file> <mtd_output_file>\n\t\t<alphabet> <algorithm> <order> \n",argv[0]);
   			fprintf(stderr,"  (train-predict mode) %s  <input_training_file> <mtd_output_file>\n\t\t       <alphabet> <algorithm> <mtd_order> <input_file> <output_file>\n\n",argv[0]);
   			
			fprintf(stderr,"  -m\t\t<mtd_file> the mtd model file.\n");
			fprintf(stderr,"  -i\t\t<input_file> the file containing the sequences,\n\t\t either DNA or protein. The file should be in FASTA format.\n");
			fprintf(stderr,"  -i_train\t<input_training_file> the file containing the sequences that\n\t\twill be used for training, either DNA or protein.The file should\n\t\tbe in FASTA format.\n");
			fprintf(stderr,"  -o\t\t<output_file> the output prediction file.\n");
			fprintf(stderr,"  -o_mtd\t<mtd_output_file> the mtd output file.\n");
			fprintf(stderr,"  -a\t\t<alphabet> the alphabet to be used\n\t\t( use: 1 for protein sequences | 2 for DNA sequences ).\n");
			fprintf(stderr,"  -or\t\t<mtd_order> the order of the mtd model.\n");
			fprintf(stderr,"  -alg\t\t<algorithm> (optional) the training algorithm\n\t\t(use: 1 for EM (default) | 2 for Viterbi | 3 for Gradient )\n");
      			fprintf(stderr,"  -p\t\twhen this flag is used an additional file will be outputed,\n\t\tcontaining the prediction together with the amino acid\n\t\tprobabilities.\n\n\n\n");    			
			exit(1);

			/**
      			 * Predict mode options
      			 */
    			case 7:
			for ( i = 1; i < argc; i++ )
			{
			      runflag=1;
			      if (!strcmp(argv[i], "-i")){
				file_name = argv[++i];
				optFlag++;
			      }
			      if (!strcmp(argv[i], "-m")){
				model_name = argv[++i];
				optFlag++;
			      }
			       if (!strcmp(argv[i], "-o")){
				out_file = argv[++i];
				optFlag++;
			      }
			}
			if(optFlag!=3){
			  fprintf(stderr,"\nSyntax Error: type %s for a short reminder.\n\n",argv[0]);
			  exit(1);
			}
      			break;

			case 8:
			for ( i = 1; i < argc; i++ )
			{
			      runflag=1;
			      if (!strcmp(argv[i], "-i")){
				file_name = argv[++i];
				optFlag++;
			      }
			      if (!strcmp(argv[i], "-m")){
				model_name = argv[++i];
				optFlag++;
			      }
			       if (!strcmp(argv[i], "-o")){
				out_file = argv[++i];
				optFlag++;
			      }
			      if (!strcmp(argv[i], "-p")){ 
				option=1;
				optFlag++;
			      }
			}
			if(optFlag!=4){
			  fprintf(stderr,"\nSyntax Error: type %s for a short reminder.\n\n",argv[0]);
			  exit(1);
			}  
      			break;
      			
      			/**
      			 * Train mode options
      			 */
      			case 9:
			for ( i = 1; i < argc; i++ )
			{
			      runflag=2;
			      if (!strcmp(argv[i], "-i_train")){
				train_file_name = argv[++i];
				optFlag++;
			      }
			      if (!strcmp(argv[i], "-o_mtd")){
				out_mtd_file = argv[++i];
				optFlag++;
			      }
			      if (!strcmp(argv[i], "-a")){
				alphabet = atoi(argv[++i]);
				optFlag++;
			      }
			      if (!strcmp(argv[i], "-or")){
				order = atoi(argv[++i]);
				optFlag++;
			      }
			}
			if(optFlag!=4){
			  fprintf(stderr,"\nSyntax Error: type %s for a short reminder.\n\n",argv[0]);
			  exit(1);
			}  
      			break;
      			
      			case 11:
			for ( i = 1; i < argc; i++ )
			{
			      runflag=2;
			      if (!strcmp(argv[i], "-i_train")){
				train_file_name = argv[++i];
				optFlag++;
			      }
			      if (!strcmp(argv[i], "-o_mtd")){
				out_mtd_file = argv[++i];
				optFlag++;
			      }
			      if (!strcmp(argv[i], "-a")){
				alphabet = atoi(argv[++i]);
				optFlag++;
			      }
			      if (!strcmp(argv[i], "-or")){
				order = atoi(argv[++i]);
				optFlag++;
			      }
			       if (!strcmp(argv[i], "-alg")){
				alg = atoi(argv[++i]);
				optFlag++;
			      }  
			}
			if(optFlag!=5){
			  fprintf(stderr,"\nSyntax Error: type %s for a short reminder.\n\n",argv[0]);
			  exit(1);
			}  
      			break;
      			
      			/**
      			 * Train-Predict mode options
      			 */
      			case 15:
			for ( i = 1; i < argc; i++ )
			{ 
			      runflag=3;
			      if (!strcmp(argv[i], "-i")){
				file_name = argv[++i];
				optFlag++;
			      }
			      if (!strcmp(argv[i], "-i_train")){
				train_file_name = argv[++i]; 
				optFlag++;
			      }
			       if (!strcmp(argv[i], "-o")){
				out_file = argv[++i];
				optFlag++;
			      }
			      if (!strcmp(argv[i], "-o_mtd")){
				out_mtd_file = argv[++i];
				model_name = out_mtd_file;
				optFlag++;
			      }
			      if (!strcmp(argv[i], "-a")){
				alphabet = atoi(argv[++i]);
				optFlag++;
			      }
			      if (!strcmp(argv[i], "-or")){
				order = atoi(argv[++i]);
				optFlag++;
			      }
			       if (!strcmp(argv[i], "-alg")){
				alg = atoi(argv[++i]);
				optFlag++;
			      }
			}
			if(optFlag!=7){
			  fprintf(stderr,"\nSyntax Error: type %s for a short reminder.\n\n",argv[0]);
			  exit(1);
			}  
      			break;
      
    			default:
      			fprintf(stderr,"\nSyntax Error: type %s for a short reminder.\n\n",argv[0]);
      			exit(1);
    		}

		printf("%d   \n",runflag);
		
		if( runflag==2 || runflag==3 )
		{
		    /**
		     * Train mode 
		     */
		    
		    if (fopen(train_file_name, "r") == NULL )
		    {
        		perror("\nError: failed to open <strings_file> for read");
          		exit(1);
		    }
			
		    if(alphabet==2){ dimension = DNA_LENGTH; }			
		    else if(alphabet==1){ dimension = PROT_LENGTH; }			
		    else{ 
			printf("\nError: wrong alphabet. Type either 1 for protein sequences or 2 for DNA sequences.\n\n");
			exit(1);
		    }
		
			
		    printf("\nThe MTD model will be saved at: %s \n", out_mtd_file);		
		
		    printf("Order of the MTD model: %d \n",order);	
	
		    if(alg==1){ printf("The algorithm used is: Expectation Maximization (EM)\n"); }
		    else if(alg==2){ printf("The algorithm used is: Viterbi\n"); }
		    else if(alg==3){ printf("The algorithm used is: Gradient\n"); }
		
		    seq_num=0;
		    trainseqArray=readFasta( train_file_name, alphabet, &seq_num);	 /* Read training file */			
		    printf("Number of sequences: %d \n\n", seq_num);
		
		    trainMTD(trainseqArray, seq_num, out_mtd_file, alphabet, dimension, order, alg);	/* Built MTD model */
		    
		}    
		
		if( runflag==1 || runflag==3 )
		{
		  
		  /**
		   * Predict mode 
		   */
		    printf("%s   \n",model_name);
		    if (fopen(model_name, "r") == NULL)  
		    { 
        		perror("\nError: failed to open <mtd_file> for reading!");
          		exit(1);
		    }
		    if (fopen(file_name, "r") == NULL)
		    {
        		perror("\nError: failed to open <strings_file> for reading!");
          		exit(1);
		    }
		
	
		    printf("\nThe results of the prediction are saved at: %s \n", out_file);
			  	  		  	  
		    find_Dimensios(model_name, &order, &I ,&J);		
		    printf("Order of the MTD model: %d \n",order);		
	
		    double lambda[order], a[order][I][J], null_lambda[order], null_a[order][I][J];
		    read_modelFile(model_name, order, I, J, lambda, a);	/* Read MTD model */
		
		    if(order<=20)
		    {
				int j,k;
				if(order==1){ 
				  
				 for(i=0; i<order; i++)
				   for(j=0; j<J; j++)
				     for(k=0; k<I; k++)
				       null_a[i][j][k]=a_1[i][j][k];
				  
				  for(i=0; i<order; i++) 
				    null_lambda[i] = lambda_1[i];
				} 
				if(order==2){ 
	
				 for(i=0; i<order; i++)
				   for(j=0; j<J; j++)
				     for(k=0; k<I; k++)
				       null_a[i][j][k]=a_2[i][j][k];
				  
				  for(i=0; i<order; i++) 
				    null_lambda[i] = lambda_2[i];
				}
				if(order==3){ 
				  printf("%lf\n\n", a_3[1][0][0]); 
				 for(i=0; i<order; i++){
				   for(j=0; j<J; j++){
				     for(k=0; k<I; k++){
				       null_a[i][j][k]=a_3[i][j][k];
				       printf("%.8lf ", null_a[i][j][k]);
				     }puts("");
				   }puts("");
				 }
				   for(i=0; i<order; i++) 
				    null_lambda[i] = lambda_3[i]; 
				}
				if(order==4){ 
	
				 for(i=0; i<order; i++)
				   for(j=0; j<J; j++)
				     for(k=0; k<I; k++)
				       null_a[i][j][k]=a_4[i][j][k];
				  
				  for(i=0; i<order; i++) 
				    null_lambda[i] = lambda_4[i];
				}
				if(order==5){ 
	
				 for(i=0; i<order; i++)
				   for(j=0; j<J; j++)
				     for(k=0; k<I; k++)
				       null_a[i][j][k]=a_5[i][j][k];
				  
				  for(i=0; i<order; i++) 
				    null_lambda[i] = lambda_5[i];
				}
				if(order==6){ 
	
				 for(i=0; i<order; i++)
				   for(j=0; j<J; j++)
				     for(k=0; k<I; k++)
				       null_a[i][j][k]=a_6[i][j][k];
				  
				  for(i=0; i<order; i++) 
				    null_lambda[i] = lambda_6[i];
				}
				if(order==7){ 
	
				 for(i=0; i<order; i++)
				   for(j=0; j<J; j++)
				     for(k=0; k<I; k++)
				       null_a[i][j][k]=a_7[i][j][k];
				  
				  for(i=0; i<order; i++) 
				    null_lambda[i] = lambda_7[i];
				}
				if(order==8){ 
	
				 for(i=0; i<order; i++)
				   for(j=0; j<J; j++)
				     for(k=0; k<I; k++)
				       null_a[i][j][k]=a_8[i][j][k];
				  
				  for(i=0; i<order; i++) 
				    null_lambda[i] = lambda_8[i];
				}
				if(order==9){ 
	
				 for(i=0; i<order; i++)
				   for(j=0; j<J; j++)
				     for(k=0; k<I; k++)
				       null_a[i][j][k]=a_9[i][j][k];
				  
				  for(i=0; i<order; i++) 
				    null_lambda[i] = lambda_9[i];
				}
				if(order==10){ 
	
				 for(i=0; i<order; i++)
				   for(j=0; j<J; j++)
				     for(k=0; k<I; k++)
				       null_a[i][j][k]=a_10[i][j][k];
				  
				  for(i=0; i<order; i++) 
				    null_lambda[i] = lambda_10[i];
				}
				if(order==11){ 
	
				 for(i=0; i<order; i++)
				   for(j=0; j<J; j++)
				     for(k=0; k<I; k++)
				       null_a[i][j][k]=a_11[i][j][k];
				  
				  for(i=0; i<order; i++) 
				    null_lambda[i] = lambda_11[i];
				}
				if(order==12){ 
	
				 for(i=0; i<order; i++)
				   for(j=0; j<J; j++)
				     for(k=0; k<I; k++)
				       null_a[i][j][k]=a_12[i][j][k];
				  
				  for(i=0; i<order; i++) 
				    null_lambda[i] = lambda_12[i];
				}
				if(order==13){ 
	
				 for(i=0; i<order; i++)
				   for(j=0; j<J; j++)
				     for(k=0; k<I; k++)
				       null_a[i][j][k]=a_13[i][j][k];
				  
				  for(i=0; i<order; i++) 
				    null_lambda[i] = lambda_13[i];
				}
				if(order==14){ 
	
				 for(i=0; i<order; i++)
				   for(j=0; j<J; j++)
				     for(k=0; k<I; k++)
				       null_a[i][j][k]=a_14[i][j][k];
				  
				  for(i=0; i<order; i++) 
				    null_lambda[i] = lambda_14[i];
				}
				if(order==15){ 
	
				 for(i=0; i<order; i++)
				   for(j=0; j<J; j++)
				     for(k=0; k<I; k++)
				       null_a[i][j][k]=a_15[i][j][k];
				  
				  for(i=0; i<order; i++) 
				    null_lambda[i] = lambda_15[i];
				}
				if(order==16){ 
	
				 for(i=0; i<order; i++)
				   for(j=0; j<J; j++)
				     for(k=0; k<I; k++)
				       null_a[i][j][k]=a_16[i][j][k];
				  
				  for(i=0; i<order; i++) 
				    null_lambda[i] = lambda_16[i];
				}
				if(order==17){ 
	
				 for(i=0; i<order; i++)
				   for(j=0; j<J; j++)
				     for(k=0; k<I; k++)
				       null_a[i][j][k]=a_17[i][j][k];
				  
				  for(i=0; i<order; i++) 
				    null_lambda[i] = lambda_17[i];
				}
				if(order==18){ 
	
				 for(i=0; i<order; i++)
				   for(j=0; j<J; j++)
				     for(k=0; k<I; k++)
				       null_a[i][j][k]=a_18[i][j][k];
				  
				  for(i=0; i<order; i++) 
				    null_lambda[i] = lambda_18[i];
				}
				if(order==19){ 
	
				 for(i=0; i<order; i++)
				   for(j=0; j<J; j++)
				     for(k=0; k<I; k++)
				       null_a[i][j][k]=a_19[i][j][k];
				  
				  for(i=0; i<order; i++) 
				    null_lambda[i] = lambda_19[i];
				}
				if(order==20){ 
	
				 for(i=0; i<order; i++)
				   for(j=0; j<J; j++)
				     for(k=0; k<I; k++)
				       null_a[i][j][k]=a_20[i][j][k];
				  
				  for(i=0; i<order; i++) 
				    null_lambda[i] = lambda_20[i];
				}
				  
				  for(i=0; i<order; i++){
				    null_lambda[i] = lambda_2[i]; printf("%lf ",null_lambda[i]);  }puts("");
				
				//char null_model_name[5000];     
		       		//sprintf(null_model_name,"/home/teo/Desktop/MTD_code_latest_version/negative_models/uniprot_sprot.mtd_%d",order);
		       		
		   		//read_modelFile(null_model_name, order, I, J, null_lambda, null_a);	/* Read null MTD model */	
		    }
		    		
	
		
		
		    /*int k;
		    for(k=0; k<order; k++){        printf("%lf ",lambda[k]);  }puts("");*/

	  	  
		    if(I==4 && J==4){ /* Check if the alphabet is definied correctly */
			alphabet=2;
			dimension = DNA_LENGTH; 
			printf("Alphabet used: \tDNA (A, C, G, T)\n");
		    }			
		    else if(I==20 && J==20){ 
			alphabet=1; 
			dimension = PROT_LENGTH;
			printf("Alphabet used: \tPROTEIN (A, R, N, D, C, E, Q, G, H, I, L, K, M, F, P, S, T, W, Y, V)\n");
		    }				
		    else{  
			perror("\nError: wrong MTD model used. Choose another MTD file or built a new MTD model.");
			exit(1);
		    }

		    seq_num=0;
		    testseqArray=readFasta( file_name, alphabet, &seq_num);   /* Read test file */

		    printf("Number of sequences: %d \n\n", seq_num);

		
		    predict_fastaSet(file_name, testseqArray, out_file, option, alphabet, dimension, order, seq_num, lambda, a, null_lambda, null_a); /* Predict */
		
		}

		return(0);
}

