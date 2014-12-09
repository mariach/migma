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

#include "file_handle.c"
#include "train_algs.c"

#define PROT_LENGTH 20
#define DNA_LENGTH 4


int main ( int argc, char *argv[] )
{  

		int alphabet, order, dimension, algorithm=1, seq_num=0;
		char *file_name, **trainseqArray, model[200];


		switch (argc)
    		{
    			case 1:  /* short reminder */
      			//system("clear");
			fprintf(stderr,"\nMTD program. Copyright (C) 2012. \nCreated by Maria Chatzou and Pantelis Bagos.\n");
    			fprintf(stderr,"This program is free software: you can redistribute it and/or modify\n");
    			fprintf(stderr,"it under the terms of the GNU General Public License as published by\n");
    			fprintf(stderr,"the Free Software Foundation, either version 3 of the License, or\n");
    			fprintf(stderr,"(at your option) any later version.For more details see\n");
			fprintf(stderr,"the GNU General Public License <http://www.gnu.org/licenses/>.\n\n");
  
      			fprintf(stderr,"Usage:\n");
   			fprintf(stderr,"  %s <strings_file> <alphabet> <order> <algorithm>\n\n",argv[0]);
			fprintf(stderr,"<strings_file>\t\t file containing the sequences, either DNA or protein, that are going to be used for training.<strings_file> should be in FASTA format.\n");
			fprintf(stderr,"<alphabet>\t\t should be either 1 for protein sequences or 2 for DNA sequences and <order> is the order of the MTD model to be built.\n");   			
			fprintf(stderr,"<algorithm>\t\t should be either 1 for Expectation Maximization, 2 for Viterbi or 3 for Gradient (not necessary to define/ default Expectation Maximization).\n\n\n");  
			exit(1);


    			case 4:
      			file_name=argv[1]; 
			alphabet=atoi(argv[2]);
			order=atoi(argv[3]);			
      			break;

			case 5:
      			file_name=argv[1]; 
			alphabet=atoi(argv[2]);
			order=atoi(argv[3]);
			algorithm=atoi(argv[4]);			
			  if(algorithm<1 || algorithm>3){ 
			    printf("\nError: wrong choise of parameter estimation algorithm. Type either 1 for Expectation Maximization, 2 for Viterbi or 3 for Gradient.\n\n");
			    exit(1);
			  }
      			break;

    			default:
      			fprintf(stderr,"\nSyntax Error: type %s for a short reminder\n\n",argv[0]);
			exit(1);
		}

		if (fopen(file_name, "r") == NULL)
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
		
		
		sprintf(model,"%s.mtd",file_name);	
		printf("\nThe MTD model will be saved at: %s \n",model);		
		
		printf("Order of the MTD model: %d \n",order);	
	
		if(algorithm==1){ printf("The algorithm used is: Expectation Maximization (EM)\n"); }
		else if(algorithm==2){ printf("The algorithm used is: Viterbi\n"); }
		else if(algorithm==3){ printf("The algorithm used is: Gradient\n"); }
		
		trainseqArray=readFasta( file_name, alphabet, &seq_num);	 /* Read training file */			
		printf("Number of sequences: %d \n\n", seq_num);
		
		trainMTD(trainseqArray, seq_num, model, alphabet, dimension, order, algorithm);	/* Built MTD model */
		return(0);
}




