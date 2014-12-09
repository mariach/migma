/* MTD program. Copyright (C) 2012. Created by Maria Chatzou and Pantelis Bagos.
   This program is free software: you can redistribute it and/or modify 
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version. For more details see
   the GNU General Public License <http://www.gnu.org/licenses/>. */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>	
#include <math.h>
#include <ctype.h>

#include "predict_algs.h"


void score_pred(char **testseqArray, int alphabet, int dimension, int order, int seq_num, double Lambda[order],double a[order][dimension][dimension], double Score[seq_num]){
		
	int SIDi=0, SIDj=0, i, k, flag=0, seqcount=0;
	double Score1=0.0, Score2=0.0;			
	char zebgos[3];

	FILE *tmpfile = fopen("tmp","w");
			
	while(*testseqArray!='\0')
	{
	   if(**testseqArray=='>' ){  fprintf(tmpfile,"%s\n",*testseqArray);  }
	   else
	   {
		Score2=0;
		for(i=1; i < strlen(*testseqArray); i++)
		{		
			if(i<order){
				Score1=0;
				for(k=0; k<i; k++){		
							
					sprintf(zebgos,"%c%c", (*testseqArray)[i-k-1], (*testseqArray)[i]); 																							
					abc(&SIDi,&SIDj,alphabet,zebgos);             	    
								 					 
					Score1=Score1+Lambda[k]*a[k][SIDi][SIDj]; 									  													 					  										
				}
		     	}
		     	else{	
		     		Score1=0;																				
				for(k=0; k<order; k++){		
							
					sprintf(zebgos,"%c%c", (*testseqArray)[i-k-1], (*testseqArray)[i]); 																									
					abc(&SIDi,&SIDj,alphabet,zebgos);             	    
								 					 
					Score1=Score1+Lambda[k]*a[k][SIDi][SIDj];   										  													 					  										
				}
			}
			
			
			if(Score1< 0.0000001){ 
				Score2=(Score2)-4; 
				fprintf(tmpfile,"%c\t0\t-4\n",(*testseqArray)[i]);
			}
			else{	
				Score2= Score2 +( log(Score1) ); 
				fprintf(tmpfile,"%c\t%lf\t%lf\n",(*testseqArray)[i],Score1,log(Score1));
			} 	   	   	   	 																							  
			flag=1;   
		}

		if(flag==1){
			Score[seqcount]=Score2;	  																			
			flag=0;
		} 			
		
		seqcount++;
	   }
	   testseqArray++;	 		
	}

	fclose(tmpfile);	
}


void null_model_score(char **testseqArray, int alphabet, int dimension, int order, int seq_num, double Lambda[order],double a[order][dimension][dimension],double null_Lambda[order],double null_a[order][dimension][dimension], double Score[seq_num], double null_Score[seq_num]){
		
	int SIDi=0, SIDj=0, i, k, seqcount=0;
	double Score1=0.0, Score2=0.0, Score1_null=0.0, Score2_null=0.0;				
	char zebgos[3];

	FILE *tmpfile = fopen("tmp","w");

			
	while(*testseqArray!='\0')
	{
	   if(**testseqArray=='>' ){  fprintf(tmpfile,"%s\n",*testseqArray);  }
	   else
	   {
		Score2=0;
		Score2_null=0;
		for(i=1; i < strlen(*testseqArray); i++)
		{		
			if(i<order){
				Score1=0;
				Score1_null=0;
				for(k=0; k<i; k++){		
							
					sprintf(zebgos,"%c%c", (*testseqArray)[i-k-1], (*testseqArray)[i]); 																							
					abc(&SIDi,&SIDj,alphabet,zebgos);             	    
								 					 
					Score1=Score1+Lambda[k]*a[k][SIDi][SIDj]; 
					Score1_null=Score1_null+null_Lambda[k]*null_a[k][SIDi][SIDj]; 
				}
		     	}
		     	else{	
		     		Score1=0;
				Score1_null=0;
				for(k=0; k<order; k++){		
							
					sprintf(zebgos,"%c%c", (*testseqArray)[i-k-1], (*testseqArray)[i]); 																									
					abc(&SIDi,&SIDj,alphabet,zebgos);             	    
								 					 
					Score1=Score1+Lambda[k]*a[k][SIDi][SIDj]; 
					Score1_null=Score1_null+null_Lambda[k]*null_a[k][SIDi][SIDj]; //printf("%lf %lf\n", Score1_null, null_Lambda[k] );
				}
			}
			
			 
			//null model score
			if(Score1_null< 0.0000001)
				Score2_null=(Score2_null)-4; 
			else 
				Score2_null= Score2_null +( log(Score1_null) ); 	  
			
			//positive model score
			if(Score1< 0.0000001){ 
				Score2=(Score2)-4; 
				fprintf(tmpfile,"%c\t0\t-4\n",(*testseqArray)[i]);
			}
			else{	
				Score2= Score2 +( log(Score1) ); 
				fprintf(tmpfile,"%c\t%lf\t%lf\n",(*testseqArray)[i],Score1,log(Score1));
			} 	 	   	   	   	 					  
		}

		Score[seqcount]=Score2;	 
		null_Score[seqcount]=Score2_null;			
		
		seqcount++;
	   }
	   testseqArray++;	 		
	}

	fclose(tmpfile);	
}



void null_Score( char **testseqArray, int alphabet, int dimension, int seq_num, double null_Score[seq_num]){
	

	double Score;
  	int seq_count=0, i;
	char *pos;

	char amino[]="ACDEFGHIKLMNPQRSTVWYBZ";
	double aapos[22]={ 0.0817, 0.014 , 0.0542 , 0.0674 , 0.0388 , 0.0704 , 0.0228 , 0.0594 , 0.0587 , 0.0967 , 0.0241 , 0.0407 , 0.0474 , 0.0395 , 0.055 , 0.0662 , 0.0534 , 0.0682 , 0.0109 , 0.0293, 0.0542, 0.0704  }; 

	char nucleo[]="AaCcGgTt";
	double nupos[4]={ 0.25, 0.25, 0.25, 0.25 }; 
				

	while(*testseqArray!='\0')
	{
	   if(**testseqArray!='>')
	   {
		Score=0;
		for(i=0; i < strlen(*testseqArray); i++)
		{		
			if(alphabet==1){
				pos = strchr( amino, ((*testseqArray)[i]) );
				if( pos ){  Score = Score + log(aapos[pos-amino]); }
			}

			if(alphabet==2){
				pos = strchr( nucleo, ((*testseqArray)[i]) );
				if( pos ){  
					if( (pos-nucleo)%2==0 ){ Score = Score + log(nupos[(pos-nucleo)/2]); }
					else{ Score = Score + log(nupos[(pos-1-nucleo)/2]); }
				}
			}		  	  	  	  	 		
		}
		null_Score[seq_count]=Score;				
		seq_count++;
	   }
	   testseqArray++;	 		
	}
}




void predict_fastaSet(char *file_name, char **testseqArray, char *res_name, int option, int alphabet, int dimension, int order, int seq_num, double lambda[order], double a[order][dimension][dimension], double null_lambda[order], double null_a[order][dimension][dimension] ){
		
	int seq_count=0;
	double Score[seq_num], nullScore[seq_num];
	
	FILE *res_file=fopen(res_name,"w");	
	
	/**
	 * If order is smaller than 20 it uses the pre-computed null models to calculate the probability of a sequence belonging to the null model
	 * else it does a normal probability calculation
	 */
	if(order>20)
	{
	    score_pred(testseqArray, alphabet, dimension, order, seq_num, lambda, a, Score );		
	    null_Score(testseqArray, alphabet, dimension, seq_num, nullScore);  
	}
	else
	    null_model_score(testseqArray, alphabet, dimension, order, seq_num, lambda, a, null_lambda, null_a, Score, nullScore);

		
		 
	/*_______________________________PRINT RESULTS_______________________________________*/

	fprintf(res_file,"\tName\t\tScore\t\t\tNull_Score\t\t\t\tTotal_Score\n");
	fprintf(res_file,"_________________________________________________________________________________________________________\n\n");	

		
	while(*testseqArray!='\0')
	{
	   if(**testseqArray=='>')
	   {	
		fprintf( res_file,"%s\t%lf\t\t%lf\t\t\t%lf\n", *testseqArray, Score[seq_count], nullScore[seq_count], (Score[seq_count]-nullScore[seq_count]) );																
		seq_count++;
	   }
	   testseqArray++;	
	}			

		
	if(option==0){  remove("tmp");  }
	else{ 
		char aaProb[100];
		sprintf(aaProb,"%s.aaProbabilities",res_name);
		rename("tmp",aaProb); 
	}

	fclose(res_file);		
}
