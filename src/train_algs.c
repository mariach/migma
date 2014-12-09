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
#include <omp.h>

#include "file_handle.h"
#include "train_algs.h"




void initialize( char **trainseqArray,  int seq_num, int alphabet, int diastasi, int order, int lamda_option, double a[order][diastasi][diastasi], double Lambda[order] ){
	
	double A[order][diastasi][diastasi], sum=0;
	int i, j, k;
	
		
	for(i=0; i<order; i++){	/* Initialize A[i][j][k], estim_a[i][j][k] = 0.00001 */
		for(j=0; j<diastasi; j++){
			for(k=0; k<diastasi; k++){
				A[i][j][k]=0.000000000000001;
				a[i][j][k]=0.000000000000001;				
	} } }	 	 
		
	   	   
	if(lamda_option==1){	/* Initialize Lambda */
		for(i=0; i<order; i++){  Lambda[i]=(1/(double)order);  }
	}
	else if(lamda_option==2){		
		double par=0;
		for(i=0; i<order; i++){  par=par+i+1;  }
		for(i=0; i<order; i++){  Lambda[i]=((double)(order-i)/par);  }
	} 
	else{	printf("Error: wrong choise of Lambda. Write 1:Lambda=(1/order)  2:Lambda=((order-i)/par)  ");	}
	
	
	int chunk=1;
	/**
	 * Create A array
	 */
	#pragma omp parallel for private(i,j,k) shared(trainseqArray,A) //schedule(dynamic,chunk)
	for (i = 1; i < seq_num*2; i+=2){
	    for (j = order; j < strlen(trainseqArray[i]); j++){
	      for (k = 0; k < order; k++)
	      {
		int SIDi=0, SIDj=0;
		char zebgos[3];	
		
		sprintf(zebgos,"%c%c", trainseqArray[i][j-k-1], trainseqArray[i][j]);
		abc(&SIDi,&SIDj,alphabet,zebgos); 

		#pragma omp critical 
		A[k][SIDi][SIDj]+=1;
		
		//printf("i=%d, j=%d, k=%d, thread = %d, pair=%s \n", i, j, k, omp_get_thread_num(),zebgos);
		//printf("%s %d %d %lf\n",zebgos, SIDi,SIDj,A[k][SIDi][SIDj] );
	      }
	    }
	} 
    
	
	for(i=0; i<order; i++){  	/* Create a array */ 
		for(j=0; j<diastasi; j++){
			for(k=0; k<diastasi; k++){  sum=sum+ A[i][j][k];  } 
			for(k=0; k<diastasi; k++){ //printf("%.0lf\t",A[i][j][k]);
			  if(sum<0.000001)
			      a[i][j][k] = 0.000000000000001; 
			  else
			      a[i][j][k] = (double)A[i][j][k]/(double)sum;  
			}
			sum=0;	      	   	    
		}
	}	
}





void emAlg(char **trainseqArray, int seq_num, int alphabet, int diastasi, int order, double a[order][diastasi][diastasi], double lambda[order], double Lambda[order], double estim_a[order][diastasi][diastasi] ){
	
	char zebgos[3];
	double A[order][diastasi][diastasi], em_sum=0, p=0, em=0, lamda_sum=0, sum=0;
	int i, j, k, SIDi=-1, SIDj=-1;
	

		
	for(i=0; i<order; i++){           /* Initialize Lambda, A[i][j][k], estim_a[i][j][k] = 0.00001 */
		Lambda[i]=0.0000000000001;
		for(j=0; j<diastasi; j++){
			for(k=0; k<diastasi; k++){	  	  	  	  	  
				A[i][j][k]=0.0000000000001;
				estim_a[i][j][k]=0.0000000000001;
	}}}  	   	 

	#pragma omp parallel for private(i,j,k,zebgos,SIDi,SIDj) shared(trainseqArray,a,Lambda) reduction(+:em_sum, p, em)
	for (j = 1; j < seq_num*2; j+=2){
		for (i = order; i < strlen(trainseqArray[j]); i++)
		{
			em_sum=0;
	    			
			for(k=0; k<order; k++)
			{													
				sprintf(zebgos,"%c%c",trainseqArray[j][i-k-1], trainseqArray[j][i]); 
				abc(&SIDi,&SIDj,alphabet,zebgos);           		 			 	 	 	 	 	 	 	 
									 				
				p=a[k][SIDi][SIDj]*lambda[k]; 									 	
				em_sum+=p;			  
			}
			
			for(k=0; k<order; k++)
			{	   	   	   	   	   	   	   	
				sprintf(zebgos,"%c%c",trainseqArray[j][i-k-1], trainseqArray[j][i]);   
				abc(&SIDi,&SIDj,alphabet,zebgos);           	 	 	 	 	 	 	 	 	 									 					 	  	  	  	  	  	  	  
									 												  										
				p=a[k][SIDi][SIDj]*lambda[k];  			
				em=p/em_sum;
				
				#pragma omp critical
				A[k][SIDi][SIDj]+=em;	
				
				#pragma omp critical
				Lambda[k]+=em;	 
				
				em=0;
			}
		}																																								
	}	
		 	 	

			
	for(i=0; i<order; i++){  lamda_sum=lamda_sum+ Lambda[i];  }	/* Create Lambda array */				
	for(i=0; i<order; i++){  Lambda[i]=(Lambda[i]/lamda_sum);  } 

						
	for(i=0; i<order; i++){	/* Create a array */																				 	    
		for(j=0; j<diastasi; j++){
			for(k=0; k<diastasi; k++){  sum=sum+ A[i][j][k];  } 
			for(k=0; k<diastasi; k++){													
				if(sum<0.000001)
				    estim_a[i][j][k] = 0.000000000000001;
				else
				    estim_a[i][j][k] = (double)A[i][j][k]/(double)sum;    	 	 	 	 	 	   	   	   	   	   	   	   	   	   
			} 	  	  	  	  	  	  															 
			sum=0;										      	   	    															 
		}			   																											 	 	   
	}
}


void viterbiAlg(char **trainseqArray, int seq_num, int alphabet, int diastasi, int order, double a[order][diastasi][diastasi], double lambda[order], double Lambda[order], double estim_a[order][diastasi][diastasi] ){
	
	char zebgos[3];
	double A[order][diastasi][diastasi], p=0.0, lamda_sum=0.0, sum=0.0,  max, eqp[order];
	int i, j, k, SIDi=-1, SIDj=-1, thesi, eqcount, eqSIDiA[order], eqSIDjA[order], equalthesi[order], SIDiA=-1, SIDjA=-1;
	

		
	for(i=0; i<order; i++){           /* Initialize Lambda, A[i][j][k], estim_a[i][j][k] = 0.00001 */
		Lambda[i]=0.0000000000001; 					//printf("lag %d:\n",i);
		for(j=0; j<diastasi; j++){
			for(k=0; k<diastasi; k++){	  	  	  	  	  
				A[i][j][k]=0.0000000000001;
				estim_a[i][j][k]=0.000000000000001;		//printf("%lf ",a[i][j][k]);
	}}}  	   	 

	
	for (j = 1; j < seq_num*2; j+=2){
		for (i = 1; i < strlen(trainseqArray[j]); i++)
		{
			max=0.0; thesi=0; eqcount=0; 					
				
			if(i<order){
				for(k=0; k<i; k++)
				{					
					sprintf(zebgos,"%c%c",trainseqArray[j][i-k-1], trainseqArray[j][i]); 
					abc(&SIDi,&SIDj,alphabet,zebgos);           		 			 	 	 	 	 	 	 	 
									 				
					p=a[k][SIDi][SIDj]*lambda[k]; 		//printf("%s %lf  %lf\n",zebgos, a[k][SIDi][SIDj],lambda[k] );							 
					
					if(p>max ){
					      max=p; 
					      thesi=k; 	
					      SIDiA= SIDi;
					      SIDjA= SIDj;
					}
					else if(p==max && p>0.0){ 		//printf("E: %s %lf %lf\n",zebgos, max, p);
					      equalthesi[eqcount]=k;
					      eqp[eqcount]=p;
					      eqSIDiA[eqcount]= SIDi;
					      eqSIDjA[eqcount]= SIDj;
					      eqcount++;
					} 									  													 					  										
				}
		     	}
		     	else{	
				for(k=0; k<order; k++)
				{													
					sprintf(zebgos,"%c%c",trainseqArray[j][i-k-1], trainseqArray[j][i]); 
					abc(&SIDi,&SIDj,alphabet,zebgos);           		 			 	 	 	 	 	 	 	 
									 				
					p=a[k][SIDi][SIDj]*lambda[k]; 		//printf("%s %lf  %lf\n",zebgos, a[k][SIDi][SIDj],lambda[k] );							 
					
					if(p>max ){
					      max=p; 
					      thesi=k; 	
					      SIDiA= SIDi;
					      SIDjA= SIDj;
					}
					else if(p==max && p>0.0){  		//printf("E: %s %lf %lf\n",zebgos, max, p);
					      equalthesi[eqcount]=k;
					      eqp[eqcount]=p;
					      eqSIDiA[eqcount]= SIDi;
					      eqSIDjA[eqcount]= SIDj;
					      eqcount++;
					}
				}
				
				//A[thesi][SIDiA][SIDjA]= A[thesi][SIDiA][SIDjA]+1;	 		//printf(" %lf \n",A[thesi][SIDiA][SIDjA]);								  	  	  	  	  	  	  	  	  	  	  	  	  	  
				//Lambda[thesi] = Lambda[thesi]+1;
			}
			
			 	 	 	
			A[thesi][SIDiA][SIDjA]= A[thesi][SIDiA][SIDjA]+1;	 		//printf(" %lf \n",A[thesi][SIDiA][SIDjA]);								  	  	  	  	  	  	  	  	  	  	  	  	  	  
			Lambda[thesi] = Lambda[thesi]+1;
			
			if(eqcount!=0){ 
			  for(k=0; k<eqcount; k++){ 
				if( max == eqp[k] ){
					A[equalthesi[k]][eqSIDiA[k]][eqSIDjA[k]]= A[equalthesi[k]][eqSIDiA[k]][eqSIDjA[k]]+1;	 		//printf(" %lf \n",A[thesi][SIDiA][SIDjA]);																	  	  	  	  	  	  	  	  	  	  	  	  	  	  
			   		Lambda[equalthesi[k]] = Lambda[equalthesi[k]]+1;       							//printf("e %d \n", equalthesi[k]); 
				}
				equalthesi[k]=0;
				eqp[k]=0;
			  }
			}
		}																																								
	}	
		 	 	
			
	for(i=0; i<order; i++){  lamda_sum=lamda_sum+ Lambda[i];  }	/* Create Lambda array */				
	for(i=0; i<order; i++){  Lambda[i]=(Lambda[i]/lamda_sum);  } 

	puts("");					
	for(i=0; i<order; i++){	/* Create a array */		//printf("lag %d:\n",i);																				 	    
		for(j=0; j<diastasi; j++){
			for(k=0; k<diastasi; k++){  sum=sum+ A[i][j][k];  } 
			for(k=0; k<diastasi; k++){													
				if(sum<0.00001)
				    estim_a[i][j][k] = 0.000000000000001;
				else
				    estim_a[i][j][k] = (double)A[i][j][k]/(double)sum;    
				//printf("%lf| %lf %lf - ", (double) estim_a[i][j][k], A[i][j][k], sum);
			}//puts("");	  															 
			sum=0.0;										      	   	    															 
		}//puts("");			   																											 	 	   
	}
}



void gradientAlg(char **trainseqArray, int seq_num, int alphabet, int diastasi, int order, double a[order][diastasi][diastasi], double lambda[order], double estim_Lambda[order], double estim_a[order][diastasi][diastasi] ){
	
	char zebgos[3];
	int i, j, k, SIDi, SIDj;
	double A[order][diastasi][diastasi], Z[order][diastasi][diastasi], sumA[order][diastasi], em_sum, n=0.00000001, nl=0.1, p=0, em=0, lamda_sum=0, sum=0;
	
		
	for(i=0; i<order; i++){           /* Initialize Lambda, A[i][j][k], estim_a[i][j][k] = 0.00001 */
		//Lambda[i]=0.00001;
		for(j=0; j<diastasi; j++){
			for(k=0; k<diastasi; k++){
				Z[i][j][k]=0.00001;	  	  	  	  
				A[i][j][k]=0.00001;
				estim_a[i][j][k]=0.00001;
	}}}  	
	
	int t=0;
	for (j = 1; j < seq_num*2; j+=2){
		for (i = order; i < strlen(trainseqArray[j]); i++)
		{
			em_sum=0;
						
			for(k=0; k<order; k++){													
				sprintf(zebgos,"%c%c",trainseqArray[j][i-k-1], trainseqArray[j][i]);     //printf("%s: \t",zebgos);
				abc(&SIDi,&SIDj,alphabet,zebgos);           		 		 //printf(" %d %d\t",SIDi,SIDj);	 	 	 	 	 	 	 	 
									 				
				p=a[k][SIDi][SIDj]*lambda[k]; 						 //printf("a: %lf\t l: %lf\t p: %lf\n",a[k][SIDi][SIDj],lambda[k],p);			 													
			     						  
				em_sum=em_sum+p;   							//printf("EM=%lf\n", em_sum);																									  			  
			}										//printf(" max=%lf thesi=%d SIDiA=%d SIDjA=%d \t",max,thesi,SIDiA, SIDjA);
											
			for(k=0; k<order; k++){	   	   	   	   	   	   	   	 // printf("YPOLOGISMOS PI8ANOTHTAS LETTER---EM=%lf\n", em_sum);
				sprintf(zebgos,"%c%c",trainseqArray[j][i-k-1], trainseqArray[j][i]);   //   printf("%s: \t",zebgos);
				abc(&SIDi,&SIDj,alphabet,zebgos);           	 //  printf(" %d %d\n",SIDi,SIDj);	 	 	 	 	 	 	 	 									 					 	  	  	  	  	  	  	  

									 												  										
				p=a[k][SIDi][SIDj]*lambda[k];  						  			 // printf("!!!a: %lf\t l: %lf\t p: %lf\n",a[k][SIDi][SIDj],lambda[k],p);
				em=p/em_sum;
			       								//  printf("!!!em=%lf\t",em);					 
				A[k][SIDi][SIDj]= A[k][SIDi][SIDj]+em;	           //SWsto Σ(αλ/Σαλ) printf("A[%d][%d][%d]: %lf \t",k, SIDi, SIDj, A[k][SIDi][SIDj]);
				em=0;
			}
		}																																						
	}

	for(i=0; i<order; i++){	/* Create sumA array */															 	    
		for(j=0; j<diastasi; j++){
			for(k=0; k<diastasi; k++){  sumA[i][j]=sumA[i][j]+A[i][j][k]; }       	   	    															 
		}																											 	 	   
	}
	t=0; 
	while(trainseqArray[t])
	{  
	   if(trainseqArray[t][0]!='>')
	   { 
		for(i=order; i < strlen(trainseqArray[t]); i++)
		{										
			for(k=0; k<order; k++){	   	   	   	   	   	   	   	 // printf("YPOLOGISMOS PI8ANOTHTAS LETTER---EM=%lf\n", em_sum);
				sprintf(zebgos,"%c%c",trainseqArray[t][i-k-1], trainseqArray[t][i]);     //   printf("%s: \t",zebgos);
				abc(&SIDi,&SIDj,alphabet,zebgos);           				 //  printf(" %d %d\n",SIDi,SIDj);
																 
				Z[k][SIDi][SIDj]= Z[k][SIDi][SIDj] + (A[k][SIDi][SIDj]-a[k][SIDi][SIDj]*sumA[k][SIDi]);	 //SWsto aλ/Σαλ-a*Σ(αλ/Σαλ)     printf("A[%d][%d][%d]: %lf \t",k, SIDi, SIDj, A[k][SIDi][SIDj]);
				estim_Lambda[k]= estim_Lambda[k]+ (A[k][SIDi][SIDj]-lambda[k]);	 	 	 	            		 //printf("  lamda[%d]: %lf \n",k,Lambda[k]);	   
				em=0; 
			}
			
		}
	   }
	   t++;																																						
	}
		 	 	

       /*_________________Creation of Lambda_____________________*/		
	for(i=0; i<order; i++){  lamda_sum=lamda_sum + lambda[i]*(nl*estim_Lambda[i]);  }					
	for(i=0; i<order; i++){  estim_Lambda[i]=( (lambda[i]*(nl*estim_Lambda[i]))/lamda_sum);  } 

						
	for(i=0; i<order; i++){	/* Create a array */															 	    
		for(j=0; j<diastasi; j++){
			for(k=0; k<diastasi; k++){  sum=sum+ a[i][j][k]*exp(n*Z[i][j][k]);  } 
			for(k=0; k<diastasi; k++){													
				if(sum<0.000000001){ estim_a[i][j][k] = 0.0001; }
				else{ estim_a[i][j][k] = ( (double) a[i][j][k]*exp(n*Z[i][j][k]) )/(double) sum; } 	//printf("%lf ",estim_a[i][j][k]);  	 	 	 	 	 	   	   	   	   	   	   	   	   	   
			} 	  	  	  	  	  	  															 
			sum=0;										//printf("\n");      	   	    															 
		}											//printf("\n");																				 	 	   
	}
	
}



void score( char **trainseqArray, int seq_num, int alphabet, int diastasi, int order, double Lambda[order],double a[order][diastasi][diastasi], double *Score ){
		
	double Score1=0.0, Score2=0.0, Score3=0.0;
	int i, j, k;
			
	#pragma omp parallel for private(i,j,k) shared(trainseqArray,a,Lambda) reduction(+:Score1,Score2,Score3)
	for (j = 1; j < seq_num*2; j+=2){
		Score2=0;
		for (i = 1; i < strlen(trainseqArray[j]); i++)
		{
			Score1=0.0;
			if(i<order){
			     for(k=0; k<i; k++){
				int SIDi=-1, SIDj=-1;
				char zebgos[3];	
								
				sprintf(zebgos,"%c%c",trainseqArray[j][i-k-1], trainseqArray[j][i]);
				abc(&SIDi,&SIDj,alphabet,zebgos);           		 	 	 	 	 	 	 	 
				 
				Score1=Score1+Lambda[k]*a[k][SIDi][SIDj]; 
			    }
			}
			else{
			    for(k=0; k<order; k++){
				int SIDi=-1, SIDj=-1;
				char zebgos[3];	
								
				sprintf(zebgos,"%c%c",trainseqArray[j][i-k-1], trainseqArray[j][i]);
				abc(&SIDi,&SIDj,alphabet,zebgos);           		 	 	 	 	 	 	 	 
				 
				Score1=Score1+Lambda[k]*a[k][SIDi][SIDj]; 
			    }
			}											  													  						  										
			
			if(Score1< 0.0001)
			      Score2=(Score2)-4;
			else
			      Score2=(Score2)+(log(Score1));
		}
		
		Score3=Score3+Score2;
	} 
	(*Score)=Score3;		
}





void trainMTD( char **trainseqArray, int seq_num, char *model, int alphabet, int dimension, int order, int alg_option ){		
		
	FILE *mtd_model = fopen(model,"w");
		
	double Score=0, Old_Score=0;
	double init_a[order][dimension][dimension], a[order][dimension][dimension], estim_a[order][dimension][dimension];
	double Lambda[order], lambda[order], estim_Lambda[order];

	int i,j,k,flag=2;
		
	

	initialize(trainseqArray, seq_num, alphabet, dimension, order, 1, init_a, lambda); 	/* 1:lamda=(1/order)  2:lamda=((order-i)/par) */
	score(trainseqArray, seq_num, alphabet, dimension, order, lambda, init_a, &Old_Score ); 	
	printf("Initial : %lf\n", Old_Score); 	

	//find_lambda(trainseqArray, alphabet, dimension, order, init_a, lambda, Lambda, a);
	if(alg_option==1){ emAlg(trainseqArray, seq_num, alphabet, dimension, order, init_a, lambda, Lambda, a); }	/* 1:EM 2:Viterbi */
	else if(alg_option==2){ viterbiAlg(trainseqArray, seq_num, alphabet, dimension, order, init_a, lambda, Lambda, a); }
	else if(alg_option==3){ gradientAlg(trainseqArray, seq_num, alphabet, dimension, order, init_a, lambda, Lambda, a);  }   	/* 3:Gradient */
	else{ emAlg(trainseqArray, seq_num, alphabet, dimension, order, init_a, lambda, Lambda, a); }				/* default:EM */
	
	score(trainseqArray, seq_num, alphabet, dimension, order, Lambda, a, &Score );  	 	 	 
	printf("Iteration 1: %lf\n",Score);

	
	while( (Score-Old_Score)>0.0001 )		/* Iterations */
	{			 
		Old_Score = Score; 				 							
		Score=0;
				 
		for(i=0; i<order; i++){	
			 estim_Lambda[i]=0;	 
			 estim_Lambda[i]= Lambda[i];			
			 Lambda[i]=0;
		}									
					
		for(i=0; i<order; i++){	   	   	   	   	   					
			for(j=0; j<dimension; j++){
				for(k=0; k<dimension; k++){
					estim_a[i][j][k]=0;
					estim_a[i][j][k]= a[i][j][k]; 	
					a[i][j][k]=0;										  
				}				      														  
			}															  
		}
		
		if(alg_option==1){            emAlg( trainseqArray, seq_num, alphabet, dimension, order, estim_a, estim_Lambda, Lambda, a ); }	/* 1:EM 2:Viterbi */
	        else if(alg_option==2){  viterbiAlg( trainseqArray, seq_num, alphabet, dimension, order, estim_a, estim_Lambda, Lambda, a ); }
		else if(alg_option==3){ gradientAlg( trainseqArray, seq_num, alphabet, dimension, order, estim_a, estim_Lambda, Lambda, a ); }   		/* 3:Gradient */
		else{ 			      emAlg( trainseqArray, seq_num, alphabet, dimension, order, estim_a, estim_Lambda, Lambda, a ); }				/* default:EM */
		
		score(trainseqArray, seq_num, alphabet, dimension, order, Lambda, a, &Score ); 	 	 	
		printf("Iteration %d: %lf\n",flag, Old_Score); 	 	 	    
				
		flag++;	   
	}
		

	for(i=0; i<order; i++){		/* Create mtd model file */
		fprintf(mtd_model,"Lag%d:\n",i+1);  															 	 	  	  	  	  
		for(j=0; j<dimension; j++){
			for(k=0; k<dimension; k++)
			{
				if(k==(dimension-1)){  fprintf(mtd_model,"%1.8lf", a[i][j][k]);  }
				else{  fprintf(mtd_model,"%1.8lf ", a[i][j][k]);  }
			}
			fprintf(mtd_model,"\n");	 																			 	 	 	 	 	 	 	 	 	 	 
		}															
	}		  	  
	fprintf(mtd_model,"\n");
		
	for(i=0; i<order; i++){	
		if( i==(order-1) ){  fprintf(mtd_model,"%lf", Lambda[i]);  }
		else{  fprintf(mtd_model,"%lf ", Lambda[i]);  } 	   	   	   	     
	}	

	fclose(mtd_model);						
}


