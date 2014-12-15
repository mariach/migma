/* MTD program. Copyright (C) 2012. Created by Maria Chatzou and Pantelis Bagos.
   This program is free software: you can redistribute it and/or modify 
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.For more details see
   the GNU General Public License <http://www.gnu.org/licenses/>. */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>	
#include <ctype.h>

#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "file_handle.h"

#define PROT_LENGTH 20
#define DNA_LENGTH 4
#define LINE_SIZE 5000


char **readFasta( char *in, int alphabet, int *seqNum)
{
	gzFile fp;
	kseq_t *seq;
	int l, seqCount=0;
	char **seqArray=malloc(1 * sizeof *seqArray);
	*seqNum=0;

	/*char *AB;
	if(alphabet==1){ AB="ARNDCEQGHILKMFPSTWYVBZ\n";  }
	if(alphabet==2){ AB="AaCcGgTt\n";  }*/


	if ((fp = gzopen(in, "r")) == NULL)
	{
		perror("Error: failed to open file for read\n");
		exit(1);
     	}	
     
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		
		/* Add name */
		seqArray = (char**) realloc ( seqArray, (seqCount+1)*sizeof(*seqArray) );
		if(seqArray==NULL){  exit(1);  }
		
		size_t idlen=0;
		if (seq->comment.l){   idlen= strlen(seq->name.s) + strlen(seq->comment.s) + 2;   }
		else{  idlen= strlen(seq->name.s) + 1; }
		
		char id[idlen];
		if (seq->comment.l){ sprintf(id,">%s %s",seq->name.s,seq->comment.s); }
		else{ sprintf(id,">%s",seq->name.s); }
		  
		seqArray[seqCount] = (char*)malloc( (idlen+1) * sizeof(char) );
		if( strcpy(seqArray[seqCount], "")==0 ){  exit(1);  }
		
		strcpy(seqArray[seqCount],id);  // printf("name: %d %s \n",seqCount, seqArray[seqCount]);
		
		seqCount++;
		
		
		/* Add seq */
		seqArray = (char**) realloc ( seqArray, (seqCount+1)*sizeof(*seqArray) );
		if(seqArray==NULL){  exit(1);  }
			 
		seqArray[seqCount] = (char*)malloc( (strlen(seq->seq.s)+1) * sizeof(char) );
		if(strcpy(seqArray[seqCount], "")==0){  exit(1);  }
	 
		strcpy(seqArray[seqCount],seq->seq.s); //printf("name: %d %s \n",seqCount, seqArray[seqCount]);
		
		seqCount++;
		
		(*seqNum)++;
	}
	seqArray = (char**) realloc ( seqArray, (seqCount+1)*sizeof(*seqArray) );
	if(seqArray==NULL){  exit(1);  }
			 
	seqArray[seqCount] = (char*)malloc(2 * sizeof(char) );
	if( strcpy(seqArray[seqCount], "")==0 ){  exit(1);  }
	
	seqArray[seqCount]='\0';
	
	kseq_destroy(seq);
	gzclose(fp);
	return seqArray;
}



void aa_search(char px, int *SID){
			
	char AB[]="AaRrNnDdCcEeQqGgHhIiLlKkMmFfPpSsTtWwYyVv";				
	char *pos = strchr(AB, px);	
		
	if( pos ){
	    if( (pos-AB)%2==0 ){ *SID= (pos-AB)/2; }
	    else{ *SID= ((pos-1-AB)/2);  }	
	}
	else if(px=='B'){  (*SID)=3;  }
	else if(px=='Z'){  (*SID)=7;  }
	else{
		printf("Error: unknown letter ' %c ' exists.\n ",px);
		exit(1);
	}
}



void dna_search(char px, int *SID){
			
	char AB[]="AaCcGgTt";				
	char *pos = strchr(AB, px);	
		
	if( pos ){   
	    if( (pos-AB)%2==0 ){ *SID= (pos-AB)/2; }
	    else{ *SID= ((pos-1-AB)/2);  }
	}
	else{
		printf("Error: unknown letter ' %c ' exists.\n ",px);
		exit(1);
	}
}



void abc(int *SIDi,int *SIDj, int x, char px[25]){
				
		int SID; 
		
		if(x==1){							
			aa_search(px[0],&SID);						
			(*SIDi)=SID;											
			
			aa_search(px[1],&SID);						
			(*SIDj)=SID;													
		}			
		else if(x==2){						
			dna_search(px[0],&SID);						
			(*SIDi)=SID;											
			
			dna_search(px[1],&SID);						
			(*SIDj)=SID;					
		}			
		else{ printf("\nError: wrong choise of alphabet.\n");}
}



void find_Dimensios(char *in, int *order, int *I, int *J){
	
	char line[LINE_SIZE], *token=NULL;
	const char dlm[] =" ";
	int lag_count=0, I_count=0, J_count=0, flag=0, count=0;

	FILE *file=fopen(in,"r");

	while (fgets(line, LINE_SIZE, file)!=NULL)
	{  
		if(strncmp("Lag",line,3) == 0)
		{	     	 	  	  	     	  	  
			lag_count++;			
			flag=1;	 
	
			if(J_count != I_count){ 
				perror("\nError: wrong MTD model used. Choose another MTD file or built a new MTD model.");
				exit(1);
			}
			J_count=0; 
			I_count=0;			 	   	   
		}
		else{ 
			if(line[0]!='\n')
			{	
				J_count++; 

				token=strtok(line,dlm);
				while(token!=NULL){						    	
					count++; 															
					token=strtok(NULL,dlm);	 	 	 	 	 	 											   	 	 	    			   	    								   	    					 	 	 	 	 
				}
			      	
				if(flag==1){ I_count=count;   flag=0;  }
				if( (J_count<=I_count && count!=I_count) || (J_count>I_count && count!=lag_count) ){ 
					printf("%d %d %d\n\n", I_count, J_count, lag_count);
					perror("\nError: wrong MTD model used. Choose another MTD file or built a new MTD model.");
					exit(1);
				}
			 	count=0;
			}	
		}		 	 	 	 	 	 
	}
	fclose(file);
	
	(*order)=lag_count;
	(*J)=J_count-1;		// -1 because it also counts the last line (Lambda line)
	(*I)=I_count;
}




void read_modelFile(char *in, int order, int I, int J, double lambda[order], double a[order][I][J]){
																											
	FILE *file=fopen(in,"r");
	if( file == NULL )
	{
	    printf("Error while opening the file: %s\n", in);
	    exit(EXIT_FAILURE);
	}
	
	char line[LINE_SIZE];	 
	int count1=0, count2=0, count3=0, lcount=0;	   	   	   	   	   	   
	char *token;
	const char dlm[] =" ";	  
	double nNumber;
	
	int k,i,j;                                          
	for(k=0; k<order; k++){                                               
		for(i=0; i<I; i++){                               
			for(j=0; j<J; j++){	 
				 a[k][i][j]=0; 	   	   	             	 	 	 	 	
			}                     	            
		}                                                                                
	}       
	
				
	while (fgets(line, LINE_SIZE, file)!=NULL)	
	{		 
		if(strncmp("Lag",line,3) == 0){	     			 												    
			count1++; 
			count2=0; 
		}
		else{	 
			    token=strtok(line,dlm);																			 
			
			    while(token!=NULL){	   
					
					nNumber=atof(token);                                                        			
					
					if(order>(count1-1) && J>count3 && I>count2){  a[count1-1][count2][count3]=nNumber;	}	
					if(I<count2){
							if((strlen(line)==1) && (strchr( line, '\n' )!=NULL)){ break;}
							else{
								lambda[lcount]=nNumber;	   	   	   	   	   	   	   	   	   	   			 	
								lcount++;  
							}	 	 
					}
					token=strtok(NULL,dlm);						
					count3++; 	  	  	  	  	    			   	   	   	      	   	   	   	   	   	   	   	   	      	   	   	   	   	   	    				
			    }																																							  	  	  
			    count3=0;
			    count2++;
		}                                                                                   
	}    
	fclose(file); 
}

