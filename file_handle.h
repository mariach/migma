/* MTD program. Copyright (C) 2012. Created by Maria Chatzou and Pantelis Bagos.
   This program is free software: you can redistribute it and/or modify 
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.For more details see
   the GNU General Public License <http://www.gnu.org/licenses/>. */


#ifndef _file_handle_h
#define _file_handle_h

char **readFasta( char *in, int alphabet, int *seqNum);


void aa_search(char px, int *SID);


void dna_search(char px, int *SID);


void abc(int *SIDi,int *SIDj, int x, char px[25]);


void find_Dimensios(char *in, int *order, int *I, int *J);


void read_modelFile(char *in, int order, int I, int J, double lambda[order], double a[order][I][J]);


#endif /*_file_handle_h */
