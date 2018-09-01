// program that finds the maximum length of the string and the dictionary used
// Alessio Milanese <milanese@embl.de>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <iostream>

//----------------------------------------------------------------------------//
//                                   MAIN                                     //
//----------------------------------------------------------------------------//
int main(int argc, char **argv){
  // input 1
  char *file_name = argv[1];
  // search the alphabet = 1; do not search the alphabet = 2
  // note that to find the alphabet the time complexity is O(len_str * len_alphabet)
  int search_alph = std::stoi (argv[2]);

  // record the maximum length
  int max_length = 0;

  // Alphabet
  char alphabet[100];
  for(int i = 0;i<100;i++){
    alphabet[i] = ' ';
  }
  int cont_alp = 0;

	// create and read the FILE
	FILE *f;
	f = fopen(file_name, "r");

	if (f){
    // check first character, it should be '>'
    if(fgetc(f) != '>'){
      fprintf( stderr,"Error. \nThe file %s is not a fasta file.\n",file_name);
  		exit(1);
    }
		//the first line is description
		while(fgetc(f) != '\n'){
		}
    // END OF FIRST LINE -------------------------------------------------------
		//after the first '\n' I start to read the sequence
		int c;
		int current_length=0;
		c = fgetc (f);
		while(c != EOF){
      		if (c != '\n'){
      			current_length++;
      		}
          // if it is a fasta header
          if (c == '>'){
            // update max length
            if (current_length > max_length){
              max_length = current_length;
            }
            // read all the fasta header
      			current_length = 0;
            c = fgetc (f);
            while(fgetc(f) != '\n'){
          		}
          // end of fasta header
        }else{ // if it is not '>' we add the char to the alphabet
           if (search_alph == 1){
             if (c != '\n'){
             int found = 0;
             for(int i = 0;i<100;i++){
               if(alphabet[i] == c){
                 found = 1;
                 break;
               }
             }
             if(found == 0){ // if the char is new
               alphabet[cont_alp] = c;
               cont_alp++;
             }
           }

           }
        }
      		c = fgetc (f);
		}
		fclose(f);
	}else{
		fprintf( stderr,"Error. \nThe file %s cannot be open.\n",file_name);
		exit(1);
	}

  fprintf( stdout,"%d",max_length);
  if (search_alph == 1){
    fprintf( stdout,"\t%c",alphabet[0]);
    for(int i = 1;i<100;i++){
      if(alphabet[i] == ' '){
        break;
      }
      fprintf( stdout,"\t%c",alphabet[i]);
    }
  }
  fprintf( stdout,"\n");
  exit(0);

}
