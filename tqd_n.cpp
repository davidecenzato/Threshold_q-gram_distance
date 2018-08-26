// program that calculate the threshold q-gram distance between two nucleotide strings

//Alessio Milanese
//13 November 2016

#include "Profile.cpp"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>


//----------TO SET:
int threshold = 1;
int max_string_length = 36000000; // you cannot read strings longer than this value



//----------------------------------------------------------------------------------//
//                                   GLOBAL VARIABLES                               //
//----------------------------------------------------------------------------------//




Profile profile(threshold);
int alphSize = 4; //alphabet size, here is tqd_n for nucleotides
int q; //length of the q-grams
char NucleicAlphabet [4] = {'a','t','c','g'};
char NucleicAlphabetUpperCase [4] = {'A','T','C','G'};


//----------------------------------------------------------------------------//
//                               FUNCTIONS DECLARATION                        //
//----------------------------------------------------------------------------//
// convert from string to an array of numbers
// (we doesn't work with string, but with arrays of integer)
void str2num (int num[], int NumDim, char str[]);

//load a fasta file
char* load_fasta(char file_name[]);

// calculate the profiles of two strings and save it in a RB-tree
void createProfile (int str1[], int strlen1, int str2[], int strlen2);

// calculate the position in the profile of a given q-gram, and calculate the
// q-gram given a position in the profile
unsigned long long int qposition (int vec[],int size_vec);
void qgramFromPosition (int vec[], unsigned long long int pos);

//useful function
unsigned long long int myPow(int a ,int b);





//----------------------------------------------------------------------------//
//                                   MAIN                                     //
//----------------------------------------------------------------------------//
int main(int argc, char **argv){

  //----------------------------------------------------------------------------
  // INPUT ---------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // q
  q = std::stoi (argv[1]);
  // threshold
  int threshold_val = std::stoi (argv[2]);
  // type_input: 1 strings, 2 files
  int type_input = std::stoi (argv[3]);
  // input 1
  char *input1 = argv[4];
  // input 2
  char *input2 = argv[5];
  // verbose level
  int verbose = std::stoi (argv[6]);

  //check the value for q
	if (q > 32){
    printf("Error. Maximum value for q is 32.\n");
    exit(1);
	}

	//----------------------------------------------------------------------------
	// INIZIALIZING THE VARIABLES
	//----------------------------------------------------------------------------
  if (verbose > 3){printf("Inizialise the variables \n");}
	char *str1;
	char *str2;

	//create a variable to calculate the time to load the fasta files
  time_t now1 = time(NULL);

  if (type_input == 1){  // the input are strings
    str1 = argv[4];
    str2 = argv[5];
  }

  if (type_input == 2){
    str1 = load_fasta(input1);
    str2 = load_fasta(input2);
  }

  time_t m1 = difftime(time(NULL), m1);


  //----------------------------------------------------------------------------
  //BEGINNING OF THE COMPUTATION
  //----------------------------------------------------------------------------
  //calculate the length of the strings ----------------------------------------
  if (verbose > 3){printf("Calculate the length of the strings.\n");}
  int len1 = strlen(str1);
  int len2 = strlen(str2);
  if (len1<q){
        printf("Error. the first string is shorter then q.\n");
        exit(1);
  }
  if (len2<q){
        printf("Error. the second string is shorter then q.\n");
        exit(1);
  }
  if (verbose > 3){printf("First string:  %d\n",len1);}
  if (verbose > 3){printf("Second string: %d\n",len2);}


	//convert the strings in arrays of numbers -----------------------------------
	if (verbose > 3){printf("Convert the strings in numbers.\n");}
	int *numericStr1;
	int *numericStr2;
	numericStr1 = (int *)malloc(len1*sizeof(int));
	numericStr2 = (int *)malloc(len2*sizeof(int));
	if (numericStr1 == NULL || numericStr2 == NULL){
		printf("\n Error.\n Out of memory.\n");
	}

	str2num(numericStr1,len1,str1);
	str2num(numericStr2,len2,str2);

  if (type_input == 2){
    free(str1); // free memory of the char string
    free(str2); // free memory of the char string
  }

	//----------calculate the profile

	if (verbose > 3){printf("Calculate the profile\n");}
	createProfile (numericStr1, len1, numericStr2, len2);

	//calculate the distance
	if (verbose > 3){printf("Calculate the distance\n");}
	unsigned long long int dist = profile.calculateDistance();

	//calculate the total pair-status
	if (verbose > 3){printf("Calcolate the pair statuses\n");}

	int *status;
	status=profile.get_total_status();

	//------------------------------------------------------------------------
	//PRINT THE RESULTS
	//------------------------------------------------------------------------
	printf("------------------------------------------\n");
	printf("RESULTS\n");
	printf("------------------------------------------\n");

	printf("=== INPUT: ===============================\n");

	printf("alphabet: a, t, c, g, A, T, C, G.\n");
	cout << "q: "<< q <<"\n";
	if(type_input==2){
		cout<< "first string: "<< input1 << "\n";
		cout<< "second string: "<< input2 << "\n\n";
	}


	printf("\n=== TIME: ================================\n");

	if(type_input==2){
		printf("time to load the fasta files (seconds): %ld\n",m1);
	}
//	m = difftime(time(NULL), now);
  //  printf("calculate time (seconds): %ld\n\n",m);


	printf("\n=== THRESHOLD Q-GRAM DISTANCE: ===========\n");

	cout << "number of nodes: " << profile.getNumber_of_nodes()<<"\n";
	//printf("percentage: %.2f\n\n",(float(dist)/n_nodes));

	cout<<"counting pair status in the profiles:\n";
	for (int i=0;i<3;i++){
		for (int j=0;j<3;j++){
			printf("%d - %d: %d \n",i,j,status[i*(threshold+2) + j]);
		}
	}
	printf("\n");

	cout << "\nTQD: "<<dist<<"\n";

}





//----------------------------------------------------------------------------------//
//                                      FUNCTIONS                                   //
//----------------------------------------------------------------------------------//

//////////////////////  LOAD FASTA  -------------------------------------------------------------
// load a fasta file

char* load_fasta(char file_name[]){
	// we create the string str that will contain the string that we are going to read
	char* str = (char *)malloc(max_string_length*sizeof(char));
	if (str == NULL){
		printf("\n Error.\n Out of memory.\n");
	}

	// create and read the FILE
	FILE *f;
	f = fopen(file_name, "r");

	if (f){
		//the first line is description
		while(fgetc(f) != '\n'){
		}
		//after the first '\n' I start to read the sequence
		int c;
		int cont=0;
		c = fgetc (f);
		while(c != EOF){
      		if (c != '\n' && c!='N'){ //--------------------------------------------we take out the N's
      			str[cont]=c;
      			cont++;
      		}
      		c = fgetc (f);
		}
		fclose(f);
	}else{
		printf("Error. \nThe file %s cannot be open.\n",file_name);
		exit(1);
	}
	return str;
}

//////////////////////  STRING => ARRAY OF NUMBERS  -------------------------------------------------------------
// convert a string in an array of number (according to the alphabet)
void str2num (int num[], int NumDim, char str[]){
	int correct=0;

	for(int i=0;i<NumDim;i++){
		for (int j=0;j<alphSize;j++){
			if (str[i]==NucleicAlphabet[j] || str[i]==NucleicAlphabetUpperCase[j] ){
				num[i]=j;
				break;
			}
			if(j==alphSize-1 && str[i]!=NucleicAlphabet[alphSize-1]  && str[i]!=NucleicAlphabetUpperCase[alphSize-1] ){ // if we have look at every character in the alphabet and we have not find the character
			   num[i]=-1;                                             // then there is an error
			   printf("\nError.\nThe character '%c' does not belong to the alphabet.",str[i]);
			   exit(1);
			}
		}
	}
}
//////////////////////  CREATE PROFILE  --------------------------------------------------------------------------
// calculate the profile of a given string

void createProfile (int str1[], int strlen1, int str2[], int strlen2){

	unsigned long long int pos=0;
	int qgram[q];

	// find the first q-gram of the FIRST string------------------------------
	for (int j=0;j<q;j++){
			qgram[j]=str1[j];
	}

	// calculate the position of the first q-gram and update the profile
	pos = qposition(qgram,q);
	profile.insert(pos,1);


	// calculate the entire profile
	for(int i=1;i<strlen1-q+1;i++){
		// calculate the next position using the previous one
		pos=(pos-( str1[i-1] *  myPow(alphSize,q-1) )) * alphSize + str1[i+q-1];
		//insert
		profile.insert(pos,1);
	}

	// calculate the profile for the second string
	pos=0;

	// find the first q-gram of the SECOND string------------------------------
	for (int j=0;j<q;j++){
			qgram[j]=str2[j];
	}

	// calculate the position of the first q-gram and update the profile
	pos = qposition(qgram,q);
	profile.insert(pos,2);

	// calculate the entire profile
	for(int i=1;i<strlen2-q+1;i++){
		// calculate the next position using the previous one
		pos=(pos-( str2[i-1] *  myPow(alphSize,q-1) )) * alphSize + str2[i+q-1];
		//insert
		profile.insert(pos,2);

	}
    return;
}
//////////////////////  Q-GRAM => PROFILE POSITION  -------------------------------------------------------------
// calculate the position of the q-gram passed as input, in the q-gram profile

unsigned long long int qposition (int vec[],int size_vec){
	unsigned long long int pos=0;
	unsigned long long int molt=1;
	for(int i=size_vec-1;i>=0;i--){
		pos=pos+vec[i]*(molt);
		molt=molt*alphSize;
	}
	return pos;
}
//////////////////////  PROFILE POSITION => Q-GRAM  -------------------------------------------------------------
// calculate the q-gram given the position of the q-gram in the profile


void qgramFromPosition (int vec[], unsigned long long int pos){
	for(int i=q-1;i>=0;i--){
		vec[i]=(int)pos%alphSize;
		pos=pos/alphSize;
	}
}

//////////////////////  MYPOW  -------------------------------------------------------------
// calculate the power a^b, and then saves the result in a long long int


unsigned long long int myPow(int a, int b){
	unsigned long long int res = a;
	int i = 1;
	for (i=1;i<b;i++){
		res = res*a;
	}
	return res;
}
