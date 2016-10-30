/* Match every k-character snippet of the query_doc document
	 among a collection of documents doc1, doc2, ....

	 ./rkmatch snippet_size query_doc doc1 [doc2...]

*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <strings.h>
#include <assert.h>
#include <time.h>
#include <ctype.h>
#include <string.h>


#include "bloom.h"

enum algotype { SIMPLE = 0, RK, RKBATCH};

/* a large prime for RK hash (BIG_PRIME*256 does not overflow)*/
long long BIG_PRIME = 5003943032159437; 

/* constants used for printing debug information */
const int PRINT_RK_HASH = 5;
const int PRINT_BLOOM_BITS = 160;

/* modulo addition */
long long
madd(long long a, long long b)
{
	return ((a+b)>BIG_PRIME?(a+b-BIG_PRIME):(a+b));
}

/* modulo substraction */
long long
mdel(long long a, long long b)
{
	return ((a>b)?(a-b):(a+BIG_PRIME-b));
}

/* modulo multiplication*/
long long
mmul(long long a, long long b)
{
	return ((a*b) % BIG_PRIME);
}

/* read the entire content of the file 'fname' into a 
	 character array allocated by this procedure.
	 Upon return, *doc contains the address of the character array
	 *doc_len contains the length of the array
	 */
void
read_file(const char *fname, char **doc, int *doc_len) 
{
	struct stat st;
	int fd;
	int n = 0;

	fd = open(fname, O_RDONLY);
	if (fd < 0) {
		perror("read_file: open ");
		exit(1);
	}

	if (fstat(fd, &st) != 0) {
		perror("read_file: fstat ");
		exit(1);
	}

	*doc = (char *)malloc(st.st_size);
	if (!(*doc)) {
		fprintf(stderr, " failed to allocate %d bytes. No memory\n", (int)st.st_size);
		exit(1);
	}

	n = read(fd, *doc, st.st_size);
	if (n < 0) {
		perror("read_file: read ");
		exit(1);
	}else if (n != st.st_size) {
		fprintf(stderr,"read_file: short read!\n");
		exit(1);
	}
	
	close(fd);
	*doc_len = n;
}


/* The normalize procedure examines a character array of size len 
	 in ONE PASS and does the following:
	 1) turn all upper case letters into lower case ones
	 2) turn any white-space character into a space character and, 
	    shrink any n>1 consecutive spaces into exactly 1 space only
			Hint: use C library function isspace() 
	 You must do the normalization IN PLACE so that when the procedure
	 returns, the character array buf contains the normalized string and 
	 the return value is the length of the normalized string.
*/

int normalize(char *buf,	/* The character array containing the string to be normalized*/
					int len			/* the size of the original character array */)
{
  int flag = 0;
	char *current=buf;
	char *checker=buf;
	int firstCharIsLetter; // Checks wheather the first character is a letter
	if (isspace(buf[0]))
		firstCharIsLetter = 0;
	else
		firstCharIsLetter = 1;
	int index = 0; //Keeps track of the number of iterations
	int length = len; //length of the char array, needed so that I can change the len variable in the while loop	
	while(index<length){
		
		/* Convert to lowercase*/
		if(isupper(*checker)){
			*checker=tolower(*checker);
		}

		/* Check whether it starts with a whitespace character*/
		if (index==0){
			if (!firstCharIsLetter){
			  while((isspace(*checker))&&(index<length)){
					checker++;
					len--;	
					index++;
			  }
			  *current = tolower(*checker);
				firstCharIsLetter=1;
				current++;
				checker++;	
				index++;
			}else{
			  index++;
			  checker++;
			  current++;
			}
		}else{
		  
			/* Check whether a character is a whitespace character
			   whenever checker moves the index gets incremented
			 */	
			if(isspace(*checker)){
				index++;
				if(isspace(*(checker++))){	
					*(current++)=' ';
					while((isspace(*checker))&&(index<length)){
					  
					  checker++;
					  len--;	
					  index++;
					        
						
					}
					*current =tolower(*checker);
					current++;
					checker++;
					index++;
					if(flag==1)
					  break;
				}
			
			}
			else{
			  *current = tolower(*checker);
				current++;
				checker++;
				index++;
			}
				
		}
	}
	//if the last character is a space increment the len
	if(buf[len-1]==' ')
		len--;

	return len;

}

/* check if a query string ps (of length k) appears 
	 in ts (of length n) as a substring 
	 If so, return 1. Else return 0
	 You may want to use the library function strncmp
	 */
int
simple_match(const char *ps,	/* the query string */
						 int k, 					/* the length of the query string */
						 const char *ts,	/* the document string (Y) */ 
						 int n						/* the length of the document Y */)
{
  for(int i=0; i<n; i++){
    if (strncmp(ts, ps, k)==0)
      return 1;
    ts++;
  }
	
  return 0;
}

/* Check if a query string ps (of length k) appears 
	 in ts (of length n) as a substring using the rabin-karp algorithm
	 If so, return 1. Else return 0
	 In addition, print the first 'PRINT_RK_HASH' hash values of ts
	 Example:
	 $ ./rkmatch -t 1 -k 20 X Y
	 605818861882592 812687061542252 1113263531943837 1168659952685767 4992125708617222 
	 0.01 matched: 1 out of 148
	 */

int
rabin_karp_match(const char *ps,	/* the query string */
								 int k, 					/* the length of the query string */
								 const char *ts,	/* the document string (Y) */ 
								 int n						/* the length of the document Y */ )
{
  long long int val=0;            //The variable that stores computed hash values
  long long int power=1;          //The variable that is used to compute power
  long long int hashForQuerry;    //Hash value for query string
  long long int hashForString;    //Hash value for document
  long long int headOfHashed;     // keeps the value of 256^(k)
  int innerK=k;                   //used to compute the power
  int charInInt;
 
  
  for(int i = 0; i < k; i++){
    power=1;
    for(int j = 1; j < (innerK); innerK--){ //Computes 256^(innerK), where innerK depends on i
      power=mmul(256,power);
	  
    }
    charInInt = ps[i];
    val=madd(val,mmul(power,charInInt));  //updates the value as i gets incremented, val+= ps[i]*power
    innerK=k-i-1;
  }

  hashForQuerry = val;
  innerK=k;
  val = 0;
  int new_line = 0;

  //Computes hash value for each of the n-k pieces of the document Y
  for(int i=0; i<n-k; i++){
    if (i==0){
      for(int j=0; j<k; j++){      
	power=1;
        for(int f = 1; f < (innerK); innerK--){
	  power=mmul(256,power);   //Computes 256^(innerK), where innerK depends gets decremented and depends on k
			
        }
	if(j==0){
	  headOfHashed=power; //Keeps track of 256^(k-1)
	}              
	charInInt = ts[j];
	val=madd(val,mmul(power,charInInt));
	innerK=k-j-1;
	
      }
      hashForString = val;
    }
    else{          // Computes the rolling hashing, for all except the first parts of the document Y
      val = 0; 
      charInInt = ts[i];
      val = mmul(headOfHashed,ts[i-1]);  //256^(k-1)*ts[i]
      val = mdel(hashForString, val);   // y-(256^(k-1)*ts[i])
      val = mmul(256,val);             //(y-(256^(k-1)*ts[i]))*256
      charInInt = ts[i+k-1];
      val = madd(val,charInInt);   
      hashForString = val;
      
    }
    //Prints the hash values
    if(i<5)
      printf("%lld ", hashForString);
    if(i==5){
      printf("\n");
      new_line = 1;
    }

    /*Checking if there is a match, 
      making sure it is an actual match and 
      returning the value based on the result*/
    
    if(hashForQuerry==hashForString){
		if(strncmp((ts+i),ps, k)==0){
		  if (new_line == 0)
		    printf("\n");
		  return 1;
		}
    }
    else{
    }
    
  }
  
  return 0;
}

/* Initialize the bitmap for the bloom filter using bloom_init().
	 Insert all m/k RK hashes of qs into the bloom filter using bloom_add().
	 Then, compute each of the n-k+1 RK hashes of ts and check if it's in the filter using bloom_query().
	 Use the given procedure, hash_i(i, p), to compute the i-th bloom filter hash value for the RK value p.

	 Return the number of matched chunks. 
	 Additionally, print out the first PRINT_BLOOM_BITS of the bloom filter using the given bloom_print 
	 after inserting m/k substrings from qs.
*/



int
rabin_karp_batchmatch(int bsz,        /* size of bitmap (in bits) to be used */
                      int k,          /* chunk length to be matched */
                      const char *qs, /* query docoument (X)*/
                      int m,          /* query document length */ 
                      const char *ts, /* to-be-matched document (Y) */
                      int n           /* to-be-matched document length*/)

{
  long long int val=0; //The variable that stores computed hash values
  long long int power=1; //The variable that is used to compute power
  long long int hashForQuerry; //Hash value for query string
  long long int hashForString; //Hash value for document
  long long int headOfHashed; // keeps the value of 256^(k)
  int innerK=k; //used to compute the power
  int charInInt;

  //Initializing bloom filter
  bloom_filter Bloom = bloom_init(bsz);

  int query_res;
  int num_of_matches=0;

  //Computes the hash value and adds them to a bloom filter for each of m/k substrings of size k
  for(int f=0; f<m; f+=k){
    innerK=k;
    val = 0;
    for(int i = 0; i < k; i++){
      power=1;
      for(int j = 1; j < (innerK); innerK--){ //Computes 256^(innerK), where innerK depends on i
	power=mmul(256,power);
	  
      }
      charInInt = qs[f+i];
      val=madd(val,mmul(power,charInInt));  //updates the value as i gets incremented, val+= ps[i]*power
      innerK=k-i-1;
    }
    hashForQuerry = val;
    
    //adding to the Bloom filter
    bloom_add(Bloom, hashForQuerry);
  }

  //prints the first PRINT_BLOOM_BITS number of bits
  bloom_print(Bloom, PRINT_BLOOM_BITS);

  //reinitializing the variables
  innerK=k;
  val = 0;

  for(int i=0; i<n-k; i++){
    if (i==0){
      for(int j=0; j<k; j++){ //Computes hash value for each of the n-k pieces of the document
	power=1;
        for(int f = 1; f < (innerK); innerK--){
	  power=mmul(256,power);   //Computes 256^(innerK), where innerK is getting decremented and reinitialized at the end of the loo 
			
        }
	if(j==0){
	  headOfHashed=power; //Keps track of 256^(innerK)
	}              
	charInInt = ts[j];
	val=madd(val,mmul(power,charInInt)); //updates the value as i gets incremented, val+= ps[i]*power
	innerK=k-j-1; // innerK is reinitialized
	
      }
      
    }
    else{          // Computes the rolling hashing, for all except the first parts of the document Y
      val = 0;
      charInInt = ts[i];     
      val = mmul(headOfHashed,ts[i-1]);      
      val = mdel(hashForString, val);      
      val = mmul(256,val);
      charInInt = ts[i+k-1];    
      val = madd(val,charInInt);      
    }
    
    hashForString = val;
    
    //Searches in the bloom filter for each of the strings 
    query_res = bloom_query(Bloom, hashForString);
    
    //if a match checks wheather the substring is in X
    if(query_res){
      for(int index=0; index<m; index+=k){
	if(strncmp((qs+index), ts+i, k)==0)
	  num_of_matches +=1;
      }
    }  
  }
  
  return num_of_matches;
}
  

int 
main(int argc, char **argv)
{
 
	int k = 100; /* default match size is 100*/
	int which_algo = SIMPLE; /* default match algorithm is simple */

	char *qdoc, *doc; 
	int qdoc_len, doc_len;
	int i;
	int num_matched = 0;
	int to_be_matched;
	int c;

	/* Refuse to run on platform with a different size for long long*/
	assert(sizeof(long long) == 8);

	/*getopt is a C library function to parse command line options */
	while (( c = getopt(argc, argv, "t:k:q:")) != -1) {
		switch (c) 
		{
			case 't':
				/*optarg is a global variable set by getopt() 
					it now points to the text following the '-t' */
				which_algo = atoi(optarg);
				break;
			case 'k':
				k = atoi(optarg);
				break;
			case 'q':
				BIG_PRIME = atoi(optarg);
				break;
			default:
				fprintf(stderr,
						"Valid options are: -t <algo type> -k <match size> -q <prime modulus>\n");
				exit(1);
			}
	}

	/* optind is a global variable set by getopt() 
		 it now contains the index of the first argv-element 
		 that is not an option*/
	if (argc - optind < 1) {
		printf("Usage: ./rkmatch query_doc doc\n");
		exit(1);
	}

	/*argv[optind] contains the query_doc argument */

	/*for(int i=0; i<qdoc_len; i++){
	  printf("%c", qdoc[i]);
	}
	printf("\n\n\n");*/

	
	read_file(argv[optind], &qdoc, &qdoc_len);
	qdoc_len = normalize(qdoc, qdoc_len);

	/*for(int i=0; i<qdoc_len; i++){
	  printf("%c", qdoc[i]);
	}
	printf("\n\n\n");*/
        
	/* argv[optind+1] contains the doc argument */
	read_file(argv[optind+1], &doc, &doc_len);
	doc_len = normalize(doc, doc_len);

	
	switch (which_algo) 
		{
			case SIMPLE:
				/* for each of the qdoc_len/k chunks of qdoc, 
					 check if it appears in doc as a substring*/
			  
				for (i = 0; (i+k) <= qdoc_len; i += k) {
					if (simple_match(qdoc+i, k, doc, doc_len)) {
						num_matched++;
					}
				}
				break;
			case RK:
				/* for each of the qdoc_len/k chunks of qdoc, 
					 check if it appears in doc as a substring using 
				   the rabin-karp substring matching algorithm */
				for (i = 0; (i+k) <= qdoc_len; i += k) {
					if (rabin_karp_match(qdoc+i, k, doc, doc_len)) {
						num_matched++;
					}
				}
				break;
			case RKBATCH:
				/* match all qdoc_len/k chunks simultaneously (in batch) by using a bloom filter*/
				num_matched = rabin_karp_batchmatch(((qdoc_len*10/k)>>3)<<3, k, qdoc, qdoc_len, doc, doc_len);
				break;
			default :
				fprintf(stderr,"Wrong algorithm type, choose from 0 1 2\n");
				exit(1);
		}
	
	to_be_matched = qdoc_len / k;
	printf("%.2f matched: %d out of %d\n", (double)num_matched/to_be_matched, 
			num_matched, to_be_matched);

	free(qdoc);
	free(doc);

	return 0;
}
