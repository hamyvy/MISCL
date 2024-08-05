#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <cstring>
using namespace std;
//------------------------------------------------------------------------------

int main(int argc, char** argv)
{
	FILE *infile,*outfile,*infile2;
	int i,na,nd,nSam,alpha;
	long site;
	float likelihood,mn,tl,*mean,*tail;
	char str[50];
	
	if(argc<4){
		printf("pls give input/output files' names\n");
		exit(-1);
	}

	
	if ( (infile = fopen( argv[1], "r" )) == NULL )
		{
			printf("Could not find input file 1\n");
			exit(-1);
		}
	if ( (infile2 = fopen( argv[2], "r" )) == NULL )
		{
			printf("Could not find input file 2\n");
			exit(-1);
		}
	if ( (outfile = fopen( argv[3], "w")) == NULL )
        {
                printf("Could not open output file \n");
                exit(-1);
        }
 	
 	fscanf(infile2, "%d%*d%*d", &nSam );
 	if(nSam <= 0){
	 	printf("\n file 2 error\n");
	 	exit(-1);
	 }
	 
	mean = new float[nSam];
	tail = new float[nSam];
	for(i=0;i<nSam;i++){
		mean[i] = 0;
		tail[i] = 0;
	}
	
	while ( !feof(infile2) ){
		fscanf(infile2, "%d%f%f", &nd,&mn,&tl );
		if(nd < nSam){
			mean[nd] = mn;
			tail[nd] = tl;
		}
		else{
			printf("\nTable error\n");
		}
	}
	
	fclose(infile2);
	
	while ( !feof(infile) ){
		fscanf(infile, "%s%ld%d%d%d%f", str,&site,&nd,&na,&alpha,&likelihood );
		if(nd<nSam && tail[nd] != 0){
			likelihood = (likelihood-mean[nd])/(tail[nd]-mean[nd]);
			fprintf(outfile,"%s %ld %d %d %d %f\n",str,site,nd,na,alpha,likelihood);
		}
		else{
			fprintf(outfile,"%s %ld %d %d %d Null\n",str,site,nd,na,alpha);
		}
	}
	
	fclose(infile);
	fclose(outfile);
	return 0;
}

