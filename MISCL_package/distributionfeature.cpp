#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <cstring>
using namespace std;
//------------------------------------------------------------------------------

int main(int argc, char** argv)
{
	FILE *infile,*outfile;
	int i,j,k,na,nd,nSam,mh;
	long *nd_long,*ntail;
	float *mean,*tail,likelihood,maxLf,percentile;
	vector<float> :: iterator p;
	
	if(argc<3){
		printf("pls give input/output files' names\n");
		exit(-1);
	}
	percentile = 0.002;
	
	for(i=1;i<argc;i++){
		if(!strcmp(argv[i], "-n")){
			if(i!=argc-1){
				if(argv[i+1][0]!='-') nSam = atoi(argv[++i]);
			}
			else{
				printf("sample size is missing\n");
				exit(-1);
			}
		}
		if(!strcmp(argv[i], "-p")){
			if(i!=argc-1){
				if(argv[i+1][0]!='-') percentile = atof(argv[++i]);
			}
		}
	}
	
	if ( (infile = fopen( argv[1], "r" )) == NULL )
		{
			printf("Could not find input file\n");
			exit(-1);
		}	
	if ( (outfile = fopen( argv[2], "w")) == NULL )
        {
                printf("Could not open output file \n");
                exit(-1);
        }

 	vector<float> Lf0[nSam];
 	
	nd_long = new long[nSam];
	ntail = new long[nSam];
	mean = new float[nSam];
	tail = new float[nSam];
	for(i=0;i<nSam;i++){
		mean[i] = 0.0;
		tail[i] = 0.0;
		nd_long[i] = 0;
		ntail[i]= 0;
	}
	
	fscanf(infile, "%*s%*s%*s%*s%*s%*s" );
	
	while ( !feof(infile) ){
		fscanf(infile, "%*s%*s%d%d%*s%f", &nd,&na,&likelihood );
		if(na+nd==nSam){
			mean[nd] += likelihood;
			nd_long[nd]++;
		}
	}
	
	for(i=0;i<nSam;i++)
		if(nd_long[i]>0){
			mean[i] = mean[i]/nd_long[i];
			ntail[i]= long(nd_long[i]*percentile);
		}
	
	rewind(infile);
	
	while ( !feof(infile) ){
		fscanf(infile, "%*s%*s%d%d%*s%f", &nd,&na,&likelihood );
		if(na+nd==nSam){
			if(likelihood>mean[nd])
				Lf0[nd].push_back(likelihood);
		}
	}
	
	fprintf(outfile,"%d 0 0\n",nSam);
	for(i=0;i<nSam;i++){
		if(nd_long[i]>0){
			for(j=0;j<ntail[i];j++){
				maxLf=0;mh=0;
				for(k=0;k<Lf0[i].size();k++)
					if(maxLf<Lf0[i][k]){
							maxLf = Lf0[i][k];
							mh = k;
					}
				tail[i] += maxLf;
				p=Lf0[i].begin();
				p += mh;
				Lf0[i].erase(p);
			}
			fprintf(outfile,"%d %f %f\n",i, mean[i], tail[i]/ntail[i]);
		}
	}
	fclose(infile);
	fclose(outfile);
	return 0;
}

