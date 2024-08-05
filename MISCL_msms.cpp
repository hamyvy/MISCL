#include "MISCL.h"

float *pHomo;
double ****Table;
vector<int *> allele;
FILE *inputdata,*outfile;
//------------------------------------------------------------------------------
void CbTable()
//calculate combination(n,k)
{
	int i,k,m,n;
	float C;
	
	Cb = new float*[2*nSam+1];
	for(n=0;n <= 2*nSam;n++) Cb[n] = new float[2*nSam+1];
	
	for(n=0;n <= 2*nSam;n++)
		for(k=0;k <= n/2;k++){
			m=n-k;
			if(m==0 || k==0) C=1.0;
			else if(m<0 ||n==0) C=0.0;
			else{
				C=1.0;
				for(i=m+1;i<=n;i++) C=C*i/(i-m);
			}
			Cb[n][k]= C;
		}
	for(n=0;n <= 2*nSam;n++)
		for(k=n/2+1;k <= 2*nSam;k++){
			if(k>n) Cb[n][k]= 0.0;
			else Cb[n][k] = Cb[n][n-k];
		} 
}
//------------------------------------------------------------------------------
void subPop_freqSpec()
//calculate neutral sampling probability of having k derived alleles in sample of size n 
//(n<nSam due to missing data)
{
	int i,j,n;
	
	subP = new float*[nSam+1]; //subP is the neutral sub_sampling probability spectrum
	for(n=1;n <= nSam;n++) subP[n] = new float[n];
	
	for(n=1;n <= nSam;n++){
		if(n==nSam) for(j=1;j<n;j++) subP[n][j]=1.0/j;
		else
			for(j=1;j<n;j++){
				float pnj = 0.0;
				for(i=j;i<nSam;i++)
					if((nSam-i-n+j)>=0)
						pnj += Cb[i][j]*Cb[nSam-i][n-j]/Cb[nSam][n]/i;
				subP[n][j]=pnj;
			}
	}
}
//------------------------------------------------------------------------------
double ProbNielsen(int nSam1,int nSam2,int Ze,int k1,int k2)
{
	int i,ntt,ktt,befSwp_n;
	float Pe;
	double p1,p2,p_i, P;
	
	Pe = Ze/100.0;
	ntt = nSam1+nSam2;
	ktt = k1+k2;
	P=pow(Pe,nSam1)*Cb[nSam1][k1]/Cb[ntt][ktt]*subP[ntt][ktt];
	p_i=0;
	
	for(i=0;i<nSam1;i++){
		befSwp_n = nSam2+i+1;
		p1=0;
		if((k1+1+i-nSam1)>0)
			p1 = subP[befSwp_n][k1+1+i-nSam1+k2]*Cb[i][nSam1-k1]/Cb[befSwp_n][ntt-ktt];	

		p2=0;
		if((1+i-k1)>0)
			p2 = subP[befSwp_n][ktt]*Cb[i][k1]/Cb[befSwp_n][ktt];
			
		p_i += (p1+p2)*Cb[nSam1][i]*pow(Pe,i)*pow(1-Pe,nSam1-i);
	}
	P += p_i;
	P = P*Cb[nSam2][k2];
	return P;
}
//------------------------------------------------------------------------------
void make_interval()
//for a given selection coefficient s, sequences are devided into intervals associating with Ze intervals
{
	float delta_alpha,alpha,R;
	
	delta_alpha = max_alpha/100;
	for(int s=1;s<11;s++){
		alpha = s*s*delta_alpha;
		AlphaTable[s]=alpha;
		interval[s][0]=1;
		for(int z=1;z<100;z++){
			R = -2*alpha*log(1-z/100.0)/log(2*alpha);
			interval[s][z] = long(R/rho);
		}
		for(int z=1;z<99;z++) interval[s][z] = (interval[s][z]+interval[s][z+1])/2;
	}
}
//------------------------------------------------------------------------------
void readmsout()
{
	int i, j, numSg;
	char Line[500];

	Sites.clear();
	for(j=0;j<allele.size();j++)
    	delete[](allele[j]);
  	allele.clear();
	numSg = 0;

	while( fscanf(inputdata,"%s", Line), !strstr(Line, "//") );
	fscanf(inputdata, "%s%d", Line, &numSg );
	if ( !strstr(Line, "segsites") ){
			printf("cannot find segsites\n");
			exit(-1);
		}

	fscanf(inputdata, "%s", Line );
	if ( !strstr(Line, "positions") ){
 			printf("cannot find positions\n");
    		exit(-1);
        }

	for (i=0; i<numSg; i++)
	{
		float position=0.0;
	 	fscanf(inputdata,"%f", &position );
		Sites.push_back((long)(seq_length*position+0.5));
		if ( Sites[i] == 0 ) Sites[i] = 1;		
		if ( i>0 && Sites[i] <= Sites[i-1] )
			Sites[i] += (Sites[i-1]-Sites[i]+1);
	}
	for (i=0; i<nSam; i++)
	{
		char sstr[numSg+2];
		fscanf(inputdata,"%s",sstr);
		int *pl = new int[numSg];
		for(j=0;j<numSg;j++) pl[j]=0;
		for(j=0;j<numSg;j++){
	  		if(sstr[j]=='1') pl[j]= 1;
	  		else if(sstr[j]=='0') pl[j]= 0;
		}
		allele.push_back(pl);
	}
	
}
//------------------------------------------------------------------------------
void InputProbability()
{
	int i,j,n,n1,n2,z;
	
//for calculating Pn[n1][n2][0][0], it is calculated by (1-probability of being SNP)
	
	pHomo = new float[n1size];
	
	for(n=0; n< n1size; n++)
	{
		n1=n+minn1;
		n2=nSam-n1;
		pHomo[n] = 0.0;
		for(i=0; i<=n1; i++)
		for(j=0; j<=n2; j++)
			if( (i+j)!=0 && (i+j)!=nSam )
		 		pHomo[n] += Cb[n1][i]*Cb[n2][j]/Cb[nSam][i+j]/(i+j);
	}
//P[n1][n2][k1][k2] for selective sweep model
	//assign space
	Table = new double***[n1size];
	for(n=0; n<n1size; n++){
		Table[n] = new double**[100];
		n1=n+minn1;
		n2=nSam-n1;
		for(z=0; z< 100; z++){
			Table[n][z] = new double*[n1+1];
			for(i=0; i<= n1; i++)
				Table[n][z][i] = new double[n2+1];
		}
	}
	//calculating
	for(n=0; n<n1size; n++){
		n1=n+minn1;
		n2=nSam-n1;
		for(z=1; z< 100; z++){
			Table[n][z][0][0] = 0.0;
			for(i=0; i<= n1; i++)
			for(j=0; j<= n2; j++)
			if( (i+j)!=0 && (i+j)!= nSam ){
				Table[n][z][i][j] = ProbNielsen(n1,n2,z,i,j);
				Table[n][z][0][0] += Table[n][z][i][j];
			}
		}
	}
	
	printf("table done ");
}
//-------------------------------------------------------------------------------
float CalcLikelihood(int n1,int s,int core,int **groupInfo)
{
	float LnL1,LnN;
	int i,k1,k2,leftzmax,rightzmax,n1t,z,ktt,n2;
	long x,numSg,leftbound,rightbound;
	
	x=Sites[core];
	LnL1 = 0.0; LnN=0.0;
	n1t = n1-minn1; n2=nSam-n1;
	leftbound = x; rightbound = seq_length -x;
	numSg = Sites.size();
	
	if(x>interval[s][99]){
		leftzmax = 99;
		leftbound = interval[s][99];
	}
	else{
		i=1;
		while(1){
			if(x <= interval[s][i]){
				leftzmax = i;  break;
			} 
			i++;
			if(i>99) {printf("'while1' error"); return 0;}
		}
	}

	if(rightbound>interval[s][99]){
		rightzmax = 99;
		rightbound = interval[s][99];
	}
	else{
		i = 1;
		while(1){
			if(rightbound<=interval[s][i]){
				rightzmax = i;  break;
			}
			i++;
		}
	}

	int index = core-1; int nSNP; int nSNPtt = 0;
	
//leftside
	for(z=1;z<=leftzmax;z++){
		nSNP = 0;
		while(index >= 0 && (x-Sites[index])<interval[s][z]){
			k1 = groupInfo[index][0];
			k2 = groupInfo[index][1];
			ktt = k1+k2;
			if( ktt==0 || ktt==nSam ){
				LnL1 += log(1.0-theta*Table[n1t][z][0][0]);
				LnN += log(1.0-theta*pHomo[n1t]);
			}
			else{
				LnL1 += log( theta*Table[n1t][z][k1][k2] );
				LnN += log(theta*Cb[n1][k1]*Cb[n2][k2]/Cb[nSam][ktt]/ktt);
			}
			index--;
			nSNP++;
		}
		if(z==leftzmax)
			LnL1 += (leftbound-interval[s][z-1]-nSNP)*log(1.0-theta*Table[n1t][z][0][0]);
		else
			LnL1 += (interval[s][z]-interval[s][z-1]-nSNP)*log(1.0-theta*Table[n1t][z][0][0]);
		nSNPtt += nSNP;
	}
//right side
	index=core+1;
	for(z=1;z<rightzmax;z++){
		nSNP=0;
		while(index < numSg && (Sites[index]-x)<interval[s][z]){
			k1 = groupInfo[index][0];
			k2 = groupInfo[index][1];
			ktt = k1+k2;
			if( ktt==0 || ktt==nSam ){
				LnL1 += log(1.0-theta*Table[n1t][z][0][0]);
				LnN += log(1.0-theta*pHomo[n1t]);
			}
			else{
				LnL1 += log( theta*Table[n1t][z][k1][k2] );
				LnN += log(theta*Cb[n1][k1]*Cb[n2][k2]/Cb[nSam][ktt]/ktt);
			}
			index++;
			nSNP++;
		}
		if(z==rightzmax)
			LnL1 += (rightbound-interval[s][z-1]-nSNP)*log(1-theta*Table[n1t][z][0][0]);
		else
			LnL1 += (interval[s][z]-interval[s][z-1]-nSNP)*log(1-theta*Table[n1t][z][0][0]);
		nSNPtt += nSNP;
	}
	LnN += (rightbound+leftbound-nSNPtt-2)*log(1-theta*pHomo[n1t]);
	return (LnL1-LnN);
}
//------------------------------------------------------------------------------
void maximumLikelihood(int ntry)
//return ...
{
	int i,j,k,s,core,numSg;
	int **groupInfo;

	theta = 0.0;
	numSg = Sites.size();
	for(j=0; j<numSg; j++){	
		i = 0;
		for(k=0; k<allele.size(); k++) i += allele[k][j];
		theta += 2.0*i*(nSam-i);
	}
		
	theta = theta/seq_length/nSam/(nSam-1);
	//printf("theta= %f  \n", theta);
	groupInfo = new int*[numSg];
	for(i=0;i<numSg;i++) groupInfo[i]= new int[2];

	for(core=2; core<(numSg-1); core++){
			
		int n1 = 0;
		for(i=0; i<nSam; i++) n1 += allele[i][core];
		if( n1 >= minn1 && n1 <= maxn)
		{
			for(i=0;i<numSg;i++)
				for(j=0;j<2;j++) groupInfo[i][j] = 0;
				
			for(i=0;i<nSam;i++){
				if(allele[i][core]==1)
					for(j=0; j<numSg; j++) groupInfo[j][0] += allele[i][j];	
				else
					for(j=0; j<numSg; j++) groupInfo[j][1] += allele[i][j];
			}
	 		float L1 = -1e5; int s1 = 0;
			for(s=1; s<11; s++){
				float L1s = CalcLikelihood(n1,s,core,groupInfo);
				if(L1<L1s){	
					L1 = L1s; 
					s1 = s;
				}  
			}
			fprintf(outfile,"%d %ld %d %d %.0f %f\n",ntry,Sites[core],n1,nSam-n1,AlphaTable[s1],L1);	
		}
	}
	for (i=0;i<numSg;i++)
		delete [] groupInfo[i];
	delete [] groupInfo;
}
//------------------------------------------------------------------------------

int main(int argc, char** argv)
{
	int i,Ntry;
	float Rn;
	char Line[100];
	
	if(argc<3){
		printf("pls give input/output files' name");
		exit(-1);
	}
	
	max_alpha = 20000;
	for(i=3;i<argc;i++){
		if(!strcmp(argv[i], "-a")){
			if(i!=argc-1){
				if(argv[i+1][0]!='-') max_alpha = atof(argv[++i]);
			}
		}
	}
	
	if ( (inputdata = fopen( argv[1], "r" )) == NULL ){
			printf("Could not find input file\n");
			exit(-1);
		}
	if ( (outfile = fopen( argv[2], "w")) == NULL ){
        	printf("Could not open output file \n");
         	exit(-1);
        }
        
 	Rn=0.0; seq_length=0; Ntry=0;nSam=0;
	fgets(Line,100,inputdata);
	char *str = strtok(Line," ");
	str = strtok(NULL," ");
	nSam=atoi(str);
	str = strtok(NULL," ");
	Ntry = atoi(str);
	while(str != NULL){	
		if(strstr(str,"-r")){
			str = strtok(NULL," ");
			Rn=atof(str);
			str = strtok(NULL," ");
			seq_length = atoi(str);
		}
		str = strtok(NULL," ");
	}

	//printf("nSample=%d, Ntry=%d, Rn=%f, seq_length=%ld\n",nSam,Ntry,Rn,seq_length);
	
	rho = Rn/seq_length;
	minn1 = int(nSam*0.3+0.5);
	maxn = int(nSam*0.7+0.5);
	n1size = maxn-minn1+1;
	//start calculating
	CbTable();
	subPop_freqSpec();
	InputProbability();
	make_interval();
	
	fprintf(outfile,"REP Position nd na alpha L\n");
	for(i=0;i<Ntry;i++){
		readmsout();
		maximumLikelihood(i);
	}
	
	fclose(inputdata);
	fclose(outfile);
	return 0;
}

