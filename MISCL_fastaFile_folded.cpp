#include "MISCL.h"

int minncore,maxmiss_core,*minn2,*n2size;
float **pHomo;
double *****Table;
vector<char *> allele;
FILE *outfile,*infile;
//------------------------------------------------------------------------------
void CbTable()
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
void InputProbability()
{
	int i,j,n,m,n1,n2,ntt,z;
	
//for calculating Pn[n1][n2][0][0], it is calculated by (1-probability of being SNP)
	
	pHomo = new float*[n1size];
	for(n=0;n<n1size;n++) pHomo[n] = new float[n2size[n]];
	

	for(int n=0; n< n1size; n++)
	for(int m=0; m< n2size[n]; m++){
		n1=n+minn1;
		n2=m+minn2[n];
		ntt=n1+n2;
		pHomo[n][m] = 0.0;
		for(i=0; i<=n1; i++)
		for(j=0; j<=n2; j++)
			if( (i+j)!=0 && (i+j)!=ntt )
		 		pHomo[n][m] += Cb[n1][i]*Cb[n2][j]/Cb[ntt][i+j]*subP[ntt][i+j];
	}
	
//P[n1][n2][k1][k2] for selective sweep model
	//assign space
	Table = new double****[n1size];
	for(n=0; n<n1size; n++){
		Table[n] = new double***[n2size[n]];
		for(m=0; m<n2size[n]; m++){
			Table[n][m] = new double**[100];
			n1=n+minn1;
			n2=m+minn2[n];
			for(z=0; z< 100; z++){
				Table[n][m][z] = new double*[n1+1];
				for(i=0; i<= n1; i++)
					Table[n][m][z][i] = new double[n2+1];
			}
		}
	}

	//calculating
	for(n=0; n<n1size; n++)
	for(m=0; m<n2size[n]; m++){
		n1=n+minn1;
		n2=m+minn2[n];
		for(z=0; z< 100; z++){
			Table[n][m][z][0][0] = 0.0;
			for(i=0; i<= n1; i++)
			for(j=0; j<= n2; j++)
			if( (i+j)!=0 && (i+j)!= (n1+n2) ){
				Table[n][m][z][i][j] = ProbNielsen(n1,n2,z,i,j);
				Table[n][m][z][0][0] += Table[n][m][z][i][j];
			}
		}
	}

}
//--------------------------------------------------------------------------------
void readrawdata()
{
	int i, j,k,Lseg,l,max_miss;
	char Line[500];
	vector<char *>chro_segment;
	   
	long count=1; char c=fgetc(infile);
	if(c == '>'){
		while(c != '\n') c=fgetc(infile);
		c = fgetc(infile);
	}
	while(c != '\n' && !feof(infile)){
		c = fgetc(infile);
		count++;
	}
	
	seq_length = count;
	//printf("seq_length: %ld\n",seq_length);
	Lseg = 1000000;
	int maxii = seq_length/Lseg+1;
	int numSg = 0;
	for(int ii=0;ii<maxii;ii++){
		
		for(j=0;j<chro_segment.size();j++)
    		delete[](chro_segment[j]);
  		chro_segment.clear();
  		
  		i=0; long bgseg = ii*Lseg; long enseg = bgseg+Lseg;
		if(enseg>seq_length) enseg = seq_length;
  		
		rewind(infile);
 		
		while(!feof(infile)){
			char *pl= new char[Lseg];
			memset(pl,'0',Lseg);
			l=0; char c=fgetc(infile);
			if(c == '>'){
				while(c != '\n') c=fgetc(infile);
				c = fgetc(infile);
			}
			while(c != '\n' && !feof(infile)){
				if(l >= bgseg && l<enseg) pl[l-bgseg]=c;
				c=fgetc(infile);
				l++;
			}
			chro_segment.push_back(pl);
			i++;
		}
		
		if(ii==0){
			nSam = i-1; //printf("\n nSam: %d %ld",nSam,chro_segment.size());
			max_miss = int(nSam*0.5+0.5);
		}
		for(i=0;i<(enseg-bgseg);i++){
		 	int ht=1;
 			int checksnp = 0;
 			int nm = 0;
 			char fb = chro_segment[0][i];
 		
 			if(fb =='N'){
 				for(j=1;j<nSam;j++){
			 		if(chro_segment[j][i] != 'N'){
		 				fb = chro_segment[j][i];
	 					ht = j+1;
	 					break;
	 				}
			 	}
			 	nm = ht-1;
 			}
 			
 			for(j=ht;j<nSam;j++){
			 	if(chro_segment[j][i] == 'N') nm++;
			 }
			 
 			if(fb == 'N' || nm>max_miss){
			 	char *pl= new char[1];
			 	pl[0]='M';
			 	allele.push_back(pl);
			 } 
 			else{
		 		for(j=ht;j<nSam;j++){
 					if(chro_segment[j][i] != 'N' && chro_segment[j][i] != fb){
			 			checksnp = 1;
			 			break;
			 		}
 				}
		 	 
				 if(checksnp==1){
				 	char *pl= new char[nSam+2];
				 	memset(pl,'0',nSam+2);
				 	pl[0]='P';
				 	pl[1]=fb;
				 	for(k=2;k<(nSam+2);k++) pl[k]= chro_segment[k-2][i];
		 			allele.push_back(pl);
 					numSg++;
 				  }
 				 else{
 			 		char *pl= new char[1];
				 	pl[0]='H';
				 	allele.push_back(pl);
	 			  }
 			  }
		}
	}
	
	float Sum = 0.0;
	//printf("\n numSg %d, allele.size %ld, nSam %d, L %ld\n",numSg,allele.size(),nSam,seq_length);
	fclose(infile);
}
//------------------------------------------------------------------------------

void CalcLnL(int s,long Pos,int nSam1,int nSam2, double *LnL)
{
	int z,j,n1t,n2t,n1tE,n2tE,ntt,ktt,n1Table,n2Table,n1TableE,n2TableE;
	long goleft, goright;
	
	LnL[0] = 0.0; LnL[1]=0.0; LnL[2]=0.0;
	n1Table = nSam1-minn1; n2Table = nSam2-minn2[n1Table];
	n1TableE = nSam2-minn1; n2TableE = nSam1-minn2[n1TableE];
	goleft=Pos-1; goright=Pos+1;
	
	for(z=1; z<100 ; z++){
			int nrHSite = 0;
			while(goright<seq_length && (goright-Pos)<interval[s][z]){
				if(allele[goright][0]=='P'){
					int k1 = 0; int k2 = 0; int n1 = 0; int n2 = 0;
					for(j=2;j<nSam+2;j++){
						if(allele[Pos][j]==allele[Pos][1]){
							if(allele[goright][j] != 'N'){
								n1++;
								if(allele[goright][j]==allele[goright][1]) k1++;
							}
						}else if(allele[Pos][j] != 'N') {
							if(allele[goright][j] != 'N'){
								n2++;
								if(allele[goright][j]==allele[goright][1]) k2++;
							}
						}
					}
					n1t = n1-minn1; n1tE = n2-minn1;
					if(n1>=minn1 && n2>=minn2[n1t]){
						n2t = n2-minn2[n1t]; n2tE = n1-minn2[n1tE]; ntt =n1+n2; ktt = k1+k2;
						if( ktt==0 || ktt==ntt ){
							LnL[0] += log(1.0-theta*pHomo[n1t][n2t]);
							LnL[1] += log(1.0-theta*Table[n1t][n2t][z][0][0]);
							LnL[2] += log(1.0-theta*Table[n1tE][n2tE][z][0][0]);
						}
						else{
							LnL[0] += log(theta*Cb[n1][k1]*Cb[n2][k2]/Cb[ntt][ktt]*(subP[ntt][ktt]+subP[ntt][ntt-ktt])/2);
							LnL[1] += log(theta*(Table[n1t][n2t][z][k1][k2]+Table[n1t][n2t][z][n1-k1][n2-k2])/2 );
							LnL[2] += log(theta*(Table[n1tE][n2tE][z][k2][k1]+Table[n1tE][n2tE][z][n2-k2][n1-k1])/2 );
						}
					}
				}
				else if(allele[goright][0] == 'H') nrHSite++;
				goright++;
			}
		
			int nlHSite = 0;
			while(goleft > 0 && (Pos-goleft)<interval[s][z]){
				if(allele[goleft][0]=='P'){
					int k1 = 0; int k2 = 0; int n1 = 0; int n2 = 0;
					for(j=2;j<nSam+2;j++){
						if(allele[Pos][j]==allele[Pos][1]){
							if(allele[goleft][j] != 'N'){
								n1++;
								if(allele[goleft][j]==allele[goleft][1]) k1++;
							}
						}else if(allele[Pos][j] != 'N') {
							if(allele[goleft][j] != 'N'){
								n2++;
								if(allele[goleft][j]==allele[goleft][1]) k2++;
							}
						}
					}
					n1t = n1-minn1; n1tE = n2-minn1;
					if(n1>=minn1 && n2>=minn2[n1t]){
						n2t = n2-minn2[n1t]; n2tE = n1-minn2[n1tE]; ntt =n1+n2; ktt = k1+k2;
						if( ktt==0 || ktt==ntt  ){
							LnL[0] += log(1.0-theta*pHomo[n1t][n2t]);
							LnL[1] += log(1.0-theta*Table[n1t][n2t][z][0][0]);
							LnL[2] += log(1.0-theta*Table[n1tE][n2tE][z][0][0]);
						}
						else{
							LnL[0] += log(theta*Cb[n1][k1]*Cb[n2][k2]/Cb[ntt][ktt]*(subP[ntt][ktt]+subP[ntt][ntt-ktt])/2);
							LnL[1] += log(theta*(Table[n1t][n2t][z][k1][k2]+Table[n1t][n2t][z][n1-k1][n2-k2])/2 );
							LnL[2] += log(theta*(Table[n1tE][n2tE][z][k2][k1]+Table[n1tE][n2tE][z][n2-k2][n1-k1])/2 );
						}
					}
					
				}
				else if(allele[goleft][0] == 'H') nlHSite++;
				goleft--;
			}
		LnL[0] += (nrHSite+nlHSite)*log(1.0-theta*pHomo[n1Table][n2Table]);
		LnL[1] += (nrHSite+nlHSite)*log(1.0-theta*Table[n1Table][n2Table][z][0][0]);
		LnL[2] += (nrHSite+nlHSite)*log(1.0-theta*Table[n1TableE][n2TableE][z][0][0]);
	}
	
}
//------------------------------------------------------------------------------
void Likelihood(long Pos)
{
	int i,n1,n2,nmiss;
	double *LnL,deltaLs[2];
	
	LnL = new double[3];
	n1 = 0; nmiss = 0;
	for(i=2; i<nSam+2; i++)
		if(allele[Pos][i]== 'N') nmiss += 1;
		else if(allele[Pos][i]==allele[Pos][1]) n1 += 1;
	n2 = nSam-n1-nmiss;
	if( n1 >= minncore && nmiss <= maxmiss_core && n2 >= minncore )
	{
		
			double deltaL = -1e5; int sc=0; int nd=0; int na=0;
			
			for(int s=1;s<11;s++){
				LnL[0] = 0.0; LnL[1]=0.0; LnL[2]=0.0;
				CalcLnL(s,Pos,n1,n2,LnL);
				deltaLs[0] = LnL[1]-LnL[0];
				deltaLs[1] = LnL[2]-LnL[0];
				
				if(deltaLs[0]>deltaLs[1]){
					if(deltaL < deltaLs[0]){
						deltaL = deltaLs[0];
						sc=s;
						nd=n1;
						na=n2;
					}
				}else{
					if(deltaL < deltaLs[1]){
						deltaL = deltaLs[1];
						sc=s;
						nd=n2;
						na=n1;
					}
				}
			}
			
			fprintf(outfile, "F %ld %d %d %.0f %f\n",Pos,nd,na,AlphaTable[sc],deltaL);
	}
}
//--------------------------------------------------------------------------------------
int main(int argc, char** argv)
{
	int i;
	
	if(argc<3){
		printf("pls give input/output files' name");
		exit(-1);
	}
		
	if ( (infile = fopen( argv[1], "r")) == NULL )
        {
                printf("Could not open input file \n");
                exit(-1);
        }
     
	if ( (outfile = fopen( argv[2], "w")) == NULL )
        {
                printf("Could not open output file \n");
                exit(-1);
        }
	
	max_alpha = 20000;
	for(i=1;i<argc;i++){
		if(!strcmp(argv[i], "-a")){
			if(i!=argc-1){
				if(argv[i+1][0]!='-') max_alpha = atof(argv[++i]);
			}
			else{
				printf("\n Population size is missing\n");
				exit(-1);
			}
		}
		if(!strcmp(argv[i], "-r")){
			if(i!=argc-1){
				if(argv[i+1][0]!='-') rho = atof(argv[++i]);
			}
			else{
				printf("\nPopulation recombination rate is missing\n");
				exit(-1);
			}
		}
	}
	if(rho <= 0 || max_alpha<100){
		printf("\nmaximum alpha and recombination rate can not be smaller than 0\n");
		exit(-1);
	}
	
	readrawdata();
	theta =0;
	for(i=0;i<allele.size();i++){
		if(allele[i][0]=='P'){
			int nm = 0; int ni = 0;
			for(int j=2; j<nSam+2; j++){
				if(allele[i][j]== 'N') nm += 1;
				else if(allele[i][j]==allele[i][1]) ni += 1;
			}
			nm = nSam-nm;
			theta += 2.0*ni*(nm-ni)/nm/(nm-1);
		}
	}
	theta = theta/seq_length;
	//printf(" theta: %f ",theta);
	
	minncore = int(nSam*0.3+0.5); maxmiss_core = int(nSam*0.1+0.5);
	minn1 = int(nSam*0.2+0.5); maxn = int(nSam*0.7+0.5);
	n1size = maxn-minn1+1; //printf("\nmaxalpha: %f,rho: %f,nn1size %d, ",max_alpha,rho,n1size);
	minn2 = new int[n1size]; n2size = new int[n1size];
	for(i=0;i<n1size;i++){
		int n1 = i+minn1;
		minn2[i] = 0.6*nSam>n1 ? int(0.8*nSam-n1+0.5) : int(0.2*nSam+0.5); //printf("\nminn2[%d] %d ",i,minn2[i]);
		int maxn2 = 0.3*nSam>n1 ? int(0.7*nSam+0.5) : (nSam-n1); //printf("%d ",maxn2);
		n2size[i] = maxn2-minn2[i]+1;
	}
	
	CbTable();
	subPop_freqSpec();
	InputProbability();
	make_interval();

	fprintf(outfile,"T Position nd na alpha L\n");
	for(i=0;i<allele.size();i++)
		if(allele[i][0]=='P') Likelihood(i);
	
	fclose(outfile);
	return 0;
}
