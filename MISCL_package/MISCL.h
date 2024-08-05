
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string.h>
using namespace std;

int nSam;
int minn1;
int maxn;
int n1size;
long interval[11][100];
long seq_length;
float theta;
float max_alpha;
float AlphaTable[11];
float rho;
float **Cb;
float **subP;

vector<long> Sites;

void CbTable();
void subPop_freqSpec();
double ProbNielsen(int nSam1,int nSam2,int Ze,int k1,int k2);
void make_interval();
