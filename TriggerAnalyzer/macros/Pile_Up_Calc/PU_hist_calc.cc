#include<iostream>
#include<stdio.h>
#include<math.h>
#include<fstream>
#include<string.h>
#include<stdlib.h>

//#include "TFile.h"
//#include "TH1.h"

using namespace std;

int main()
{
	//TFile *f = new TFile("output.root");
	//TH1D * h1 = (TH1D*)f->Get("pileup");	


	ifstream fin;
	ofstream fout;
	int n;
	double PU_MC[100], PU_Data[100], PU_weight[100];
	int x[100];
	double sum;
	double norm;
	
	n=75;
	sum=0;
	
	fin.open("PU_MC.txt");
	for(int i=0; i<n; i++)
	{
		x[i]=i;
		fin>>PU_MC[i];
		//sum+=PU_MC[i];
	}
	fin.close();

	fin.open("PU_Data.txt");
	for(int i=0; i<n; i++)
	{
		fin>>PU_Data[i];
		sum+=PU_Data[i];
	}
	fin.close();
	
	norm=sum;
	sum=0;
	fout.open("PU_weights.txt");
	
	for(int i=0; i<n; i++)
	{
		PU_Data[i] = PU_Data[i]/norm;	
		sum+=PU_Data[i];
		if(PU_MC[i]==0)
			PU_weight[i]=0;
		else	
			PU_weight[i] = PU_Data[i]/PU_MC[i];
		fout<<x[i]<<"   "<<PU_weight[i]<<"\n";
	}
	
	fout.close();




	return 0;
}
