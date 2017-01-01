#include<iostream>
#include<cmath>
#include<fstream>
#include<stdlib.h>

using namespace std;

int main()
{
	ifstream fin;
	int lumi[2000];
	int N=0;
	int temp;
	int flag;
	int i;
	
	fin.open("lumi_list.txt");
	while(!fin.eof())
	{
		fin>>lumi[N++];
	}
	N--;
	fin.close();

	i=87;
	while(i<=1880)
	//for(int i=87; i<=216; i++)
	{	
		if(i==217) i=219;
		else if (i==588) i=589;
		else if (i==681) i=683;
		else if (i==1201) i=1217;
		else if (i==1411) i=1413;
		else if (i==1849) i=1880;
	
		flag=0;
		for(int j=0; j<N; j++)
		{
			if(i==lumi[j]) {
				flag=1;
				break;
			}
		}
		if(flag==0)
			std::cout<<i<<"\n";
		i++;
	}
	
	return 0;
}
