#include<iostream>
#include<cmath>
#include<fstream>
#include<stdlib.h>

using namespace std;

int main()
{
	ifstream fin;
    char file[1000], a[100];
    int j=0;
    double list[1000];
    int n=0;
    int flag;

    fin.open("lumiF.txt");
    while(!fin.eof()){
        fin>>file;
        j=0;
        for(int i=17; file[i]!='\0'; i++){
            if(file[i]=='.')
                break;
            a[j++] = file[i];
        }
        a[j]='\0';
        list[n] = atof(a);
        n++;
    }
    n--;
    fin.close();

    for(double k=1; k<=904; k++){
        flag=0;
        for(int i=0; i<n; i++){
            if(k==list[i]){
                flag = 1;
                break;
            }
        }
        if(flag==0)
            std::cout<<k<<"\n";
    }


	return 0;
}
