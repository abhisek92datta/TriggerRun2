#include<iostream>
#include<cmath>
#include<fstream>
#include<stdlib.h>

using namespace std;

int max(int array[1000], int size){

    int max=array[0];
    int j=0;
    for(int i=0; i<size; i++){
        if(array[i]>=max){
            max = array[i];
            j=i;
        }
    }
    //std::cout<<max<<"  "<<j<<"\n";
    return max;
}


int main()
{

    system("awk '{if($1>=280919 && $1<=284044){print $1, $2}}' json_test.txt > test1.txt");

    int n=0;
    double run[1000], start[1000], end[1000], omit[1000][1000];
    double prev, next;
    int omit_size[1000];
    int omit_max;


    char temp[1000];
    char temp2[1000];
    int j=0;
    int k=0;
    int m=0;

    ifstream fin;
    ofstream fout;
    fin.open("test1.txt");
    while(!fin.eof()){
        fin>>temp>>temp;
        n++;
    }
    n--;
    fin.close();
    fin.open("test1.txt");
    fout.open("list.txt");

    for(int i=0; i<n; i++)
        omit_size[i]=0;

    for(int i=0; i<n; i++){
        fin>>run[i]>>temp;
        j=0;
        k=2;
        while(temp[k]!=','){
            temp2[j++]=temp[k];
            k++;
        }
        temp2[j]='\0';
        start[i] = atof(temp2);
        j=0;
        k++;

        while(temp[k]!='\0'){
            //std::cout<<temp[k]<<"  ";

            if(temp[k]!=']' && temp[k]!=','){
                temp2[j++]=temp[k];
                k++;
            }

            else if(temp[k]==']' && temp[k+1]==']'){
                temp2[j]='\0';
                end[i] = atof(temp2);
                //std::cout<<end[i]<<"\n";
                break;
            }
            else if(temp[k]==']' && temp[k+1]!=']'){
                temp2[j]='\0';
                prev = atof(temp2);
                j=0;
                k+=3;
            }
            else if(temp[k]==','){
                temp2[j]='\0';
                next = atof(temp2);
                j=0;
                k++;
                for(double z=prev+1; z<next; z++)
                    omit[i][omit_size[i]++]=z;
            }



        }
        //std::cout<<"\n";
        omit[i][omit_size[i]++] = 0;
    }
    fin.close();

    omit_max = max(omit_size,n);

    for(int i=0; i<n; i++){
        fout<<run[i]<<", ";
    }
    fout<<"\n\n";
    for(int i=0; i<n; i++){
        fout<<start[i]<<", ";
    }
    fout<<"\n\n";
    for(int i=0; i<n; i++){
        fout<<end[i]<<", ";
    }
    fout<<"\n\n";

    for(int i=0; i<n; i++){
        //std::cout<<omit_size[i]<<"\n";
        fout<<"{ ";
        for(int ii=0; ii<omit_size[i]; ii++){
            fout<<omit[i][ii]<<" ";
            if(ii!=(omit_max-1))
                fout<<", ";
        }
        for(int ii=omit_size[i]; ii<omit_max; ii++){
            fout<<"0 ";
            if(ii!=(omit_max-1))
                fout<<", ";
        }
        fout<<"},\n";
    }




    cout<<n<<"\n";
    fout.close();

    system("rm -rf test1.txt");

	return 0;
}
