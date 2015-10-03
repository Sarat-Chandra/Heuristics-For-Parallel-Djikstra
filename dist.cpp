//General DS for Shortest Path computation
#include <LEDA/graph/graph.h>
#include <LEDA/graph/node_pq.h>

//Used for co-ordinate generation
#include <LEDA/geo/point.h>
#include <LEDA/geo/random_point.h>

//LEDA Namespace is used generally
using namespace leda; 
//C++ Headers n Libs
#include<iostream>
#include<string.h>
#include<stdio.h>
//#include<conio.h>

using std::cin;
using std::cout;
using std::endl;

#include<iomanip>
#include<stdlib.h>
int main()
{
FILE *fp,*fp1;
int numb=0;
int srch=0;
char num[200],place[100][20];
int i=0,k=1;
char line[200],line1[100];
fp=fopen("dist1.txt","r");
//clrscr();
while(!feof(fp))
{
fgets(line,100,fp);
cout<<line;
}
fclose(fp);
cout<<"Enter number:"<<endl;
cin>>numb;
cout<<numb<<endl;
fp1=fopen("dist1.txt","r");
cout<<srch;
while(!feof(fp1) && (srch<numb))
{
srch++;
fgets(line1,100,fp1);
}
cout<<srch;
cout<<line1;
fclose(fp1);
/*
for(k=0;k<=i;k++)
printf("%s\n",num[k]);*/
//printf("%s",num[3]);
//getch();
return 0;
}

