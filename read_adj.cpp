#include<iostream>
//#include<conio.h>
#include<stdio.h>
#include<string.h>
using namespace std;
int main()
{
char lp[100],ss[10]=" ",c;
int i=0,j=0,k=0;
FILE *fp;
fp=fopen("adj_list.txt","r");
while(!feof(fp)){
if((c=fgetc(fp))!='\n'){
if(c!='\t'){
k++;
//cout<<"Space"<<endl;}
ss[i++]=c;}
else{
ss[i]=' ';
cout<<ss<<endl;
j++; 
i=0;}
}}
/*while(!feof(fp))
{
fgets(lp,100,fp);

while(true){
sscanf(lp,"%s",ss);
if(strcmp(ss,NULL)==0)
break;
cout<<ss<<endl;
}
cout<<lp<<endl;
}*/
/*
while(!feof(fp)){
fscanf(fp,"%s",lp);
if(strcmp(lp,"\n")==0)
cout<<"EOL reached"<<endl;}*/
//getch();
}
