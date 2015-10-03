#include<iostream>
//#include<conio.h>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
using namespace std;
int main()
{
char lp[100],ss[10]=" ",c;
int i=0,j=0,col_num=1,u=0,v=0,dist=0;
int val[50][50];
FILE *fp;
fp=fopen("adj_list.txt","r");
for(u=1;u<=17;u++){
for(v=1;v<=17;v++){
val[u][v]=0;
cout<<val[u][v]<<"\t";}
cout<<endl;}
u=1;
v=1;
while(!feof(fp)){//while starts
if((c=fgetc(fp))!='\n'){//if1 starts
	if(c!='\t'){//if2 starts
	//if(col_num==0){
	//col_num++;
	//cout<<"Space"<<endl;}
	ss[i++]=c;}//if2 ends
	else{//else2 starts
	ss[i]=' ';
	if(col_num==1){//if3 starts
	col_num++;
	u=atoi(ss);
	}//if3 ends
	else if((col_num%2)==0)//if3
	{dist=atoi(ss);
	col_num++;}//if3
	else
	{//if3
	if((col_num%2)==1){//if3
	v=atoi(ss);
	col_num++;
	val[u][v]=dist;}//if3
	}//if3
	cout<<ss<<endl;
	j++; 
	i=0;}//else2 ends
}//if1 ends
else
{col_num=1;}
}
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
for(u=1;u<=17;u++){
for(v=1;v<=17;v++){
//val[u][v]=0;
cout<<val[u][v]<<"\t";}
cout<<endl;}
}
