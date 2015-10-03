#include<iostream>
//#include<conio.h>
#include<stdio.h>
using namespace std;
int main()
{
FILE *fp,*fp1;
int i,j;
int array[36][36];
char line[50];
int c1,c2,c3;
int cnt=0;
//clrscr();
cout<<"initializing"<<endl;
for(i=1;i<=35;i++)
{
for(j=1;j<=35;j++)
{
array[i][j]=0;
//printf("%d\t",array[i][j]);
}
}
cout<<"initialization complete"<<endl;
fp=fopen("distance_number+9.txt","r");
if (fp==NULL){
puts("ERROR");}
else{
cout<<"file ready for reading"<<endl;
while(!feof(fp)){
fscanf(fp,"%d%d%d",&c1,&c2,&c3);
array[c1][c2]=c3;
array[c2][c1]=c3;
printf("%d\t%d\t%d\n",c1,c2,c3);}}
fclose(fp);
fp1=fopen("adj_list.txt","w");
for(i=1;i<=35;i++){
fprintf(fp1,"%d\t",i);
for(j=1;j<=35;j++){
if(array[i][j]!=0)
{
fprintf(fp1,"%d\t%d\t",array[i][j],j);
cnt++;                }
}
fprintf(fp1,"\n");
}
fclose(fp1);
printf("Total count: %d\n",cnt);
return 0;
}
