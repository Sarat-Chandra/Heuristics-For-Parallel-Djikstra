#include<stdio.h>
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
printf("%s",line);
}
fclose(fp);
printf("Enter number:\n");
scanf("%d",&numb);
printf("%d\n", numb);
fp1=fopen("dist1.txt","r");
printf("%d\n",srch);
while(!feof(fp1) && (srch<numb))
{
srch++;
fgets(line1,100,fp1);
}
printf("%d\n",srch);
printf("%s",line1);
fclose(fp1);
/*
for(k=0;k<=i;k++)
printf("%s\n",num[k]);*/
//printf("%s",num[3]);
//getch();
return 0;
}

