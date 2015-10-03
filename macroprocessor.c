#include<stdio.h>
#include<conio.h>
void main()
{
char *str,*name,*arg1,*arg2,*arg3;
int st=0x00,end=0x00;
char *str1,*str2;
int count=1;
int ptrno=1;
FILE *fmac,*fdef,*fnam,*farg;
clrscr();
fmac=fopen("myin.txt","r");
fdef=fopen("deftab.txt","w+");
fnam=fopen("namtab.txt","w+");
farg=fopen("argtab.txt","w+");
while(!feof(fmac))
  {
	fgets(str,79,fmac);
	sscanf(str,"%s%s",str1,str2);
	if(strcmp(str2,"macro")==0)
		{
		   fprintf(fnam,"%s\t%x\t",str1,ptrno);
		   ptrno++;
			do
			{
			fprintf(fdef,"%s",str);
			fgets(str,79,fmac);
			sscanf(str,"%s%s",str1,str2);
			ptrno++;
			}while(strcmp(str2,"mend")!=0);
		   fprintf(fdef,"*\tmend\t*\n\n");
		   fprintf(fnam,"%x\n",ptrno-1);
		}
	rewind(fnam);
	while(!feof(fnam))
		{
		fscanf(fnam,"%s%x%x",name,&st,&end);
		printf("%s\n",str);
		if(strcmp(str2,name)==0)
		{
		sscanf(str,"%s%s%s%s%s",str1,str2,arg1,arg2,arg3);
		fprintf(farg,"Arguments for invocation(%s) %d:\t%s\t%s\t%s\n\n",name,count++,arg1,arg2,arg3);
		printf("%s%x%x\n",name,st,end);
		}
		}
  }
fcloseall();
clrscr();
gotoxy(15,15);
printf("REFER FILES\tdeftab.txt\n\t\t\t\targtab.txt\n\t\t\t\tnamtab.txt\n");
getch();
}
