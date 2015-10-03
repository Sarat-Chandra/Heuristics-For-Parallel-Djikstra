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

using std::cout;
using std::endl;

#include<iomanip>
#include<stdlib.h>

int main()
{
FILE *fp=fopen("ad.txt","r");
graph G;
int cnt=0;
char ch[20];
complete_graph(G,100);
edge_array<double> ea(G);
edge e;
node src, dest;
forall_edges(e,G) {
src=source(e);
dest=target(e);
if(!feof(fp))
{
fscanf(fp,"%s",ch);
G.print_node(src);
cout<<"-->";
G.print_node(dest);
cout<<"\t"<<ch<<endl;
ea[e]=(double)atol(ch);
cnt++;
}
}
cout<<cnt;
fclose(fp);
return 0;
}
