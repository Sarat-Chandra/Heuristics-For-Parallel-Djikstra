#include<LEDA/graph/graph.h>
#include<LEDA/graph/node_pq.h>
#include<LEDA/graph/node_list.h>
#include<stdio.h>

//Custom Values
#define MAX_DOUBLE   21474
#define MIN_DOUBLE  -21474
#define MAXSAMPLERUNS 10

//Assuming a 6 x 6 rectangular grid
#define GRIDCOUNT 36
#define MAXCOORD 32767
#define MINCOORD -32768

#if defined(LEDA_NAMESPACE)
using namespace leda;
using namespace std;
#endif
#if defined(LEDA_STD_IO_HEADERS)
//using std::cout;
//using std::endl;
#endif

void level_part(graph &G, int number, node_array<double>& node_num, FILE *fp)
{
node n;
edge e;
node tg;
node_list perm;
forall_nodes(n,G)
{
if(node_num[n]<=number)
{
if(!perm.member(n))
perm.append(n);
forall_out_edges(e,n)
{
tg=G.target(e);
if(node_num[tg]>number){
if(!perm.member(tg)){
perm.append(tg);}
//G.print_edge(e);
//cout<<endl;}
}
}
}
while(!perm.empty()){
node p;
int num;
p=perm.pop();
num=node_num[p];
fprintf(fp,"%d\t",num);
}
fprintf(fp,"\n");
}

int main(){
graph G;
FILE *fp;
fp=fopen("Levels.txt","w");
int numb=0;
int lvl_num=1;
edge e;
node n;
random_graph(G,100,rand_int(5,10) * 100 , true, true, true);
Make_Connected( G );
G.make_directed();
node_array<double> node_num(G);
forall_nodes(n,G)
{
	node_num[n]=numb++;
}
for(lvl_num=1;lvl_num<=3;lvl_num++){

cout<<"LEVEL NUMBER:\t"<<lvl_num<<endl;
numb=numb+100*lvl_num/3;
level_part(G,numb,node_num,fp);
numb=0;}
}
