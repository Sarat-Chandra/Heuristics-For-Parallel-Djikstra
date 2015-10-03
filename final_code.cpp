//General DS for Shortest Path computation
#include <LEDA/graph/graph.h>
#include <LEDA/graph/node_pq.h>

//Used for co-ordinate generation
#include <LEDA/geo/point.h>
#include <LEDA/geo/random_point.h>

//LEDA Namespace is used generally
using namespace leda;

#include<iostream>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
using namespace std;

int dijk_orig(const graph &G,const edge_array<int>& cost, node s, node t) {
  node_pq<double> PQ(G);
  node v; edge e;
  node_array<double> dist(G);
node_list perm;
  forall_nodes(v,G) {
    if (v == s) dist[v] = 0; else dist[v] = MAXDOUBLE;
      PQ.insert(v,dist[v]);
  }

  int vertexVisitCnt = 0;

  while ( !PQ.empty() ) {
    node u = PQ.del_min();
	//G.print_node(u);
    //Oops: All the vertices are not reachable from given source
    if( dist[u] == MAXDOUBLE ) {
        PQ.clear();
        break;
    }

    vertexVisitCnt++;
    //Dijk Mod: If target is reached break
    if( /*u == t ||*/ perm.member(u) ) {
        break;
    }
    //G.print_node(u);
	perm.append(u);
    forall_out_edges(e,u) {
        v = target(e);        
        double c = dist[u] + cost[e];

        if ( c < dist[v] ) {
	  G.print_node(v);
          PQ.decrease_p(v,c);  dist[v] = c;  
        }
    } 
  }
  node u;
  cout<<"DISTANCE VECTOR:\n";
  forall_nodes(v,G){
  cout<<dist[v]<<"\n";}
  return vertexVisitCnt;
}








int main()
{
graph G;
char cities[100];
edge e;
node nk;
node st,tg;
char lp[100],ss[10]=" ",c;
int i=0,j=0,col_num=1,u=0,v=0,dist=0,node_num=1,src,tgt;
int val[50][50];
FILE *fp,*fp_cities;
fp_cities=fopen("cities.txt","r");
while(!feof(fp_cities))
{
fgets(cities,100,fp_cities);
cout<<cities<<endl;}
fclose(fp_cities);
fp=fopen("adj_list.txt","r");
for(u=1;u<=17;u++){
for(v=1;v<=17;v++){
val[u][v]=0;
//cout<<val[u][v]<<"\t";
}
//cout<<endl;
}
u=1;
v=1;
complete_graph(G,17);
edge_array<int> d(G);
node_array<int> number(G);
forall_edges(e,G){
d[e]=0;}
forall_nodes(nk,G){
number[nk]=node_num++;}
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
	val[u][v]=dist;
	forall_edges(e,G){
	st=G.source(e);
	tg=G.target(e);
	if((number[st]==u) && (number[tg]==v))
	{
	d[e]=dist;}}
	}//if3
	}//if3
//	cout<<ss<<endl;
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

//cout<<"EDGE DISTANCES"<<endl;
forall_edges(e,G){
if(d[e]!=0){
//st=G.source(e);
//tg=G.target(e);
//G.print_node(st);
//cout<<"-->";
//G.print_node(tg);
//cout <<"\t";
//cout<<d[e]<<"\n";
}
else
G.del_edge(e);
}

int vertexCnt_orig;
float runTime_orig=0;
cout<<"ENTER SOURCE:";
cin>>src;
cout<<"ENTER TARGET:";
cin>>tgt;
forall_nodes(nk,G){
if(number[nk]==src)
st=nk;
if(number[nk]==tgt)
tg=nk;
}
//G.print_node(st);
//G.print_node(tg);

float T = used_time();
  vertexCnt_orig= dijk_orig(G, d, st, tg);
          runTime_orig = used_time( T );//DIJKSTRA(G,G.first_node(),cost,dist);
  cout << "\n\nThe shortest path computation took " << 
          runTime_orig << " seconds.\n\n";
forall_edges(e,G){
st=G.source(e);
G.print_node(st);
tg=G.target(e);
G.print_node(tg);
cout<<d[e]<<"\n";}
}
