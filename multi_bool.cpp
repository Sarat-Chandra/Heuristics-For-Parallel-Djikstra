#include<LEDA/graph/graph.h>
#include<LEDA/graph/node_pq.h>
#include<LEDA/graph/node_list.h>
#include<stdio.h>
#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif
#if defined(LEDA_STD_IO_HEADERS)
using std::cout;
using std::endl;
#endif

int DIJKSTRA1(const graph &G, node s, const edge_array<double> &cost, node_array<double> &dist,edge_array<bool> &b)
{
node_pq<double> PQ(G);
node_pq<double> rPQ(G);
node_array<double> rdist(G);
//node_array<double><double> try1(G);
node_list perm;
int vvc=0;

node v;
edge e;

node t=G.choose_node();
cout<<"\n\nTarget is:";
G.print_node(t);

forall_nodes(v,G)
{
if(v==s)
dist[v]=0;
else
dist[v]=MAXDOUBLE;
PQ.insert(v,dist[v]);
}

while(!PQ.empty()){
node u=PQ.del_min();
if(u==t || perm.member(u))
{
cout<<"\n\nTarget reached in forward search"<<endl;
cout<<"common node:";
G.print_node(u);
break;
}
perm.append(u);
forall_out_edges(e,u)
{
v=target(e);
double c=dist[u] + cost[e];
if(c<dist[v] && b[e]==1)
{
vvc++;
PQ.decrease_p(v,c); dist[v]=c;}
}

}
return vvc;
}


void DIJKSTRA(const graph &G, node s, const edge_array<double> &cost, node_array<double> &dist,edge_array<bool> &b)
{
node_pq<double> PQ(G);
node_pq<double> rPQ(G);
node_array<double> rdist(G);
node_list perm;

node v;
edge e;

node t=G.choose_node();

forall_nodes(v,G)
{
if(v==s)
dist[v]=0;
else
dist[v]=MAXDOUBLE;
PQ.insert(v,dist[v]);
}

while(!PQ.empty()/* && !rPQ.empty()*/){
node u=PQ.del_min();
if(u==t || perm.member(u))
{
break;
}

perm.append(u);
forall_out_edges(e,u)
{
v=target(e);
double c=dist[u] + cost[e];
if(c<dist[v])
{
b[e]=1;
PQ.decrease_p(v,c); dist[v]=c;}
}


}

}

int main()
{
//FILE *fp=fopen("ol1.txt","w+");
//FILE *fp1;
node nd;
char str[15];
int n;//=read_int("number of nodes=");
int m;//=read_int("number of edges=");
//graph G;
int cnt=0;
int cnt1;
//random_simple_undirected_graph(G,n,m);
      //random_planar_graph(G,n);
//edge_array<double> cost(G);
//node_array<double> dist(G);

edge e;

//edge_array<bool> b(G);
for(n=100;n<=1000;n=n+100)
{
graph G;
random_planar_graph(G,n);
edge_array<double> cost(G);
node_array<double> dist(G);
edge_array<bool> b(G);
forall_edges(e,G){
	cost[e]=((double) rand_int(0,100));
	b[e]=0;
	//cout<<b[e]<<"\n";
	}

//float T=used_time();
DIJKSTRA(G,G.first_node(),cost,dist,b);
float T=used_time();
cnt1=DIJKSTRA1(G,G.first_node(),cost,dist,b);

cout<<"\n\nThe shortest path computation took" << used_time(T)<<"seconds.\n\n";
cout<<"VVC:"<<cnt1;
forall_edges(e,G)
{
//cout<<b[e]<<"\n";
if(b[e]==1)
cnt++;
}
//cout<<cnt;
}
return 0;
}

