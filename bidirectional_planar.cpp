#include<LEDA/graph/graph.h>
#include<LEDA/graph/node_pq.h>
#include<LEDA/graph/node_list.h>
#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif
#if defined(LEDA_STD_IO_HEADERS)
using std::cout;
using std::endl;
#endif

void DIJKSTRA(const graph &G, node s, const edge_array<double> &cost, node_array<double> &dist)
{
node_pq<double> PQ(G);
node_pq<double> rPQ(G);
node_array<double> rdist(G);
node_list perm;

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

forall_nodes(v,G)
{
if(v==t)
rdist[v]=0;
else rdist[v]=MAXDOUBLE;
rPQ.insert(v,rdist[v]);
}


while(!PQ.empty() && !rPQ.empty()){
node u=PQ.del_min();
/*cout<<"\nF Search V:";
G.print_node(u);
perm.append(u);*/
if(u==t || perm.member(u))
{
cout<<"\n\nTarget reached in forward search"<<endl;
cout<<"common node:";
G.print_node(u);
break;
}

cout<<"\nF Search V:";
G.print_node(u);
perm.append(u);

forall_out_edges(e,u)
{
v=target(e);
double c=dist[u] + cost[e];
if(c<dist[v])
{PQ.decrease_p(v,c); dist[v]=c;}
}

node ru=rPQ.del_min();
if(ru==s||perm.member(ru))
{
cout<<"\n\nTarget reached in reverse search"<<endl;
cout<<"common node:"; G.print_node(ru);
break;
}

cout<<"\tB Search V:";
G.print_node(ru);
perm.append(ru);

forall_out_edges(e,ru)
{
v=target(e);
double c= rdist[ru]+cost[e];
if(c<rdist[v])
{
rPQ.decrease_p(v,c);
rdist[v]=c;
}
}

}

}

int main()
{
int n=read_int("number of nodes=");
int m=read_int("number of edges=");
graph G;
//random_simple_undirected_graph(G,n,m);
      random_planar_graph(G,n);
edge_array<double> cost(G);
node_array<double> dist(G);

edge e;

forall_edges(e,G)
cost[e]=((double) rand_int(0,100));

float T=used_time();
DIJKSTRA(G,G.first_node(),cost,dist);
cout<<"\n\nThe shortest path computation took" << used_time(T)<<"seconds.\n\n";
return 0;
}
