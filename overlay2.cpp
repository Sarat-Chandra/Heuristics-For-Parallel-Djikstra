#include<LEDA/graph/graph.h>
#include<LEDA/graph/node_pq.h>
#include<LEDA/graph/node_list.h>
#include<stdio.h>
#if defined(LEDA_NAMESPACE)
using namespace leda;
using namespace std;
#endif
#if defined(LEDA_STD_IO_HEADERS)
//using std::cout;
//using std::endl;
#endif

int dijk_orig(const graph &G,const edge_array<double>& cost, node s, node t, edge_array<bool>& multi, int level, node_array<int>& level_num) {
  node_pq<double> PQ(G);
  node v; edge e;
  node_array<double> dist(G);
  int v_cnt=0;

  forall_nodes(v,G) {
if(level_num[v]<=level){
	v_cnt++;
    if (v == s) dist[v] = 0; else dist[v] = MAXDOUBLE;
      PQ.insert(v,dist[v]);}
  }
cout<<"Number of nodes inserted:"<<v_cnt<<endl;
  int vertexVisitCnt = 0;

  while ( !PQ.empty() ) {
    node u = PQ.del_min();
    //Oops: All the vertices are not reachable from given source
    if( dist[u] == MAXDOUBLE ) {
        PQ.clear();
        break;
    }

    vertexVisitCnt++;
    //Dijk Mod: If target is reached break
    if( u == t ) {
        break;
    }

    forall_out_edges(e,u) {
        v = target(e);
     //if(level_num[v]>level)
	//continue;   
        double c = dist[u] + cost[e];

        if ( c < dist[v] ) {
          PQ.decrease_p(v,c);  dist[v] = c; 
	multi[e]=1; 
        }
    } 
  }
return vertexVisitCnt;
}


int main()
{
graph G;
int vertexCnt_orig=0;
int num;
cout<<"Enter number of nodes:"<<endl;
cin>>num;
random_graph(G, num, rand_int(5,10) * 100 , true, true, true);
edge e,e1;
Make_Connected( G );
      G.make_directed();
      edge_array<double> cost(G);
	edge_array<bool> multi(G);
	node_array<double> node_num(G);
	node_array<int> level_num(G);
      int edgeCount = 0;
      //edge e;  
      forall_edges(e,G) {
          cost[e] = ((double) rand_int(1,100));
	multi[e]=0;
          edgeCount++;
      }

	node n;
	int numb=1;
	forall_nodes(n,G){
	if(numb>num/3 && numb<=2*num/3)
	level_num[n]=2;
	else if(numb<=num/3)
	level_num[n]=1;
	else
	level_num[n]=3;
	node_num[n]=numb++;}	
cout<<"Total edges:"<<edgeCount<<endl;
forall_nodes(n,G)
{
G.print_node(n);
cout<<"\t"<<level_num[n]<<endl;
}
int level=1;
int vcount=0;
int ecnt=0;
node s;
node t;
for(level=1;level<=3;level++){
vcount=0;
ecnt=0;
cout<<"LEVEL:"<<level<<endl;
forall_nodes(n,G){
if(level_num[n]==level)
s=n;
}
forall_nodes(n,G){
t=n;
if(level_num[n]>level)
break;
}
vcount=dijk_orig(G,cost,s,t,multi,level,level_num);
//forall_nodes(n,G)
//{
//if(level_num[n]<=level)
//{
forall_edges(e,G){
if(multi[e]==1)
ecnt++;
//}
//}
}
cout<<"Total edges in overlay graph:"<<ecnt<<endl;
}

//node s; //= G.first_node();
forall_nodes(s,G){
	vertexCnt_orig=0;
      node t = G.last_node();
//vertexCnt_orig = dijk_orig(G, cost, s, t, multi);
//cout<<vertexCnt_orig<<endl;
}
int cnt=0;
forall_edges(e,G)
{
if(multi[e]!=1)
{cnt++; 
G.del_edge(e);}
}
cout<<cnt<<endl;
return 0;
}
