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


int dijk_orig(const graph &G,const edge_array<double>& cost, node s, node t) {
  node_pq<double> PQ(G);
  node v; edge e;
  node_array<double> dist(G);

  forall_nodes(v,G) {
    if (v == s) dist[v] = 0; else dist[v] = MAXDOUBLE;
      PQ.insert(v,dist[v]);
  }

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
        
        double c = dist[u] + cost[e];

        if ( c < dist[v] ) {
          PQ.decrease_p(v,c);  dist[v] = c;  
        }
    } 
  }

  return vertexVisitCnt;
}


int dijk_multi(const graph &G,const edge_array<double>& cost, node s, node t, edge_array<bool>& multi) {
  node_pq<double> PQ(G);
  node v; edge e;
  node_array<double> dist(G);

  forall_nodes(v,G) {
    if (v == s) dist[v] = 0; else dist[v] = MAXDOUBLE;
      PQ.insert(v,dist[v]);
  }

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
        
        double c = dist[u] + cost[e];

        if ( c < dist[v] ) {
          PQ.decrease_p(v,c);  dist[v] = c; 
	multi[e]=1; 
        }
    } 
  }
return vertexVisitCnt;
}

//Class for storing point values
class CPoint {
  public:
    double x;
    double y;
};

class Flags {
  public:
    bool region[GRIDCOUNT];

    Flags() {
      for(int i = 0; i < GRIDCOUNT; i++)
          region[i] = false; //Set every region reachable from a given arc
    }

    bool markRegion( const CPoint &p) {
        //Co ordinate range -32768 to 32767
        double region_range = (MAXCOORD - MINCOORD) / 2;
        int region_y = (p.y - MINCOORD) / region_range;
        int region_x = (p.x - MINCOORD) / region_range;

        region[ region_y * 6 + region_x ] = true;

        return true;
    }
};

void assignFlags(const graph &G,const edge_array<double>& cost, 
           edge_array<Flags>& ea, node_array<CPoint>& ncoord)
{ node_pq<double> PQ(G);
  node v,s; edge e;

  node_array<double> dist(G);
  node_array<edge> na(G);    
  int nodecount = 0;

  forall_nodes(s,G) {

    forall_nodes(v,G) {
      if (v == s) dist[v] = 0; else dist[v] = MAXDOUBLE;
        PQ.insert(v,dist[v]);
    }

  while ( !PQ.empty() ) {
    node u = PQ.del_min();

    //Oops: All the vertices are not reachable from given source
    if( dist[u] == MAXDOUBLE ) {
        PQ.clear();
        break;
    }

    if(u != s) {
      ea[ na[u] ].markRegion( ncoord[u] );
    }          

    forall_out_edges(e,u) {
        v = target(e);
        double c = dist[u] + cost[e];

        if( c < dist[v] ) {
          PQ.decrease_p(v,c);  dist[v] = c;  
         
          if(u==s)
            na[v]=e;
          else
            na[v]=na[u];
        }
    } 
  } //Container for edges with one source complete

  } //loop for all the nodes in the graph
}


int dijk(const graph &G,const edge_array<double>& cost, 
           edge_array<Flags>& ea, node_array<CPoint>& ncoord, node& s, node& t)
{ node_pq<double> PQ(G);
  node v; edge e;
  node_array<double> dist(G);

  forall_nodes(v,G)
  { if (v == s) dist[v] = 0; else dist[v] = MAXDOUBLE;
      PQ.insert(v,dist[v]);
  }

  int vertexVisitCnt = 0;

  //Find the region to which the target node belongs
  double region_range = (MAXCOORD - MINCOORD) / 2;
  int region_y = (ncoord[t].y - MINCOORD) / region_range;
  int region_x = (ncoord[t].x - MINCOORD) / region_range;

  int target_region = region_y * 6 + region_x;

  while ( !PQ.empty() )
  { node u = PQ.del_min();
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

    forall_out_edges(e,u)
    {
        //Check if the edge leads to the target node region
        if( ea[ e ].region[ target_region] == false )
            continue;

        v = target(e);
        
        double c = dist[u] + cost[e];

        if ( c < dist[v] ) {
          PQ.decrease_p(v,c);  dist[v] = c;  
        }
    } 
  }

  return vertexVisitCnt;
}

	 
int main()
{
graph G;
int vertexCnt_multi=0,vertexCnt_orig=0;
int num;
cout<<"Enter number of nodes:"<<endl;
cin>>num;
for(num=100;num<=1000;num=num+100){
random_graph(G, num, rand_int(5,10) * num , true, true, true);
edge e,e1;
float T;
Make_Connected( G );
      G.make_directed();
      edge_array<double> cost(G);
	edge_array<bool> multi(G);
	node_array<double> node_num(G);
	node_array<int> level_num(G);
	node_array<CPoint> coord(G);
      
      //Recollect that a 6 x 6 rectangular grid partition is assumed
      edge_array<Flags> ea(G); 
	node nd;
      int edgeCount = 0;
      //edge e;  
      forall_edges(e,G) {
          cost[e] = ((double) rand_int(1,100));
	multi[e]=0;
          edgeCount++;
      }
	forall_nodes(nd,G) {
          coord[nd].x = (double) rand_int(-32768, 32767);
          coord[nd].y = (double) rand_int(-32768, 32767);
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
cout<<"Total nodes:"<<num<<endl;
cout<<"Total edges:"<<edgeCount<<endl;
/*forall_nodes(n,G)
{
G.print_node(n);
cout<<"\t"<<level_num[n]<<endl;
}*/
node s=G.first_node();
node t1=G.last_node();
float d=used_time();
vertexCnt_orig = dijk_orig(G, cost, s, t1);
cout<<"orig runtime:"<<used_time(d)<<endl;
cout<<"Orig VVC:"<<vertexCnt_orig<<endl;

int level=1;
T = used_time();
 //= G.first_node();
forall_nodes(s,G){
	vertexCnt_multi=0;
       node t = G.last_node();
vertexCnt_multi = dijk_multi(G, cost, s, t, multi);
//cout<<vertexCnt_orig<<endl;
}
int cnt=0;
forall_edges(e,G)
{
if(multi[e]!=1)
{cnt++; 
G.del_edge(e);}
}
//cout<<cnt<<endl;
assignFlags(G, cost, ea, coord);
cout<<"Preprocess Time:"<<used_time(T)<<endl;
int vertexCnt=0,vertexCnt1=0;
//shortest path starts here

node t=G.last_node();
forall_nodes(s,G){
vertexCnt += dijk(G, cost, ea, coord, s, t );}
s=G.first_node();
float r;
r=used_time();
vertexCnt1 = dijk(G, cost, ea, coord, s, t );
cout<<"MULTILEVEL RUNNING TIME:"<<used_time(r)<<endl;
cout<<"MULTILEVEL ARCFLAGS VERTEX COUNT:"<<vertexCnt/num<<endl;
cout<<"--------------------------------------"<<endl;
}
return 0;
}
