                                                                     
                                                                     
                                                                     
                                             
//General DS for Shortest Path computation
#include <LEDA/graph/graph.h>
#include <LEDA/graph/node_pq.h>

//Used for co-ordinate generation
#include <LEDA/geo/point.h>
#include <LEDA/geo/random_point.h>
#include<omp.h>
//LEDA Namespace is used generally
using namespace leda;

//C++ Headers n Libs
#include <iostream>
using std::cout;
using std::endl;

#include <iomanip>
using std::left;
using std::setw;

//Custom Values
#define MAX_DOUBLE   21474
#define MIN_DOUBLE  -21474
#define MAXSAMPLERUNS 10

//Assuming a 6 x 6 rectangular grid
#define GRIDCOUNT 36
#define MAXCOORD 32767
#define MINCOORD -32768

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
#pragma omp parallel sections default(shared) private(v)
  {
      #pragma omp section
      {
  forall_nodes(s,G) {
#pragma omp section
      {
    forall_nodes(v,G) {
      if (v == s) dist[v] = 0; else dist[v] = MAXDOUBLE;
        PQ.insert(v,dist[v]);
    }
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
}}
  } //loop for all the nodes in the graph
}
	 
int dijk(const graph &G,const edge_array<double>& cost, 
           edge_array<Flags>& ea, node_array<CPoint>& ncoord, node& s, node& t, node_array<double>& dist)
{ node_pq<double> PQ(G);
  node_pq<double> rPQ(G);
  node v; edge e;
  node_array<double> rdist(G);
  node_list perm;

  forall_nodes(v,G)
  { if (v == s) dist[v] = 0; else dist[v] = MAXDOUBLE;
      PQ.insert(v,dist[v]);
  }

   forall_nodes(v,G)
	{
	if(v==t)
	rdist[v]=0;
	else
	rdist[v]=MAXDOUBLE;
	rPQ.insert(v,rdist[v]);
	}

  int vertexVisitCnt = 0;

  //Find the region to which the target node belongs
  double region_range = (MAXCOORD - MINCOORD) / 2;
  int region_y = (ncoord[t].y - MINCOORD) / region_range;
  int region_x = (ncoord[t].x - MINCOORD) / region_range;

  int target_region = region_y * 6 + region_x;

  while ( !PQ.empty()&& !rPQ.empty() )
  { node u = PQ.del_min();
    //Oops: All the vertices are not reachable from given source
    if( dist[u] == MAXDOUBLE ) {
        PQ.clear();
        break;
    }

    vertexVisitCnt++;
    //Dijk Mod: If target is reached break
    if( u == t || perm.member(u)) {
        break;
    }

   perm.append(u);
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
node ru=rPQ.del_min();
	if(ru==s || perm.member(ru)){
	//cout<<"\n\nTarget reached in reverse search"<<endl;
	//cout<<"Common node:";
	//G.print_node(ru);
	break;
	}

	//cout<<"\tBackward Search V:";
	//G.print_node(ru);
	perm.append(ru);

	forall_out_edges(e,ru)
	{
	v=target(e);
	double c=rdist[ru]+cost[e];
	if(c<rdist[v]){
	rPQ.decrease_p(v,c);
	rdist[v]=c;}
	}

  }

  return vertexVisitCnt;
}

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

int main() {
//  int start_n = read_int("Number of nodes (Starting Count) : ");
//  int end_n   = read_int("Number of nodes (Ending Count)   : ");
//  int step    = read_int("Step Value for Node Increment    : ");

  cout.precision(3);
  cout << left;
  cout << endl;
  cout << setw(10) << "NodeCnt"       << setw(10) << "EdgeCnt" 
       << setw(25) << "Arc Pre Time" 
       << setw(10) << "Arc Rt"        << setw(20) << "Arc VVC"   
       << setw(10) << "Dijk Rt"       << setw(20) << "Dijk VVC"   
       << endl;
  FILE *fp;
  int num=26,val[27][27],u=0,v=0;
fp=fopen("adj_list.txt","r");
for(u=1;u<=num;u++){
for(v=1;v<=num;v++){
val[u][v]=0;
//cout<<val[u][v]<<"\t";
}
//cout<<endl;
}
u=1;
v=1;
node nk;      
 // for(int curNodeCnt = start_n; curNodeCnt <= end_n; curNodeCnt += step ) {
      //Generate a planar graph
      graph G;

      //random_graph(G, curNodeCnt, rand_int(5,10) * curNodeCnt , true, true, true);
   //   random_planar_graph(G, curNodeCnt);
   complete_graph(G,num);
	    node_array<int> number(G);
	int dis;
   //   Make_Connected( G );
    //  G.make_directed();
      int node_num=1,j=0,col_num=1,i=1;  
       char c,ss[10]=" ";
  	node tg,st;
      cout << setw(10) << num; 

      edge_array<double> d(G);
      node_array<CPoint> coord(G);
      node_array<double> dist(G);     
      //Recollect that a 6 x 6 rectangular grid partition is assumed
      edge_array<Flags> ea(G);
      

      //Cost Assignment to edges
      int edgeCount = 0;
      edge e;  
       forall_edges(e,G) { 
          d[e] = 0;//((double) rand_int(0,100));
          //edgeCount++;
      }
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
	{dis=atoi(ss);
	col_num++;}//if3
	else
	{//if3
	if((col_num%2)==1){//if3
	v=atoi(ss);
	col_num++;
	val[u][v]=dis;
	forall_edges(e,G){
	st=G.source(e);
	tg=G.target(e);
	if((number[st]==u) && (number[tg]==v))
	{
	d[e]=dis;}}
	}//if3
	}//if3
//	cout<<ss<<endl;
	j++; 
	i=0;}//else2 ends
}//if1 ends
else
{col_num=1;}
}
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
     forall_edges(e,G){
edgeCount++;
}
 cout << setw(10) << edgeCount;

      //Coordinate generation (Layout)
      node nd; 
//      int i;
      forall_nodes(nd,G) {
          coord[nd].x = (double) rand_int(-32768, 32767);
          coord[nd].y = (double) rand_int(-32768, 32767);
      }    

      //Preprocessing
      float T = used_time();
      assignFlags(G, d, ea, coord);
      float preprocessTime = used_time(T);

      cout << setw(25) << preprocessTime;
 double runTime = 0, runTime_orig = 0;
      int vertexCnt = 0, vertexCnt_orig = 0;
node s;
     forall_nodes(s,G){ //Shortest path computation
     // node s = G.first_node();
      node t = G.last_node();
      
        //  for(int sampleRun = 1; sampleRun <= MAXSAMPLERUNS; sampleRun++ ) {
          T = used_time();
          vertexCnt += dijk(G, d, ea,coord, s, t, dist );
          runTime += used_time( T );

          vertexCnt_orig += dijk_orig(G, d, s, t);
          runTime_orig += used_time( T );

         // s = G.choose_node(); t = G.choose_node();
     // }
}
     /* vertexCnt /= MAXSAMPLERUNS;
      runTime   /= MAXSAMPLERUNS;

      vertexCnt_orig /= MAXSAMPLERUNS;
      runTime_orig   /= MAXSAMPLERUNS;
*/
      cout << setw(10) << runTime
           << setw(20) << vertexCnt/35
           << setw(10) << runTime_orig
           << setw(20) << vertexCnt_orig/35 << endl;
//  }//Runtime calculations over Containers 
     
  return 0;
}


