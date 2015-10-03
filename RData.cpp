//General DS for Shortest Path computation
#include <LEDA/graph/graph.h>
#include <LEDA/graph/node_pq.h>
#include<stdlib.h>
#include<math.h>
#include <fstream>
using std::ifstream;
#include <cstdlib> 
//#include<fstream.h>
#include<string.h>
//Used for co-ordinate generation
#include <LEDA/geo/point.h>
#include <LEDA/geo/random_point.h>

//LEDA Namespace is used generally
using namespace leda;
#if defined(LEDA_STD_IO_HEADERS)
using std::cout;
using std::endl;
#endif

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


class Flags {
  public:
    //bool region[GRIDCOUNT];
	Flags(){}
};

void reachvalue(const graph &G,const edge_array<double>& cost,node& s,node& d)
{
 
  node_pq<double> PQ(G);
  node v;
  edge e;

  node_array<double> dist(G);
  node_array<edge> na(G);    
  int nodecount = 0;
  
  forall_nodes(d,G)
      { 
          if (d == s) dist[d] = 0; else dist[d] = MAXDOUBLE;
          PQ.insert(d,dist[d]);
      }
  int vertexVisitCnt = 0;
  int a[100];
  int i=0;
  while ( !PQ.empty() ) 
  {
    node u = PQ.del_min();
	//Oops: All the vertices are not reachable from given source
    if( dist[u] == MAXDOUBLE ) 
    {
        PQ.clear();
        break;
    }

    vertexVisitCnt++;
    //Dijk Mod: If target is reached break
    if( u == d)
    {
        break;
    }
   int js = 0,a[200],b[200];
   FILE *f1;
   f1=fopen("ReachBased.txt","w");
  // cout<<"hi";
   fprintf(f1,"Reach of nodes\n");
  
    forall_out_edges(e,u) 
    {
        v = target(e);
        //cout <<"Reach of 10 nodes are :";
         
        for(int i=0; i< 2000; i++)
        {
        int c = dist[u] + cost[e];
        if ( c < dist[v] ) 
            {
          PQ.decrease_p(v,c);
		  dist[v] = c;		 
          //cout <<c<<endl;
          a[i]=c;
          fprintf(f1," Reach of path %d:%d ",js,c);
          js++;
          //fprintf(f1,"%d ",c);  
            }
         }
    }
  int t;
  for(int k=0;k<200;k++)
  {
      b[k]=a[k];
  }
  for(int l=0;l<199;l++)
  {
      for(int r=1;r<200;r++)
      {
        if(b[l]<b[r])
        {
            t=b[l];
            b[l]=b[r];
            b[r]=t;
        }
      }
  }
  int gp=b[0];
  
  int l=0,o=1;
  FILE *fn;
  fn=fopen("VALUESOFREACH.txt","w");
  for(int m=0;m<200;m++)
  {
      if(b[0]>=a[m])
          fprintf(fn,"%d",l);
      else
          fprintf(fn,"%d",o);
  }
  }
}
int dijk(const graph &G,const edge_array<double>& cost,node& s, node& t)
{ 
  FILE *f3;  
  node_pq<double> PQ(G);
  node v; edge e;
  int a;
  node_array<double> dist(G);

  f3=fopen("VALUESOFREACH.txt","r");
  fscanf(f3,"%d",&a);
  forall_nodes(v,G)
   { 
       if (v == s) dist[v] = 0; else dist[v] = MAXDOUBLE;
  //     while(!feof(f3))
    //    {
      //    if (a==1)
          PQ.insert(v,dist[v]);
        //}
    }

//infile.close();
  int vertexVisitCnt = 0;

  while ( !PQ.empty() )
  { node u = PQ.del_min();
    //Oops: All the vertices are not reachable from given source
    if( dist[u] == MAXDOUBLE ) {
        PQ.clear();
        break;
    }

    vertexVisitCnt++;
    if( u == t ) {
        break;
    }

    forall_out_edges(e,u)
    {
        v = target(e);
        
        double c = dist[u] + cost[e];

        if ( c < dist[v] ) {
          PQ.decrease_p(v,c);  dist[v] = c;  
        }
    } 
  }

  return vertexVisitCnt;
}

int dijk_orig(const graph &G,const edge_array<double>& cost, node s, node t) {
  node_pq<double> PQ(G);
  node v; edge e;
  node_array<double> dist(G);

  forall_nodes(v,G)
 {
    if (v == s) dist[v] = 0; else dist[v] = MAXDOUBLE;
      PQ.insert(v,dist[v]);
  }

  int vertexVisitCnt = 0;

  while ( !PQ.empty() ) 
   {
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

  
  cout.precision(3);
  cout << left;
  cout << endl;
  cout << setw(10) << "NodeCnt"       << setw(10) << "EdgeCnt" 
      << setw(25) << "Reach Pre Time"  
      << setw(10) << "Reach Rt"        << setw(20) << "Reach VVC"   
       << setw(10) << "Dijk Rt"       << setw(20) << "Dijk VVC"   
       << endl;
        
FILE *fp;
  int num=65,val[66][66],u=0,v=0;
int dis;
      //edge_array<double> d(G);
fp=fopen("adjlist_9.txt","r");
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
int node_num=1,j=0,col_num=1;  
       char c,ss[10]=" ";
  node tg,st;  //for(int curNodeCnt = start_n; curNodeCnt <= end_n; curNodeCnt += step ) {
      //Generate a planar graph
      graph G;

     // random_graph(G, curNodeCnt, rand_int(5,10) * curNodeCnt , true, true, true);
	complete_graph(G,num);
	    node_array<int> number(G);
forall_nodes(nk,G){
number[nk]=node_num++;}
      //random_planar_graph(G, curNodeCnt);
      //Make_Connected( G );
      //G.make_directed();

 cout << setw(10) <<num;     
   edge_array<double> d(G);
     //Cost Assignment to edges
      int edgeCount = 0;
      edge e;  
      forall_edges(e,G) {
          d[e] =0;// ((double) rand_int(1,100));
          //edgeCount++;
      }
  int i=1;        
if(fp==NULL)
cout<<"Can't open file"<<endl;

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
edgeCount++;}
      cout << setw(10) << edgeCount;
        node s = G.first_node();
        node t = G.last_node();
      
      double runTime = 0, runTime_orig = 0;
      int vertexCnt = 0, vertexCnt_orig = 0;
 float T = used_time();
forall_nodes(nk,G){

  T = used_time();
   vertexCnt_orig += dijk_orig(G, d, nk, t);
          runTime_orig += used_time( T );
}
 //    
   //   node t = G.last_node();
      
      //Preprocessing
      float TT = used_time();
      reachvalue(G,d,s,t);  
      float preprocessTime = used_time(TT);


      cout << setw(25) << preprocessTime;

      //Shortest path computation
      //node s = G.first_node();
     // node t = G.last_node();
      
     // double runTime = 0, runTime_orig = 0;
      //int vertexCnt = 0, vertexCnt_orig = 0;
forall_nodes(nk,G){
     // for(int sampleRun = 1; sampleRun <= MAXSAMPLERUNS; sampleRun++ ) {
          T = used_time();
         vertexCnt += dijk(G, d, nk, t );
        runTime += used_time( T );

        //  vertexCnt_orig += dijk_orig(G, d, nk, t);
          //runTime_orig += used_time( T );

          

}
          //s = G.choose_node(); t = G.choose_node();
     // }

      /*vertexCnt /= MAXSAMPLERUNS;
      runTime   /= MAXSAMPLERUNS;

      vertexCnt_orig /= MAXSAMPLERUNS;
      runTime_orig   /= MAXSAMPLERUNS;*/

      cout << setw(10) << runTime 
           << setw(20) << vertexCnt/35
           << setw(10) << runTime_orig
           << setw(20) << vertexCnt_orig/35 << endl;
 // }//Runtime calculations over Containers 
     
  return 0;
}
