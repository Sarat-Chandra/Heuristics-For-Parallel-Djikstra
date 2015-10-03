#include <LEDA/graph/graph.h>
#include <LEDA/graph/node_pq.h>
#include <LEDA/graph/node_list.h>
#include <LEDA/graph/node_array.h>

#include <omp.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

#if defined(LEDA_STD_IO_HEADERS)
using std::cout;
using std::endl;
#endif

#include <csignal>

#include <iomanip>
using std::left;
using std::setw;

void sHandler(int err){
    exit(1);
}

int DIJKSTRA(const graph &G, node s, 
              const edge_array<double>& cost,
              node_array<double>& dist, node t)
{ node_pq<double> PQ(G);
  node v; edge e;
  int vertex_visit_count = 0;

  //Initialization
  forall_nodes(v,G)
  { if (v == s) dist[v] = 0; else dist[v] = MAXDOUBLE;
    PQ.insert(v,dist[v]);
  }

  // main loop 
  while ( !PQ.empty()  ) { 
    node u = PQ.del_min();

    //Target unreachable from current source
    if( dist[u] == MAXDOUBLE ) {
      PQ.clear();
      break;
    }

    vertex_visit_count++;

    if( u == t ) {
        break;
    }

    forall_out_edges(e,u)
      { v = target(e);
        double c = dist[u] + cost[e];
        if ( c < dist[v] ) 
        {  PQ.decrease_p(v,c);  dist[v] = c;  
        }
      }
  }

  return vertex_visit_count;
}

int DIJKSTRA_BD(const graph &G, node s, 
              const edge_array<double>& cost,
              node_array<double>& dist, node t)
{ node_pq<double> PQ(G);
  node_pq<double> rPQ(G);

  node_array<double> rdist(G);
  node_list perm;

  node v; edge e;

  int vertex_visit_count = 0;

  forall_nodes(v,G)
  { if (v == s) dist[v] = 0; else dist[v] = MAXDOUBLE;
    PQ.insert(v,dist[v]);
  }

  forall_nodes(v,G)
  { if (v == t) rdist[v] = 0; else rdist[v] = MAXDOUBLE;
    rPQ.insert(v,rdist[v]);
  }

  while ( !PQ.empty() && !rPQ.empty() ) { 
    node u = PQ.del_min();

    //Target not reachable froum current source
    if( dist[u] == MAXDOUBLE ) {
       PQ.clear();
       break;       
    }

    vertex_visit_count++;

    if( u == t || perm.member( u ) )
        break;

    perm.append( u );

    node ru = rPQ.del_min();

    //Source cannot be reached from target
    if( rdist[ru] == MAXDOUBLE ) {
      rPQ.clear();
      break;
    }

    //vertex_visit_count++;

    if( ru == s || perm.member( ru ) )
        break;

    perm.append( ru );

    forall_out_edges(e,u)
    { v = target(e);
      double c = dist[u] + cost[e];
      if ( c < dist[v] ) 
      {  PQ.decrease_p(v,c);  dist[v] = c;  }
    }

    forall_in_edges(e,ru)
    { v = source(e);
      double c = rdist[ru] + cost[e];
      if ( c < rdist[v] ) 
      {  rPQ.decrease_p(v,c);  rdist[v] = c;  }
    }
  } //main while looop

  return vertex_visit_count;
}

int DIJKSTRA_BD_P(const graph &G, node s, 
              const edge_array<double>& cost,
              node_array<double>& dist, node t)
{ node_pq<double> PQ(G);
  node_pq<double> rPQ(G);

  node_array<double> rdist(G);
  node_list perm;

  node v; edge e;
  int vertex_visit_count = 0;

  //Initialization
  #pragma omp parallel sections default(shared) private(v)
  {
      #pragma omp section
      {
        forall_nodes(v,G)
        { if (v == s) dist[v] = 0; else dist[v] = MAXDOUBLE;
          PQ.insert(v,dist[v]);
        }
      }

      #pragma omp section
      {
        forall_nodes(v,G)
        { if (v == t) rdist[v] = 0; else rdist[v] = MAXDOUBLE;
          rPQ.insert(v,rdist[v]);
        }
      }
  } // Parallel region ends

  // main loop 
  while ( !PQ.empty() && !rPQ.empty() ) { 
    node u = PQ.del_min();

    //Target not reachable froum current source
    if( dist[u] == MAXDOUBLE ) {
       PQ.clear();
       break;       
    }

    vertex_visit_count++;

    if( u == t || perm.member( u ) )
        break;

    perm.append( u );

   node ru = rPQ.del_min();

    //Source cannot be reached from target
    if( rdist[ru] == MAXDOUBLE ) {
      rPQ.clear();
      break;
    }

    //vertex_visit_count++;

    if( ru == s || perm.member( ru ) )
        break;

    perm.append( ru );

    #pragma omp parallel sections default(shared) private(e, v)
    {
        #pragma omp section
        {
            forall_out_edges(e,u)
            { v = target(e);
              double c = dist[u] + cost[e];
              if ( c < dist[v] ) 
              {  PQ.decrease_p(v,c);  dist[v] = c;  }
            }
        }

        #pragma omp section
        {
            forall_in_edges(e,ru)
            { v = source(e);
              double c = rdist[ru] + cost[e];
              if ( c < rdist[v] ) 
              {  rPQ.decrease_p(v,c);  rdist[v] = c;  }
            }
        }
    } //Parallel region ends
  } //main while looop

  return vertex_visit_count;
}

int main()
{int count=0;
//  int start_n = read_int("Number of nodes (Starting Count) : ");
  //int end_n   = read_int("Number of nodes (Ending Count)   : ");
  //int step    = read_int("Step Value for Node Increment    : ");
  FILE *fp;
  int val[18][18],u=0,v=0;

  signal( SIGSEGV,sHandler );
  signal( SIGABRT,sHandler );

  cout << left;
  cout.precision(3);
  cout << setw(10) << "NodeCnt"          << setw(10) << "EdgeCnt" 
       << setw(15) << "Dijk Rt"          << setw(15) << "BiDir Rt"
       << setw(15) << "BiDirP Rt"        << setw(10) << "DijkVCnt"
       << setw(15) << "BiDirVCnt"        << endl;

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
node nk;//for(int curNodeCnt = start_n; curNodeCnt <= end_n; curNodeCnt += step ) {
    
    int edgeCnt = 0;
    double t_dijk, t_bd_dijk, t_bd_dijk_p;
    int v_dijk, v_bd_dijk, v_bd_dijk_p;

    t_dijk = t_bd_dijk = t_bd_dijk_p = 0;
    v_dijk = v_bd_dijk = v_bd_dijk_p = 0;

    cout << setw(10) << "17";

   for( int i = 1; i <= 20; i++ ) {
      graph G;
      int node_num=1,j=0,col_num=1;  
       char c,ss[10]=" ";
  node tg,st;

//      random_graph(G,curNodeCnt, curNodeCnt * 10, true, true, true);
	complete_graph(G,17);
    node_array<int> number(G);
//node_array<int> number(G);
     // G.make_directed();
      //node_array<int> number(G);
      int dis;
      edge_array<double> d(G);
      node_array<double> dist(G);
      edge e;

      forall_edges(e,G) { 
          d[e] = 0;//((double) rand_int(0,100));
          edgeCnt++;
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

      node t = G.last_node();
forall_nodes(nk,G)
{
      float T = used_time();
      v_dijk += DIJKSTRA(G,nk,d,dist,t);
      t_dijk += used_time(T);

      v_bd_dijk += DIJKSTRA_BD(G,nk,d,dist, t);
      t_bd_dijk += used_time(T);

      v_bd_dijk_p += DIJKSTRA_BD_P(G, nk, d, dist, t);
      t_bd_dijk_p += used_time(T);}

forall_edges(e,G){
count++;
}
  }

   cout << setw(10) << count     	 << setw(15) << t_dijk
        << setw(15) << t_bd_dijk         << setw(15) << t_bd_dijk_p
        << setw(10) << v_dijk/20            << setw(15) << v_bd_dijk/20 << endl;

 // }//Over with the range of vertices

cout<<count;
fclose(fp);  
  return 0;
}
