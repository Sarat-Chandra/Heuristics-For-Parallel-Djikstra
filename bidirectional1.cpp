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

    

    if( u == t || perm.member( u ) )
        break;

    perm.append( u );

	vertex_visit_count++;    
	node ru = rPQ.del_min();

    //Source cannot be reached from target
    if( rdist[ru] == MAXDOUBLE ) {
      rPQ.clear();
      break;
    }


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

    //vertex_visit_count++;

    if( u == t || perm.member( u ) )
        break;

    perm.append( u );

   node ru = rPQ.del_min();

    //Source cannot be reached from target
    if( rdist[ru] == MAXDOUBLE ) {
      rPQ.clear();
      break;
    }

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
{
  int start_n = read_int("Number of nodes (Starting Count) : ");
  int end_n   = read_int("Number of nodes (Ending Count)   : ");
  int step    = read_int("Step Value for Node Increment    : ");

  signal( SIGSEGV,sHandler );
  signal( SIGABRT,sHandler );

  cout << left;
  cout.precision(3);
  cout << setw(10) << "NodeCnt"          << setw(10) << "EdgeCnt" 
       << setw(15) << "Dijk Rt"          << setw(15) << "BiDir Rt"
       << setw(15) << "BiDirP Rt"        << setw(10) << "DijkVCnt"
       << setw(15) << "BiDirVCnt"        << endl;

  for(int curNodeCnt = start_n; curNodeCnt <= end_n; curNodeCnt += step ) {

    int edgeCnt = 0;
    double t_dijk, t_bd_dijk, t_bd_dijk_p;
    int v_dijk, v_bd_dijk, v_bd_dijk_p;

    t_dijk = t_bd_dijk = t_bd_dijk_p = 0;
    v_dijk = v_bd_dijk = v_bd_dijk_p = 0;

    cout << setw(10) << curNodeCnt;

    for( int i = 1; i <= 20; i++ ) {
      graph G;
      random_graph(G,curNodeCnt, curNodeCnt * 10, true, true, true);
      G.make_directed();

      edge_array<double> cost(G);
      node_array<double> dist(G);
      edge e;

      forall_edges(e,G) { 
          cost[e] = ((double) rand_int(0,100));
          edgeCnt++;
      }

      node t = G.last_node();

      float T = used_time();
      v_dijk += DIJKSTRA(G,G.first_node(),cost,dist, t);
      t_dijk += used_time(T);

      v_bd_dijk += DIJKSTRA_BD(G,G.first_node(),cost,dist, t);
      t_bd_dijk += used_time(T);

      //v_bd_dijk_p += DIJKSTRA_BD_P(G, G.first_node(), cost, dist, t);
      //t_bd_dijk_p += used_time(T);
   }

   cout << setw(10) << edgeCnt / 20      << setw(15) << t_dijk
        << setw(15) << t_bd_dijk         << setw(15) << t_bd_dijk_p
        << setw(10) << v_dijk            << setw(15) << v_bd_dijk << endl;

  }//Over with the range of vertices
  
  return 0;
}
