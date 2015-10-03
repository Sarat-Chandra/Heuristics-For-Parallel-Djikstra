//General DS for Shortest Path computation
#include <LEDA/graph/graph.h>
#include <LEDA/graph/node_pq.h>

//Used for co-ordinate generation
#include <LEDA/geo/point.h>
#include <LEDA/geo/random_point.h>

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
#define MAXSAMPLERUNS 10

//Assuming a 6 x 6 rectangular grid
#define GRIDCOUNT 36
#define MAXCOORD 32767
#define MINCOORD -32768

#define LANDMARKCOUNT 3

/*************************************************************************************/

class node_LandmarkInfo {
  public:
    node landmark;
    double dist;
};

class LandmarkDistance{
  public:
    node_LandmarkInfo nLInfo[LANDMARKCOUNT];
};

void landmark_selection(const graph &G,const edge_array<double>& cost, 
                        node_list &landmarks) {
  node_pq<double> PQ(G);
  node v,s; edge e;

  node_array<double> dist(G);
  node_array<edge> na(G);    
  node u;

  s = G.choose_node(); //Start from s for choosing landmarks

  for(int landmarkSelected = 0; landmarkSelected < LANDMARKCOUNT; landmarkSelected++ ) {
    forall_nodes(v,G) {
        if(v == s) dist[v] = 0; else dist[v] = MAXDOUBLE;
        PQ.insert(v,dist[v]);
    }

    while ( !PQ.empty() ) {
        u = PQ.del_min();
        
        if( dist[u] == MAXDOUBLE ) {
            PQ.clear();
            break;
        }


        forall_adj_edges(e,u) {
           v = target(e);

           double c = dist[u] + cost[e];

           if ( c < dist[v] ) {
            PQ.decrease_p(v,c);  dist[v] = c;  
           }
        }

    }

    landmarks.append( u ); //Select the farthest node from s as landmark
    s = u; //Next Source for another landmark selection
  } //Landmark Selection loop
}

void landmark_distanceCal(const graph &G,const edge_array<double>& cost, 
                           const node_list &landmarks, node_array<LandmarkDistance> &nodeLandmarkInfo) {
  node_pq<double> PQ(G);
  node v,s; edge e;
  int landmarkCount = 0;

  forall(s, landmarks)
  {
    forall_nodes(v,G) {
        nodeLandmarkInfo[v].nLInfo[landmarkCount].landmark = s;

        if (v == s) 
            nodeLandmarkInfo[v].nLInfo[landmarkCount].dist = 0; 
        else 
            nodeLandmarkInfo[v].nLInfo[landmarkCount].dist = MAXDOUBLE;

        PQ.insert(v, nodeLandmarkInfo[v].nLInfo[landmarkCount].dist);
    }

    while ( !PQ.empty() )
    { node u = PQ.del_min();

      if( nodeLandmarkInfo[u].nLInfo[landmarkCount].dist == MAXDOUBLE ) {
          PQ.clear();
          break; //vertex not reachable
      }

        forall_adj_edges(e,u) {
           v = target(e);
           double c = nodeLandmarkInfo[u].nLInfo[landmarkCount].dist + cost[e];

           if ( c < nodeLandmarkInfo[v].nLInfo[landmarkCount].dist ) {
               PQ.decrease_p(v,c);  
               nodeLandmarkInfo[v].nLInfo[landmarkCount].dist = c;  
           }
        }
    }
    
    landmarkCount++;
  }
}

/*************************************************************************************/

//Class for storing point values
class CPoint {
  public:
    double x;
    double y;
};

/*class Flags {
  public:
    bool region[GRIDCOUNT];

    Flags() {
      for(int i = 0; i < GRIDCOUNT; i++)
          region[i] = false;
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
}*/
	 
int dijk(const graph &G,const edge_array<double>& cost, 
         const node_list &landmarks, const node_array<LandmarkDistance> &nodeLandmarkInfo, 
         node& s, node& t)
{ node_pq<double> PQ(G);
  node v; edge e;
  node_array<double> dist(G);

  forall_nodes(v,G)
  { if (v == s) dist[v] = 0; else dist[v] = MAXDOUBLE;
      PQ.insert(v,dist[v]);
  }

  int vertexVisitCnt = 0;
  
  double maxdiff, diff;
  
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
        /*if( ea[ e ].region[ target_region] == false )
            continue;*/

        v = target(e);

        //Include Landmark based potentials also
        maxdiff = nodeLandmarkInfo[v].nLInfo[0].dist - nodeLandmarkInfo[t].nLInfo[0].dist;

        for(int landmarkCount = 1; landmarkCount < LANDMARKCOUNT; landmarkCount++) {
            diff = nodeLandmarkInfo[v].nLInfo[landmarkCount].dist - 
                   nodeLandmarkInfo[t].nLInfo[landmarkCount].dist; //Triangle Inequality part

            if( diff > maxdiff )
                maxdiff = diff; //Choose the max Lower Bound
        }

        if(maxdiff < 0) maxdiff = -maxdiff; //Needs Correction: Update the potential logic

        double c = dist[u] + cost[e] + maxdiff;
        
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
  int start_n = read_int("Number of nodes (Starting Count) : ");
  int end_n   = read_int("Number of nodes (Ending Count)   : ");
  int step    = read_int("Step Value for Node Increment    : ");

  cout.precision(3);
  cout << left;
  cout << endl;
  cout << setw(10) << "NodeCnt"       << setw(10) << "EdgeCnt" 
       << setw(13) /*<< "Arc PT"        << setw(12)*/ << "LandM PT" 
       << setw(10) << "LM Rt"        << setw(20) << "LM VVC"   
       << setw(10) << "Dijk Rt"       << setw(20) << "Dijk VVC"   
       << endl;
        
  for(int curNodeCnt = start_n; curNodeCnt <= end_n; curNodeCnt += step ) {
      //Generate a planar graph
      graph G;

      //Choose the graph type required, could comment out the other one

      //random_graph(G, curNodeCnt, rand_int(5,10) * curNodeCnt , true, true, true);
      random_planar_graph(G, curNodeCnt);
      Make_Connected( G );
      G.make_directed();

      cout << setw(10) << curNodeCnt; 

      edge_array<double> cost(G);
      node_array<CPoint> coord(G);
      
      //Here a 6 x 6 rectangular grid partition is assumed
      //edge_array<Flags> ea(G); 

      //Cost Assignment to edges
      int edgeCount = 0;
      edge e;  
      forall_edges(e,G) {
          cost[e] = ((double) rand_int(1,100));
          edgeCount++;
      }

      cout << setw(10) << edgeCount;

      //Coordinate generation
      node nd; 
      int i;
     /* forall_nodes(nd,G) {
          coord[nd].x = (double) rand_int(-32768, 32767);
          coord[nd].y = (double) rand_int(-32768, 32767);
      } */   

      //Preprocessing
      float T = used_time();
      //assignFlags(G, cost, ea, coord);
      //float preprocessTimeA = used_time(T);

      node_list landmarks;
      node_array<LandmarkDistance> nodeLandmarkInfo(G);

      //Landmark selection and distance calculation
      T = used_time();
      landmark_selection(G, cost, landmarks); 

      landmark_distanceCal(G, cost, landmarks, nodeLandmarkInfo );
      float preprocessTimeL = used_time( T );

      //cout << setw(13) << preprocessTimeL;
      cout << setw(12) << preprocessTimeL;

      //Shortest path computation
      node s = G.first_node();
      node t = G.last_node();
      
      double runTime = 0, runTime_orig = 0;
      int vertexCnt = 0, vertexCnt_orig = 0;

      for(int sampleRun = 1; sampleRun <= MAXSAMPLERUNS; sampleRun++ ) {
          T = used_time();
          vertexCnt += dijk(G, cost, landmarks, nodeLandmarkInfo, s, t );
          runTime += used_time( T );

          vertexCnt_orig += dijk_orig(G, cost, s, t);
          runTime_orig += used_time( T );

          s = G.choose_node(); t = G.choose_node();
      }

      vertexCnt /= MAXSAMPLERUNS;
      runTime   /= MAXSAMPLERUNS;

      vertexCnt_orig /= MAXSAMPLERUNS;
      runTime_orig   /= MAXSAMPLERUNS;

      cout << setw(10) << runTime 
           << setw(20) << vertexCnt 
           << setw(10) << runTime_orig
           << setw(20) << vertexCnt_orig << endl;
  }//Runtime calculations over Containers 
     
  return 0;
}
