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
int main()
{
graph G;
void *p;
edge n;
int cnt=0;
node e,e1;
random_planar_graph(G, 10);
e=G.first_node();
e1=G.last_node();
G.add_edge(e,e1,p);
forall_edges(n,G)
{
cnt++;
e=source(n);
e1=target(n);
G.print_edge(n);
G.del_edge(n);};
cout<<cnt<<endl;
forall_edges(n,G)
{
cout<<"edge"<<endl;}
}
