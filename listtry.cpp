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
int main()
{
node_list m;
graph G;
int n=100;
node e;
node m1=(node)100;
int i=10;
random_planar_graph(G,n);

forall_nodes(e,G){
m.append(e);
G.print_node(e);
}
//e=(node)i;
G.print_node(m1);
}
