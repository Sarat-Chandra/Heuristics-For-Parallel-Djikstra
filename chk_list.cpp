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

int main(){
graph G;
node n;
complete_graph(G,10);
node_list perm;
forall_nodes(n,G){
perm.append(n);
}

while(!perm.empty())
{
n=perm.pop();
G.print_node(n);
}
}


