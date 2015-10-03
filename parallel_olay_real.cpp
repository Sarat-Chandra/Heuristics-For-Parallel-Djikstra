#include<LEDA/graph/graph.h>
#include<LEDA/graph/node_pq.h>
#include<LEDA/graph/node_list.h>
#include<stdio.h>
#include<omp.h>
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
#pragma omp parallel sections default(shared) private(v)
  {
      #pragma omp section
      {
  forall_nodes(v,G) {
    if (v == s) dist[v] = 0; else dist[v] = MAXDOUBLE;
      PQ.insert(v,dist[v]);
  }
}}
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


int main()
{
graph G;
int vertexCnt_orig=0,vertexCnt_multi=0,vertexCnt_multi1=0;
//int num=35;
int num=100, val[101][101],u,v;
char c;
//cout<<"Enter number of nodes:"<<endl;
//cin>>num;
//for(num=100;num<=1000;num=num+100){
//random_planar_graph(G, num);
complete_graph(G,num);
edge e,e1;
float T;
//Make_Connected( G );
  //    G.make_directed();
      edge_array<double> d(G);
	edge_array<bool> multi(G);
	node_array<double> node_num(G);
	node_array<int> level_num(G);
      int edgeCount = 0;
	node_array<int> number(G);

      
FILE *fp,*fp_cities;
char cities[100];
 //FILE *fp_p= fopen("path.txt","w");
fp_cities=fopen("cities.txt","r");
while(!feof(fp_cities))
{
fgets(cities,100,fp_cities);
cout<<cities<<endl;}
fclose(fp_cities);
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


      //edge e;  
      forall_edges(e,G) {
          d[e] = 99999;//((double) rand_int(1,100));
	multi[e]=0;
          //edgeCount++;
      }


int i=0,j=0,col_num=1,dist=0,node_num1=1,src,tgt;
char ss[10]=" ";
node st,tg;
node nd;

forall_nodes(nd,G){
number[nd]=node_num1++;}

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
	{dist=atoi(ss);
	col_num++;}//if3
	else
	{//if3
	if((col_num%2)==1){//if3
	v=atoi(ss);
	col_num++;
	val[u][v]=dist;
	forall_edges(e,G){
	st=G.source(e);
	tg=G.target(e);
	if((number[st]==u) && (number[tg]==v))
	{
	d[e]=dist;}}
	}//if3
	}//if3
//	cout<<ss<<endl;
	j++; 
	i=0;}//else2 ends
}//if1 ends
else
{col_num=1;}
}

int  edgCnt=0;
forall_edges(e,G){
if(d[e]==9999)
G.del_edge(e);
else
edgCnt++;
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
float d1=used_time();
node s=G.first_node();
node t=G.last_node();
vertexCnt_orig=dijk_orig(G, d, s, t);
cout<<"ORIG DIJK RT:"<<used_time(d1)<<endl;
cout<<"ORIG DIJK VVC:"<<vertexCnt_orig<<endl;
/*forall_nodes(n,G)
{
G.print_node(n);
cout<<"\t"<<level_num[n]<<endl;
}*/
int level=1;
T = used_time();
//node s; //= G.first_node();
forall_nodes(s,G){
	vertexCnt_multi=0;
      node t = G.last_node();
vertexCnt_multi = dijk_multi(G, d, s, t, multi);
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
cout<<"Multilevel Preprocess Time:"<<used_time(T)<<endl;
s=G.first_node();
t=G.last_node();
vertexCnt_multi1=dijk_orig(G, d, s, t);
cout<<"MULTILEVEL RT:"<<used_time(T)<<endl;
cout<<"MULTILEVEL VVC:"<<vertexCnt_multi1<<endl<<endl<<endl;

int cnte=0;
forall_edges(e,G){
cnte++;}
cout<<cnte;
//}
return 0;
}
