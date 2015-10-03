#include<iostream>
#include<stdio.h>
#include <cstdlib>
#define INFINITY 999
using namespace std;

FILE *fp = fopen("input.txt","w");
FILE *fp1 = fopen("overlay1.txt","w");
//FILE *fp2 = fopen("overlay2.txt","w");
int ov,ov1;
class Graph
{
  private:
    int adjMatrix[15][15];
    int predecessor[15],distance[15];
    bool mark[15];
    int source;
    long numOfVertices;
  public:
  	void read();
  	void initialize();
  	int getClosestUnmarkedNode();
  	void dijkstra();    
    void output();
    void printPath(int);    
};
void Graph::read()
{
  //cout<<"Enter the number of vertices of the graph(should be > 0)\n";
  int k;
  //cin>>numOfVertices;
  numOfVertices=rand() % 15;
  //numOfVertices=16;
  cout<<numOfVertices<<endl;
  for (int i=1;i<=numOfVertices; i++){
	for(int j=1;j<=numOfVertices; j++) {
		k=rand() % numOfVertices;		
		if(i==j)
		adjMatrix[i][j]=999;
		else if(j==k)
		adjMatrix[i][j]=999;
		else{
		adjMatrix[i][j]=rand() % 100;}
		cout<<adjMatrix[i][j]<<"\t";
		fprintf(fp,"%d",adjMatrix[i][j]);
		fprintf(fp,"\t");
	}
	cout<<endl;
	fprintf(fp,"\n");
  }
cout<<"Enter the source vertex\n";  
  cin>>source;
  while((source<0) && (source>numOfVertices-1)) {
    cout<<"Source vertex should be between 0 and "<<numOfVertices-1<<endl;
    cout<<"Enter the source vertex again\n";
    cin>>source;
  }
ov=numOfVertices/4;
}
void Graph::initialize()
{
  for(int i=0;i<numOfVertices;i++) {
    mark[i] = false;
    predecessor[i] = -1;
    distance[i] = INFINITY;
  }
  distance[source] = 0;
}

int Graph::getClosestUnmarkedNode()
{
  int minDistance = INFINITY;
  int closestUnmarkedNode;
  for(int i=0;i<ov1;i++) {
    if((!mark[i]) && ( minDistance >= distance[i])) {
		minDistance = distance[i];
		closestUnmarkedNode = i;
    }
  }
  return closestUnmarkedNode;
}
void Graph::dijkstra()
{
  initialize();
  int minDistance = INFINITY;
  int closestUnmarkedNode;
  int count = 0;
  while(count < ov1) {
    closestUnmarkedNode = getClosestUnmarkedNode();
    mark[closestUnmarkedNode] = true;
    for(int i=0;i<ov1;i++) {
      if((!mark[i]) && (adjMatrix[closestUnmarkedNode][i]>0) ) {
		if(distance[i] > distance[closestUnmarkedNode]+adjMatrix[closestUnmarkedNode][i]) {
	  		distance[i] = distance[closestUnmarkedNode]+adjMatrix[closestUnmarkedNode][i];
	  		predecessor[i] = closestUnmarkedNode;
		}
      }
    }
    count++;
  }
}

/*void Graph::printPath(int node)
{
  if(node == source)
   { cout<<node<<"..";
     fprintf(fp1,"%d\n",node); }
  else if(predecessor[node] == -1)
    cout<<"No path from "<<source<<"to "<<node<<endl;
  else {
    printPath(predecessor[node]);
    cout<<node<<"..";
    fprintf(fp1,"%d\t",node);
  }
}*/

void Graph::output()
{
  for(int i=0;i<ov1;i++) {
    if(i == source)
      cout<<source<<".."<<source;
    else
      printPath(i);
    fprintf(fp1,"\n");
    cout<<"->"<<distance[i]<<endl;
  }
}
void Graph::printPath(int node)
{

  if(node == source)
{    cout<<node<<"..";
 fprintf(fp1,"%d\t",node);
}
  else if(predecessor[node] == -1)
    cout<<"No path from "<<source<<"to "<<node<<endl;
  else {
    printPath(predecessor[node]);
    cout<<node<<"..";
fprintf(fp1,"%d\t",node);
  }
}

int main()
{
int p=0;
//int ov,ov1;
Graph G;
G.read();
//ov=numOfVertices/4;
ov1=ov;
//char[15]; 
for(p=1;p<=4;p++)
{
cout<<"\n\nLeveL";
if(p==1)
fp1=fopen("overlay1.txt","w");
else if(p==2)
fp1=fopen("overlay2.txt","w");
else if(p==3)
fp1=fopen("overlay3.txt","w");
else
fp1=fopen("overlay4.txt","w");
fprintf(fp1,"\n\n");
G.dijkstra();
//ov1=ov1+ov;
  G.output();
ov1=ov1+ov;
}
  return 0;
}
