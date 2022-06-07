#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include<mpi/mpi.h>
#include <math.h>
#include "defs.h"
#include <cassert>
#include <algorithm>
#include <vector>
#include <numeric>
#include "seq_generation_utils.h"

char inFilename[FNAME_LEN];
char outFilename[FNAME_LEN];

using namespace std;

typedef struct Edge {
    int l;
    int r;
    double w;
} Edge;

typedef struct Node {
    int parent;
    int rank;
} Node;

// Кастомный компаратор
struct {
    bool operator()(std::pair<vertex_id_t , double> a, pair<vertex_id_t , double> b) const
    {
        if(a.first < b.first)
        {
            return 1;
        }
        else
        {
            if((a.first == b.first)&&(a.second > b.second))
                return 1;
            else
                return 0;
        }
    }
} customLess;

class boruvkaMST
{
public:
    boruvkaMST(int np, int mp);

    void sort(vertex_id_t * endV, double * weights, int N);
    void Build_derived_type(Edge* indata, MPI_Datatype* message_type_ptr);
    void readAndSendGraph( char *fileName, int *pnVtx, int *pnEdge, double **padj);
    void printUsage();
    void write_validation_file(int nVtx, Edge* tree);
    void write_output_information(forest_t *trees, char *filename);
    void print_tree(Edge* tree, int numEdges);

    void init_edge(Edge* edge, int i, int j, double w);
    Node* create_nodes(int N);
    int find(Node* nodes, int node);
    void merge_trees(Node* nodes, int root1, int root2);
    int serial_boruvka(Edge* edges, int N, int M, Edge *tree);
    void merge_lists(Edge* li, int ni, Edge* lj, int nj, Edge* output);
    void merge_into(Edge* oldBuf, int n1, Edge* newBuf, int n2);
    int receive_edges_from(int stepSize, int targetRank, Edge* edges);
    int receive_edges(int procRank, int numProcs, int stepSize, int squareSize, Edge* edges,
                      int numEdges, int numEdgesMax);
    int receive_forest(int procRank, int numProcs, int stepSize, Edge* edges);
    int receive_new_edges(
            int procRank, int numProcs, int stepSize, int squareSize, int numRows,
            Edge* edges, int numMaxEdges, Edge* forest, int forestSize);
    int add_bipartite_edges(int procRank, int stepSize, double* adj, int numRows, int N, Edge* edges);
    int add_local_edges(int procRank, double* adj, int numRows, int N, Edge* edges);
    void send_edges_to(int target, Edge* edges, int numEdges);
    void send_bipartite_forest(int procRank, int stepSize, int squareSize, int numRows,
                               double* adj, int N);
    void send_forest(int procRank, int stepSize, Edge* forest, int forestSize);
    int create_forest(int procRank, double *adj, int numRows, int N, Edge* forest);
    void parallel_boruvka(int procRank, int numProcs, double *adj, int N, int M, Edge* tree);
    void computeMST( int N, int M, double *adj, Edge *tree);

private:
    int np_, mp_;
};
