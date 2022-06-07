# include "boruvka.h"

int main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);
    int mp, np;
    MPI_Comm_rank(MPI_COMM_WORLD, &mp);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    boruvkaMST boruvka(np, mp);

    if (argc < 2) {
        boruvka.printUsage();
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int nVtx, nEdge;
    double *adj = NULL;

    int l = strlen(argv[1]);
    strncpy(inFilename,  argv[1], (l > FNAME_LEN-1 ? FNAME_LEN-1 : l) );

    boruvka.readAndSendGraph(argv[1], &nVtx, &nEdge, &adj);

    // Wait until all the process have read the graph
    MPI_Barrier(MPI_COMM_WORLD);

    double startTime = MPI_Wtime();
    Edge* tree;
    tree = (Edge* )calloc(nVtx-1, sizeof(Edge));
    boruvka.computeMST(nVtx, nEdge, adj, tree);

    if (mp == 0) {
        printf("computeMST took %e seconds.\n", MPI_Wtime() - startTime);
        boruvka.write_validation_file(nVtx, tree);
    }

    MPI_Finalize();

    free(adj);
    return 0;
}
