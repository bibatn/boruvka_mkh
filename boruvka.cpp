# include "boruvka.h"

boruvkaMST::boruvkaMST(int np, int mp)
{
    np_ = np;
    mp_ = mp;
}

void boruvkaMST::Build_derived_type(Edge* indata, MPI_Datatype* message_type_ptr)
{
    int block_lengths[3];
    MPI_Aint displacements[3];
    MPI_Aint addresses[4];

    MPI_Datatype typelist[3];

    typelist[0] = MPI_INT;
    typelist[1] = MPI_INT;
    typelist[2] = MPI_DOUBLE;
/* Определить количество элементов каждого типа */

    block_lengths[0] = block_lengths[1] = block_lengths[2] = 1;

/* Вычислить смещения элементов
* относительно indata */
    MPI_Get_address(indata, &addresses[0]);

    MPI_Get_address(&(indata->l), &addresses[1]);
    MPI_Get_address(&(indata->r), &addresses[2]);
    MPI_Get_address(&(indata->w), &addresses[3]);

    displacements[0] = addresses[1] - addresses[0];
    displacements[1] = addresses[2] - addresses[0];
    displacements[2] = addresses[3] - addresses[0];

    MPI_Type_create_struct(3, block_lengths, displacements, typelist, message_type_ptr);

    MPI_Type_commit(message_type_ptr);
}

void boruvkaMST::sort(vertex_id_t * endV, double * weights, int N)
{
    std::pair<vertex_id_t , double> * pairs;
    pairs = (std::pair<vertex_id_t , double> *) malloc(N * sizeof(pairs[0]));

    for ( int i = 0; i < N; ++i )
        pairs[ i ] = std::make_pair( endV[ i ], weights[ i ] );

    std::sort( &pairs[0], &pairs[N], customLess); //

    for ( int i = 0; i < N; ++i )
    {
        endV[ i ] = pairs[ i ].first;
        weights[ i ] = pairs[ i ].second;
    }
}

void boruvkaMST::readAndSendGraph(char *fileName, int *pnVtx, int *pnEdge, double **padj)
{
    int mp, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mp);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    graph_t g;
//    char name[] = "rmat-7";
//    int l = strlen(name);
//    strncpy(inFilename,  name, (l > FNAME_LEN-1 ? FNAME_LEN-1 : l) );
////    strncpy(inFilename, "rmat-4", 6);
    readGraph(&g, inFilename);
//    printGraph(&g);


    // отсортируем вершины;

//    sort(g.endV, g.weights, g.m);
    edge_id_t* rowsIndices_;
    vertex_id_t* endV_;
    weight_t* weights_;

    rowsIndices_ = (edge_id_t *)malloc((g.n+1) * sizeof(edge_id_t));
    for(int i = 0; i < g.n + 1; i++)
        rowsIndices_[i] = g.rowsIndices[i];
    endV_ = (vertex_id_t *)malloc(g.rowsIndices[g.n] * sizeof(vertex_id_t)); //
    for(int i = 0; i < g.rowsIndices[g.n]; i++)
        endV_[i] = 0;
    weights_ = (weight_t *)malloc(g.m * sizeof(weight_t));// ну и зачемже ты обозначаешь по разному??
    for(int i = 0; i < g.m; i++)
        weights_[i] = 0;

    int k_ = 0; // счетчик для массива endV
    int q_ = 0;  // счетчик для массива rowsIndices
    sort(&g.endV[0], &g.weights[0], g.rowsIndices[1] - g.rowsIndices[0]);
    int w = g.endV[0]; // предыдущий элемент
    for(int i = 1; i < g.m; i++)
    {
        int flag = 0;
        int w_n = g.endV[i];
        if(g.rowsIndices[q_+1] <= i)
        {
            q_++;
            sort(&g.endV[g.rowsIndices[q_]], &g.weights[g.rowsIndices[q_]], g.rowsIndices[q_+1] - g.rowsIndices[q_]);
            w_n = g.endV[i];

            //случай если мы переходим на новую вершину а w_n == w
            if((w_n == w)&&(w!=q_))//HERE!!  29 29 переход
            {
                endV_[k_] = w;
                weights_[k_] = g.weights[i-1];
                k_++;
                flag = 1;
            }
        }
        if ((w_n != w)&&((w!=q_)||(g.rowsIndices[q_] == i)))
        {
            endV_[k_] = w;
            weights_[k_] = g.weights[i-1];
            k_++;
        }
        else
        {
            if(!flag)
                for(int i = q_ + 1; i < g.n + 1; i++)
                    rowsIndices_[i]--;
        }

        if(i == (g.m - 1))
        {
            if((w_n != endV_[k_-1])&&(w_n!=q_))
            {
                endV_[k_] = w_n;
                weights_[k_] = g.weights[i];
                k_++;

            }
            else
            {
                rowsIndices_[g.n]--;
            }
        }

        w = w_n;
    }

//    for(int l = g.rowsIndices[28]; l < g.rowsIndices[29]; l++)
//        std::cout << g.endV[l] << "  ";
//    std::cout << std::endl;
//
//    for(int l = rowsIndices_[28]; l < rowsIndices_[29]; l++)
//        std::cout << endV_[l] << "  ";
//    std::cout << std::endl;
//
//    for(int l = g.rowsIndices[29]; l < g.rowsIndices[30]; l++)
//        std::cout << g.endV[l] << "  ";
//    std::cout << std::endl;
//
//    for(int l = rowsIndices_[29]; l < rowsIndices_[30]; l++)
//        std::cout << endV_[l] << "  ";
//    std::cout << std::endl;

//    std::cout << "K:  " << k_ << std::endl;
//
//    std:: cout << "_________" << std::endl;
//    for(int l = 0; l < g.n + 1; l++)
//        std::cout << g.rowsIndices[l] << "  ";
//    std::cout << std::endl;
//    for(int l = 0; l < g.rowsIndices[g.n]; l++)
//        std::cout << g.endV[l] << "  ";
//    std::cout << std::endl;
//    for(int l = 0; l < g.m; l++)
//        std::cout << g.weights[l] << "  ";
//    std::cout << std::endl << std::endl;
//    std:: cout << "_________" << std::endl;

//    std:: cout << "_________" << std::endl;
//    for(int l = 0; l < g.n + 1; l++)
//        std::cout << rowsIndices_[l] << "  ";
//    std::cout << std::endl;
//    for(int l = 0; l < g.rowsIndices[g.n]; l++)
//        std::cout << endV_[l] << "  ";
//    std::cout << std::endl;
//    for(int l = 0; l < g.m; l++)
//        std::cout << weights_[l] << "  ";
//    std::cout << std::endl << std::endl;
//    std:: cout << "_________" << std::endl;

//
    int nVtx = g.n;
    int nEdge = k_ / 2;
    double *adj;
    double *tempadj;
//    std::cout << nVtx << std::endl;

    tempadj = (double *) malloc(nVtx * nVtx * sizeof(double ));

    int nb_elements = nVtx*(int)ceil((float)nVtx/(float)size);

    // Give nm_elements number of vertices' adjacencies to P-1 processes and
    // remaining to the last one.
    if (mp == size-1)
        *padj = adj = (double *) malloc((nVtx*nVtx-(size-1)*nb_elements)*sizeof(adj[0]));
    else
        *padj = adj = (double *) malloc(nb_elements*sizeof(adj[0]));

    // Let the root directly read from file
    if (mp == 0)
    {
        for (int i = 0; i < nVtx; i++)
            for (int j = 0; j < nVtx; j++)
                tempadj[i*nVtx+j] = 0.0;
//        oM.N = iM.N; // = g.n   iM.row_index[l] = g.endv[l]     iM.Value[j] =g.weights[j]      iM.Col[j] = g.rowsIndices[j]
        for (int i = 0; i < (nVtx); i++) {
            for (edge_id_t j = rowsIndices_[i]; j <= rowsIndices_[i + 1] - 1; j++) {
                double a = weights_[j];
                vertex_id_t b = endV_[j];
                tempadj[i*nVtx + b] = a;
            }
        }
//        std::cout << nVtx << std::endl;
        for (int i = 0; i < nVtx; i++) {
            for (int j = 0; j < nVtx; j++) {
//                std::cout << tempadj[i*nVtx+j] << " : ";
                if(tempadj[i*nVtx+j]!=tempadj[j*nVtx+i])
//                    std::cout << "i " << i << " j " << j << " tempadj[i*nVtx+j] " << tempadj[i*nVtx+j] << " tempadj[j*nVtx+i] " << tempadj[j*nVtx+i] << std::endl;
                tempadj[i*nVtx+j]=max(tempadj[i*nVtx+j], tempadj[j*nVtx+i]);
            }
//            std::cout << "\n";
        }

    }

    *pnVtx = nVtx;
    *pnEdge = nEdge;

    if (size > 1)
    {
        MPI_Comm commFirst;
        MPI_Comm_split(MPI_COMM_WORLD, mp/(size-1), mp, &commFirst);

        //SCATTER THE MATRIX FOR P-1 PROCESS
        if (mp != size-1)
            MPI_Scatter(tempadj, nb_elements, MPI_DOUBLE, adj, nb_elements, MPI_DOUBLE, 0, commFirst);

        //SEND THE LAST ELEMENTS TO THE LAST PROCESS
        if (mp == 0)
            MPI_Send(tempadj+(size-1)*nb_elements, nVtx*nVtx-(size-1)*nb_elements, MPI_DOUBLE, size-1, 0, MPI_COMM_WORLD);
        else if (mp == size-1)
            MPI_Recv(adj, nVtx*nVtx-(size-1)*nb_elements, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    } else {
        memcpy(adj, tempadj, nVtx*nVtx*sizeof(tempadj[0]));
    }

    free(tempadj);
}

void boruvkaMST::printUsage() {
    int mp;
    MPI_Comm_rank(MPI_COMM_WORLD, &mp);
    if (mp == 0)
    {
        printf("Usage: mpirun -np [кол-во процессов] ./boruvka [имя графа]\n");
    }
}

void boruvkaMST::write_output_information(forest_t *trees, char *filename)
{
    FILE *F = fopen(filename, "w");
    assert(fwrite(&trees->numTrees, sizeof(vertex_id_t), 1, F) == 1);
    assert(fwrite(&trees->numEdges, sizeof(edge_id_t), 1, F) == 1);
    assert(fwrite(trees->p_edge_list, sizeof(edge_id_t), 2*trees->numTrees, F) == 2*trees->numTrees);
//    for(int l = 0; l < 2*trees->numTrees; l++)
//        std::cout << trees->p_edge_list[l] << "  ";
//    std::cout << std::endl;
    assert(fwrite(trees->edge_id, sizeof(edge_id_t), trees->numEdges, F) == trees->numEdges);
    std::sort(&trees->edge_id[0], &trees->edge_id[trees->numEdges]);
//    for(int l = 0; l < trees->numEdges; l++)
//        std::cout << trees->edge_id[l] << "  ";
//    std::cout << std::endl;
    fclose(F);
}

void boruvkaMST::write_validation_file(int nVtx, Edge* tree)
{
    int N = 0;
    while(1)
    {
        if(((tree[N].w == 0)&&(tree[N].l == 0)&&(tree[N].r == 0))||(N>=(nVtx-1)))
            break;
        N++;
    }
//    std::cout << nVtx << "  " << N << "  " << (tree[N].w != 0) << (tree[N].l != 0) << (tree[N].j != 0) <<(N<nVtx) << std::endl;
    graph_t g;
    readGraph(&g, inFilename);
    forest_t trees;
    trees.numEdges = N;
    trees.numTrees = nVtx - N;
//    std::cout <<"forest_size" << nVtx - N << std::endl;
    trees.p_edge_list = (edge_id_t *) malloc(sizeof(edge_id_t) * trees.numTrees * 2);
    trees.edge_id = (edge_id_t *) malloc (sizeof (edge_id_t)*N);

    for(int n = 0; n<N; n++)
    {
        for(int m = g.rowsIndices[tree[n].l]; m<g.rowsIndices[tree[n].l+1]; m++)
        {
            if((tree[n].r==g.endV[m])&&tree[n].w==g.weights[m])
                trees.edge_id[n] = m;
        }
    }

    //костыль
    trees.p_edge_list[0] = 0;
    trees.p_edge_list[1] = N;
    int pointer = N;
    for(int i = 1; i<trees.numTrees; i++)
    {
        trees.p_edge_list[2*i] = pointer;
        trees.p_edge_list[2*i+1] = pointer;
    }


//    trees.numTrees = nVtx-N;


    if (strlen(outFilename) == 0) sprintf(outFilename, "%s.mst", inFilename);
    std::cout << outFilename << std::endl;
    write_output_information(&trees, outFilename);

//        printGraph(&g);

}

void boruvkaMST::print_tree(Edge* tree, int numEdges)
{
    double sum = 0;
    for (int edge = 0; edge < numEdges; ++edge) {
        cout << tree[edge].l << " " << tree[edge].r << "  " << tree[edge].w << endl;
        sum += tree[edge].w;
    }
    cout << "SUM : " << sum << endl;
}

void boruvkaMST::init_edge(Edge* edge, int i, int j, double w) {
    edge->l = min(i, j);
    edge->r = max(i, j);
    edge->w = w;
}

Node* boruvkaMST::create_nodes(int N) {
    Node* nodes = (Node*)calloc(N, sizeof(Node));
    for (int node = 0; node < N; ++node) {
        nodes[node].parent = node;
        nodes[node].rank = 1;
    }
    return nodes;
}

int boruvkaMST::find(Node* nodes, int node) {
    if (nodes[node].parent != node) {
        nodes[node].parent = find(nodes, nodes[node].parent);
    }
    return nodes[node].parent;
}

void boruvkaMST::merge_trees(Node* nodes, int root1, int root2) {

    if (nodes[root1].rank > nodes[root2].rank) {
        merge_trees(nodes, root2, root1);
        return ;
    }

    nodes[root1].parent = root2;
    if (nodes[root1].rank == nodes[root2].rank) {
        nodes[root2].rank += 1;
    }
}

int boruvkaMST::serial_boruvka(Edge* edges, int N, int M, Edge *tree) {
    int numEdges = 0;

    Node* nodes = create_nodes(N);

    int *cheapest = new int[N];

    bool canMerge = true;

    while (canMerge)
    {
        canMerge = false;
        for(int i=0; i<N; i++) cheapest[i] = -1;

        for (int edge=0; edge<M; edge++)
        {
            int root1 = find(nodes, edges[edge].l);
            int root2 = find(nodes, edges[edge].r);

            if (root1 == root2)
                continue;

            else
            {
                if (cheapest[root1] == -1 ||
                    edges[cheapest[root1]].w > edges[edge].w)
                    cheapest[root1] = edge;

                if (cheapest[root2] == -1 ||
                    edges[cheapest[root2]].w > edges[edge].w)
                    cheapest[root2] = edge;
            }
        }

        for (int node=0; node<N; node++)
        {
            // Check if cheapest for current set exists
            if (cheapest[node] != -1)
            {
                int root1 = find(nodes, edges[cheapest[node]].l);
                int root2 = find(nodes, edges[cheapest[node]].r);

                if (root1 == root2)
                    continue;

                canMerge = true;

                memcpy(tree + numEdges++, edges + cheapest[node], sizeof(Edge));

                merge_trees(nodes, root1, root2);
            }
        }
    }

    free(nodes);

    return numEdges;
}

void boruvkaMST::merge_lists(Edge* li, int ni, Edge* lj, int nj, Edge* output) {
    int i = 0;
    int j = 0;
    while (i < ni) {
        memcpy(output+i, li+i, sizeof(Edge));
        i += 1;
    }
    while(j < nj) {
        memcpy(output+i+j, lj+j, sizeof(Edge));
        j += 1;
    }
}

void boruvkaMST::merge_into(Edge* oldBuf, int n1, Edge* newBuf, int n2) {
    Edge* output = (Edge* )calloc(n1+n2, sizeof(Edge));
    merge_lists(oldBuf, n1, newBuf, n2, output);
    memcpy(newBuf, output, (n1+n2)*sizeof(Edge));
    free(output);
}

int boruvkaMST::receive_edges_from(int stepSize, int targetRank, Edge* edges) {
    int numEdges = 0;

    // Receive number of edges.
    MPI_Recv(&numEdges, 1, MPI_INT, targetRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    Edge * buffer = (Edge*)calloc(numEdges, sizeof(Edge));
    MPI_Datatype message_type;
    Build_derived_type(buffer, &message_type);

    MPI_Recv(buffer, numEdges, message_type, targetRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = 0; i < numEdges; ++i) {
        edges[i] = buffer[i];
    }

    return numEdges;
}

int boruvkaMST::receive_edges(int procRank, int numProcs, int stepSize, int squareSize, Edge* edges,
                  int numEdges, int numEdgesMax) {

    Edge* tmpBuf = (Edge* )calloc(numEdgesMax, sizeof(Edge));
    Edge* buf1 = edges;
    Edge* buf2 = tmpBuf;

    Edge* newEdges = (Edge* )calloc(stepSize * squareSize, sizeof(Edge));
    int parity = 0;
    for (int i = 0; i < stepSize; ++i) {
        if (procRank + stepSize + i >= numProcs) break;
        int numNewEdges = receive_edges_from(stepSize, procRank + stepSize + i, newEdges);
        merge_lists(buf1, numEdges, newEdges, numNewEdges, buf2);
        numEdges += numNewEdges;

        Edge* tmp = buf1;
        buf1 = buf2;
        buf2 = tmp;

        parity = 1 - parity;
    }

    if (parity != 0) {
        memcpy(edges, tmpBuf, numEdges*sizeof(Edge));
    }

    free(newEdges);
    free(tmpBuf);

    return numEdges;
}

int boruvkaMST::receive_forest(int procRank, int numProcs, int stepSize, Edge* edges) {
    int target = procRank + stepSize;
    if (target >= numProcs)
        return 0;
    return receive_edges_from(stepSize, target, edges);
}

int boruvkaMST::receive_new_edges(
        int procRank, int numProcs, int stepSize, int squareSize, int numRows,
        Edge* edges, int numMaxEdges, Edge* forest, int forestSize) {

    int newForestSize = receive_forest(procRank, numProcs, stepSize, edges);

    int numNewEdges = receive_edges(procRank, numProcs, stepSize, squareSize,
                                    edges, newForestSize, numMaxEdges);

    merge_into(forest, forestSize, edges, numNewEdges);

    numNewEdges += forestSize;

    return numNewEdges;
}

int boruvkaMST::add_bipartite_edges(int procRank, int stepSize, double* adj, int numRows, int N, Edge* edges) {
    int start = ((procRank - procRank%stepSize) - stepSize)*numRows;
    int end = start + numRows*stepSize;

    int numEdges = 0;
    for (int i = 0; i < numRows; ++i) {

        int procIdx = procRank*numRows + i;

        if (procIdx >= N) break;

        for (int j = start; j < end; ++j) {
            if (adj[i*N + j] == 0) continue ;
            init_edge(edges+numEdges, procIdx, j, adj[i*N + j]);
            numEdges += 1;
        }
    }
    return numEdges;
}

int boruvkaMST::add_local_edges(int procRank, double* adj, int numRows, int N, Edge* edges) {
    int numEdges = 0;

    for (int i = 0; i < numRows; ++i) {

        int procIdx = procRank*numRows + i;

        if (procIdx >= N) break;

        for (int j = procRank*numRows; j <= procIdx; ++j) {
            if (adj[i*N + j] == 0) continue ;
            init_edge(edges+numEdges, procIdx, j, adj[i*N + j]);
            numEdges += 1;
        }
    }
    return numEdges;
}

void boruvkaMST::send_edges_to(int target, Edge* edges, int numEdges) {

    Edge * buf = (Edge*)calloc(numEdges, sizeof(Edge));
    for (int i = 0; i < numEdges; ++i) {
        buf[i] = edges[i];
    }
    MPI_Datatype message_type;
    Build_derived_type(buf, &message_type);

    MPI_Send(&numEdges, 1, MPI_INT, target, 0, MPI_COMM_WORLD);

    MPI_Send(buf, numEdges, message_type, target, 0, MPI_COMM_WORLD);

    free(buf);
}



void boruvkaMST::send_bipartite_forest(int procRank, int stepSize, int squareSize, int numRows,
                           double* adj, int N) {
    int target = (procRank - procRank%stepSize) - stepSize;

    Edge* edges = (Edge* )calloc(squareSize*stepSize*stepSize, sizeof(Edge));

    int numEdges = add_bipartite_edges(procRank, stepSize, adj, numRows, N, edges);

    Edge* forest = (Edge* )calloc((stepSize+1)*numRows-1, sizeof(Edge));

    int forestSize = serial_boruvka(edges, N, numEdges, forest);
    send_edges_to(target, forest, forestSize);

    free(forest);
    free(edges);
}

void boruvkaMST::send_forest(int procRank, int stepSize, Edge* forest, int forestSize)
{
    int target = procRank - stepSize;

    send_edges_to(target, forest, forestSize);
}

int boruvkaMST::create_forest(int procRank, double *adj, int numRows, int N, Edge* forest)
{
    Edge* edges = (Edge* )calloc(numRows*numRows, sizeof(Edge));

    int numEdges = add_local_edges(procRank, adj, numRows, N, edges);

    int forestSize = serial_boruvka(edges, N, numEdges, forest);

    free(edges);

    return forestSize;
}

void boruvkaMST::parallel_boruvka(int procRank, int numProcs, double *adj, int N, int M, Edge* tree) {
    // Изначальное количество вершина на процесс.
    int numRows = ceil((float)N / (float)numProcs);
    int squareSize = numRows * numRows;

    // Получение остовного дерева по тем вершинам что выданы.
    Edge* forest = (Edge* )calloc(numRows - 1, sizeof(Edge));
    int forestSize = create_forest(procRank, adj, numRows, N, forest);

    // Каждый процесс начинает как получатель
    int receiver = 1;


    for (int stepSize = 1, rank = procRank; stepSize*numRows < N; stepSize <<= 1, rank >>= 1) {

        if (rank & 1) {
            receiver = 0;

            if (procRank % stepSize == 0) {
                send_forest(procRank, stepSize, forest, forestSize);
            }

            send_bipartite_forest(procRank, stepSize, squareSize, numRows, adj, N);

        } else if (receiver) {
            int subnumEdges = stepSize*(numRows*(stepSize+1)-1)  + 2*(stepSize*numRows-1);
            int numMaxEdges = (M < subnumEdges) ? M : subnumEdges;
            Edge* edges = (Edge* )calloc(numMaxEdges, sizeof(Edge));
            int numNewEdges = receive_new_edges(procRank, numProcs, stepSize, squareSize,
                                                numRows, edges, numMaxEdges, forest, forestSize);
            forest = (Edge* )realloc(forest, (stepSize*numRows*2 - 1) * sizeof(Edge));

            // Create a spanning forest of all the edges received.
            forestSize = serial_boruvka(edges, N, numNewEdges, forest);
            free(edges);
        }
    }
    if (procRank == 0) {
        memcpy(tree, forest, sizeof(Edge) * forestSize);
    }
    free(forest);
}

void boruvkaMST::computeMST(
        int N,
        int M,
        double *adj,
        Edge *tree)
{
    int mp, np;
    MPI_Comm_rank(MPI_COMM_WORLD, &mp);
    MPI_Comm_size(MPI_COMM_WORLD, &np);


    parallel_boruvka(mp, np, adj, N, M, tree);

//    if (mp == 0) {
//        print_tree(tree, N-1);
//    }
}


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

    MPI_Barrier(MPI_COMM_WORLD);

    Edge* tree;
    tree = (Edge* )calloc(nVtx-1, sizeof(Edge));
    std::vector<double> times;
    int n = 100; // кол-во измерений
    for(int i = 0; i<100; i++)
    {
        double startTime = MPI_Wtime();
        boruvka.computeMST(nVtx, nEdge, adj, tree);
        times.push_back(MPI_Wtime() - startTime);
    }

    if (mp == 0) {
        printf("Time: %e seconds.\n", std::accumulate(times.begin(), times.end(), 0.0) / times.size());
        boruvka.write_validation_file(nVtx, tree);
    }

    MPI_Finalize();

    free(adj);
    return 0;
}

