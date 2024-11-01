#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <limits.h>
#include <stddef.h>

#define INF INT_MAX
//儲存0到各個點的編號與其對應的最短距離
typedef struct{
    int No;
    int dist;
}Node;

//記錄此條邊可到的點與其權重
typedef struct{
    int to;
    int weight;
}Edge;

typedef struct{
    Node *data;
    int size;
    int capacity;
}PriorityQueue;

PriorityQueue* createQueue(int capacity){
    PriorityQueue* pq = (PriorityQueue*)malloc(sizeof(PriorityQueue));
    pq -> data = (Node*)malloc(capacity * sizeof(Node));
    pq -> size = 0;
    pq -> capacity = capacity;
    return pq;
}

void swap(Node *a, Node *b){
    Node temp = *a;
    *a = *b;
    *b = temp;
}
// accroding to dist to make a minheap
void minHeapify(PriorityQueue *pq, int idx){
    int smallest = idx;
    int left = 2 * idx + 1;
    int right = 2 * idx + 2;

    if(left < pq -> size && pq -> data[left].dist < pq -> data[smallest].dist){
        smallest = left;
    }
    if(right < pq -> size && pq -> data[right].dist < pq -> data[smallest].dist){
        smallest = right;
    }
    if(smallest != idx){
        swap(&pq -> data[smallest], &pq -> data[idx]);
        minHeapify(pq, smallest);
    }
}

void push(PriorityQueue *pq, int No, int dist){
    //before inserting neew node, check the capacity. 
    //If size == capacity, then capacity is doubled. 
    //This can prevents insertion failures when the queue size reaches its maximum limit.
    if(pq -> size == pq -> capacity){ 
        pq -> capacity *= 2;
        pq -> data = (Node*)realloc(pq -> data, pq -> capacity * sizeof(Node));
    }
    int i = pq -> size++;
    pq -> data[i].No = No;
    pq -> data[i].dist = dist;
    //minheap -> index of parent node = (index of child node - 1) / 2
    //so compare with parent need to find the data in index (i - 1) / 2 
    while(i > 0 && pq -> data[i].dist < pq -> data[(i - 1) / 2].dist){
        swap(&pq -> data[i], &pq -> data[(i - 1) / 2]);
        i = (i - 1) / 2;
    }
}

Node pop(PriorityQueue *pq){
    Node root = pq -> data[0]; //in minheap, the uppest(first) node is root //will be remove
    pq -> data[0] = pq -> data[(pq -> size) - 1]; //move the latest node to first to replace the root
    pq -> size--;
    minHeapify(pq, 0); //run for build new heapify 
    return root;
}

int isEmpty(PriorityQueue *pq){
    if(pq -> size == 0) return 1; 
    else return 0;
}
int dijkstra(int n, int start, int end, Edge** graph, int* dist, int *edge_count, int rank){
    int update = 0;
    PriorityQueue *pq = createQueue(end - start);
    for(int i = start; i < end; i++){
        push(pq, i, dist[i]);
    }

    while (isEmpty(pq) == 0){ //while pq is not empty, run this
        Node node = pop(pq);
        int u = node.No;
        for(int i = 0; graph[u] != NULL && i < edge_count[u]; i++){
            Edge edge = graph[u][i];
            int v = edge.to;
            int weight = edge.weight;
            if(dist[u] != INF && (dist[u] + weight) < dist[v]){ //u can go to v, so (weight + dist[u]) vs. dist[v]
                dist[v] = dist[u] + weight;
                if(start <= v && v < end) push(pq, v, dist[v]);
                update = 1;
            }
        }
    }
    if(pq -> size != 0) free(pq -> data);
    free(pq);
    return update;
}

int main(){
    int size, rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int n;
    int *dist;
    int *cut_num = (int*)malloc(size * sizeof(int));
    int *start;
    int *end;
    Edge *gr;
    int sum;
    Edge **graph;
    int *edge_count;
    Edge **g;
    MPI_Datatype mpi_edge;
    MPI_Type_create_struct(2, (int[]){1, 1}, (MPI_Aint[]){offsetof(Edge, to), offsetof(Edge, weight)}, (MPI_Datatype[]){MPI_INT, MPI_INT}, &mpi_edge);
    MPI_Type_commit(&mpi_edge);
    
    
    if(rank == 0){
        char filename[50];
        scanf("%s", filename);
        FILE *file = fopen(filename, "r");
        if(file == NULL){
            printf("Could not open the file %s", filename);
            return 1;
        }
        fscanf(file, "%d", &n);
        graph = (Edge**)malloc(n * sizeof(Edge*)); //graph儲存各個點的邊數與其邊的相關資訊 //graph[i][j] -> i為起點
        for(int i = 0; i < n; i++){
            graph[i] = (Edge*)malloc(sizeof(Edge));
        }
        edge_count = (int*)calloc(n, sizeof(int)); //計算每個點各有幾條邊，且皆初始化為0
        int a, b, c;
        while(fscanf(file, "%d %d %d ", &a, &b, &c) != EOF){ //EOF：檔案結束符號
            //若已有k條邊，則記憶體+1乙存方當前新資訊
            graph[a] = (Edge*)realloc(graph[a], (edge_count[a] + 1) * sizeof(Edge)); //動態分配記憶體 //針對graph[a]將其大小改為其對應的edge_count + 1(因為要加新的邊)
            //目前應有k + 1條邊，但edge_count未更新，故仍為k。又k + 1應記錄在 0 ~ k，故直接使用edge_count
            graph[a][edge_count[a]].to = b;
            graph[a][edge_count[a]].weight = c;
            //此處再更新edge_count值
            edge_count[a]++;
        }
        fclose(file);
        sum = 0;
        for(int i = 0; i < n; i++) sum += edge_count[i];
        
        
        for(int i = 0; i < size; i++){
            cut_num[i] = n / size;
            if(i < (n % size)) cut_num[i]++;
        }
        start = (int*)malloc(size * sizeof(int));
        end = (int*)malloc(size * sizeof(int));
        start[0] = 0;
        end[0] = cut_num[0];
        
        for(int i = 1; i < size; i++){
            start[i] = end[i - 1];
            end[i] = end[i - 1] + cut_num[i];
        }
        free(cut_num);
    }    
    int s;
    int e;
    MPI_Bcast(&sum, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    gr = (Edge*)malloc(sum * sizeof(Edge));
    int *e_count = (int*)malloc(n * sizeof(int));
    if(rank == 0){
        int k = 0;
        for(int i = 0; i < n; i++){
            e_count[i] = edge_count[i];
            for(int j = 0; j < edge_count[i]; j++){
                gr[k].to = graph[i][j].to;
                gr[k].weight = graph[i][j].weight;
                k++;   
            }
        }
        free(graph);
        free(edge_count);
    }
    
    MPI_Bcast(e_count, n, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(start, 1, MPI_INT, &s, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(end, 1, MPI_INT, &e, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(gr, sum, mpi_edge, 0, MPI_COMM_WORLD);
    free(start);
    free(end);
    
    g = (Edge**)malloc(n * sizeof(Edge*));

    int dis = 0;
    for(int i = 0; i < n; i++){
        g[i] = (Edge*)malloc(e_count[i] * sizeof(Edge));
        for(int j = 0; j < e_count[i]; j++){
            g[i][j] = gr[dis + j];
        } 
        dis += e_count[i];
    }
    free(gr);
    
    dist = (int*)malloc(n * sizeof(int));
    if(rank == 0){
        dist[0] = 0;
        for(int i = 1; i < n; i++){
            dist[i] = INF;
        }
    }
    MPI_Bcast(dist, n, MPI_INT, 0, MPI_COMM_WORLD);
    int update = 1;
    MPI_Barrier(MPI_COMM_WORLD);
    int displs;
    while(update){
        update = dijkstra(n, s, e, g, dist, e_count, rank);
        MPI_Allreduce(MPI_IN_PLACE, &update, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if(update) MPI_Allreduce(MPI_IN_PLACE, dist, n, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    }
    free(g);
    free(e_count);

    for(int r = 0; r < size; r++){
        if(r == rank){
            for(int i = s; i < e; i++){
                printf("%d ", dist[i]);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    free(dist);
    return 0;
}