#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void exchange(int n, int m, int D, int *recv, int size, int rank, MPI_Comm comm){
    MPI_Status status;
    int up;
    if(rank == 0) up = size - 1;
    else up = rank - 1;
    int down;
    if(rank == size - 1) down = 0;
    else down = rank + 1;
    //printf("Before Sendrecv rank = %d\n", rank);

    //前 上 中 下 尾
    //把上端數據傳給上一位(up), 並在尾端接收下一位(down)的上端數據
    MPI_Sendrecv(&recv[(D - 1) / 2 * m], ((D - 1) / 2) * m, MPI_INT, up, 0, &recv[(n - ((D - 1) / 2)) * m], ((D - 1) / 2) * m, MPI_INT, down, 0, comm, &status);
    //把下端數據傳給下一位(down), 並在前端接收上一位(up)的下端數據
    MPI_Sendrecv(&recv[(n - (D - 1)) * m], ((D - 1) / 2) * m, MPI_INT, down, 1, &recv[0], ((D - 1) / 2) * m, MPI_INT, up, 1, comm, &status);
}

void run(int dn, int m, int D, int *K, int *recv, int *store1){
    for(int i = ((D - 1) / 2); i < dn - ((D - 1) / 2); i++){
        for(int j = 0; j < m; j++){
            store1[i * m + j] = 0;
            for(int di = -((D - 1) / 2); di <= (D - 1) / 2; di++){
                for(int dj =  -((D - 1) / 2); dj <= (D - 1) / 2; dj++){
                    int p_i = i + di;
                    int p_j = j + dj;
                    if(p_i < 0) p_i += dn;
                    else p_i = p_i % dn;
                    if(p_j < 0) p_j += m;
                    else p_j = p_j % m;
                    store1[i * m + j] += K[((D - 1) / 2 + di) * D + (D - 1) / 2 + dj] * recv[p_i * m + p_j];
                }
            }

            store1[i* m + j] = store1[i * m + j] / D / D;
        }
    }
}
int main(){
    int size, rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int times;
    int n, m;
    int **A;
    int D;
    int *sendcount;
    int *send_dis;
    int *recvbuf; 
    int *store0; 
    int *store1; 
    int *ans; 
    int *p;
    int *K; 
    
    if(rank == 0){
        char filename[50];
        scanf("%s", filename);
        FILE *file = fopen(filename, "r");
        if(file == NULL){
            printf("Could not open the file %s\n", filename);
            return 1;
        }

        times = 0;
        fscanf(file, "%d", &times);
        fscanf(file, "%d %d", &n, &m);
        A = (int**)malloc(n * sizeof(int*));
        for(int i = 0; i < n; i++){
            A[i] = (int *)malloc(m * sizeof(int));
            for(int j = 0; j < m; j++){
                fscanf(file, "%d", &A[i][j]);
            }
        }

        fscanf(file, "%d", &D);
        p = (int*)malloc(D * D * sizeof(int));
        for(int i = 0; i < D; i++){
            for(int j = 0; j < D; j++){
                fscanf(file, "%d", &p[i * D + j]);
            }
        }
        fclose(file);

        int *c_size = (int*)malloc(size * sizeof(int));
        int pre_size = 0;
        send_dis = (int*)malloc(size * sizeof(int));
        send_dis[0] = 0;
        
        sendcount = (int*)malloc(size * sizeof(int));
        store0 = (int*)malloc((n + size * (D - 1)) * m * sizeof(int));
        for(int k = 0; k < size; k++){
            c_size[k] = n / size;
            if((n % size) > k) c_size[k]++;
         
            if(k > 0) pre_size += c_size[k - 1]; //下一組該抓的數據(未包含預備空間)
            
            sendcount[k] = (c_size[k] + (D - 1)) * m;
            if(k > 0) send_dis[k] = send_dis[k - 1] + sendcount[k - 1];
            for(int i = 0; i < c_size[k] + (D - 1); i++){
                for(int j = 0; j < m; j++){
                    store0[send_dis[k] + i * m + j] = 0;
                }
            }
            for(int i = (D - 1) / 2; i < c_size[k] + (D - 1) / 2; i++){
                for(int j = 0; j < m; j++){
                    store0[send_dis[k] + i * m + j] = A[pre_size + i - (D - 1) / 2][j];
                }
            }
        }
        free(A);
        free(c_size);
    }
    
    MPI_Bcast(&times, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&D, 1, MPI_INT, 0, MPI_COMM_WORLD);
    K = (int*)malloc(D * D * sizeof(int));
    if(rank == 0){
        for(int i = 0; i < D; i++){
            for(int j = 0; j < D; j++){
                K[i * D + j] = p[i * D + j];
            }
        }
        free(p);
    }
    MPI_Bcast(K, D * D, MPI_INT, 0, MPI_COMM_WORLD);
    
    int cut_size = n / size;
    if((n % size) > rank) cut_size++;
    
    recvbuf = (int*)malloc((cut_size + (D - 1)) * m  * sizeof(int));
    store1 = (int*)malloc((cut_size + (D - 1)) * m * sizeof(int));
    ans = (int*)malloc(cut_size * m* sizeof(int));
    

    int *buf = (int*)malloc((cut_size + (D - 1)) * m * sizeof(int));
    for(int i = 0; i < (cut_size + (D - 1)) * m; i++){
        buf[i] = 0;
    }
    
    MPI_Scatterv(store0, sendcount, send_dis, MPI_INT, buf, (cut_size + (D - 1)) * m, MPI_INT, 0, MPI_COMM_WORLD);
    for(int i = 0; i < (cut_size + (D - 1)); i++){
        for(int j = 0; j < m; j++){
            recvbuf[i * m + j] = 0;
            recvbuf[i * m + j] = buf[i * m + j];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(buf);
    if(rank == 0){
        free(sendcount);
        free(send_dis);
        free(store0);
    }
    
    while(times > 0){
        for(int i = 0; i < (cut_size + (D - 1)); i++){
            for(int j = 0; j < m; j++){
                store1[i * m + j] = 0;
            }
        }
        exchange((cut_size + (D - 1)), m, D, recvbuf, size, rank, MPI_COMM_WORLD);
        // 同步，確保所有rank都執行完 exchange再繼續
        MPI_Barrier(MPI_COMM_WORLD);
        run((cut_size + (D - 1)), m, D, K, recvbuf, store1);
        
        for(int i = 0; i < (cut_size + (D - 1)); i++){
            for(int j = 0; j < m; j++){
                recvbuf[i * m + j] = store1[i * m + j];
            }
        }
        // 同步，確保所有rank執行完 run 和 更新recvbuf 再繼續
        MPI_Barrier(MPI_COMM_WORLD);
        times--;
    }
    free(K);
    free(recvbuf);

    for(int i = ((D - 1) / 2); i < cut_size + ((D - 1) / 2); i++){
        for(int j = 0; j < m; j++){
            ans[(i - ((D - 1) / 2)) * m + j]= store1[i * m + j];
        }
    }
    free(store1);
    int count = cut_size * m;
    //直接在各個rank output
    //若gatherv的話，記憶體不足無法全數彙整印出
    for(int r = 0; r < size; r++){
        //當r = 當前rank時，執行輸出
        if(rank == r){
            for(int i = 0; i < cut_size; i++){
                for(int j = 0; j < m; j++){
                    printf("%d ", ans[i * m + j]);
                }
            }
        }
        //每當r變動時，皆會執行barrier等待所有rank確仁完是否符合rank = r後，才執行r++
        //意即當r = 0時，若rank 3先執行，因r != rank，故停在barreier，所有rank皆到此時，rank 0已完成輸出，繼續r++
        //由此便可依序輸出各個rank的運算結果
        MPI_Barrier(MPI_COMM_WORLD);
    }

    free(ans);
    return 0;
}