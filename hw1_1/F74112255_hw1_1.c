#include <stdio.h> 
#include <stdlib.h>
#include <stdint.h>
#include <mpi.h>

uint64_t count_num(int n, int m, uint64_t *cover_part, uint64_t start, uint64_t end, uint64_t full_mask){
    uint64_t sum = 0;
    for(uint64_t i = start; i < end; i++){
        uint64_t check = 0;
        for(int j = 0; j < m; j++){
            //若只寫i & (1 << j) != 0的話，會先測試(1 << j) != 0 後再和i &，並判斷是否 > 0，故i & (1 << j)外層需括號
            if((i & (1 << j)) != 0){ // if i = 5(binary 101), j = 1, i & (1 << j) = 101 & 010 = 000，fail; if j = 2, i & (1 << j) = 101 & 100 = 100(binary) = 4, access. 
                check |= cover_part[j]; //若二進位制的i和 1 << j 的 "&" != 0,表示此測試有被選用，將拿來合併覆蓋區域
            }
        }
        if(check == full_mask){
            sum++;
        }
    }
    return sum;
}
int main(){
    int rank, size;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n, m;
    uint64_t *cover_part;
    
    if(rank == 0){ //只在rank 0讀取
        char filename[50];
        scanf("%s", filename);
        FILE *file = fopen(filename, "r");
        
        if(file == NULL){
            printf("Could not open the file %s\n", filename);
            return 1;
        }
         //n parts, m tests
        fscanf(file, "%d %d", &n, &m);
        int num[m]; //num of covered
        int cost[m];
        cover_part = (uint64_t*) malloc(m * sizeof(uint64_t));
        int store;

        for(int i = 0; i < m; i++){
            fscanf(file, "%d %d", &num[i], &cost[i]);
            cover_part[i] = 0;
            for(int j = 0; j < num[i]; j++){
                fscanf(file, "%d", &store); 
                cover_part[i] |= 1UL << (store - 1);
                // if store = 3 => 1 << (3 - 1) => 00000...100
                // if cover_part[i] = 00000...00000, than become 00000...100
                // if store = 5 => 1 << (5 - 1) => 00000...100000
                // than cover_part[i] = 00000...10100
            }
        }
        fclose(file);
    }
    
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);


    if(rank != 0){
        cover_part = (uint64_t*) malloc(m * sizeof(uint64_t));
    }

    MPI_Bcast(cover_part, m, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    
    //要公開給其他id使用，故放在if(rank == 0)外
    uint64_t full_mask = (1UL << n) - 1; // if n = 5, full_mask = 100000 - 1 = 011111 (binary)
    uint64_t test_sum = 1UL << m; //共2^m種測試
    
    //每個id要跑的數量
    uint64_t chunk_size = test_sum / size;
    uint64_t start = rank * chunk_size;
    uint64_t end;
    if(rank == (size - 1)) end = test_sum;
    else  end = (rank + 1) * chunk_size;

    
    
    uint64_t count = count_num(n, m, cover_part, start, end, full_mask);
    uint64_t total = 0;

    MPI_Reduce(&count, &total, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
   
    if(rank == 0){
        printf("%lu", total);
        free(cover_part);
    }
    MPI_Finalize();
    return 0;
}
