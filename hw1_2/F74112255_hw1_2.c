#include <stdio.h> 
#include <stdlib.h>
#include <mpi.h>
#include <stddef.h>

typedef struct{
    int No;
    int x;
    int y;
}Point;

//內積，若 < 0，表示oA向量 到 OB向量為順時針旋轉(右轉)
long long cross(Point o, Point a, Point b){
    return ((long long)(a.x - o.x) * (b.y - o.y)) - ((a.y - o.y) * (b.x - o.x));
}

//x升序排列，x相同則依y升序排序
int cmp(const void *a, const void *b){
    Point p1 = *(Point *)a;
    Point p2 = *(Point *)b;
    if(p1.x == p2.x){
        return p1.y -p2.y;
    }
    return p1.x - p2.x;
}

//凸包
int convex_hull(Point *point, int n, Point *ans){
    qsort(point, n, sizeof(Point), cmp);
    //上凸包
    int k = 0; //計算包含的點的數量
    for(int i = 0; i < n; i++){
        while(k >= 2 && (cross(ans[k - 2], ans[k - 1], point[i]) >= 0)){
            k--;
        }

        ans[k] = point[i];
        k++;
    }

    int t = k + 1; //k為上凸包含的點的數量
    for(int i = n - 2; i >= 0; i--){// i = n - 2 => 不重複比較最右上的點
        while(k >= t && (cross(ans[k - 2], ans[k - 1], point[i]) >= 0)){
            k--;
        }

        ans[k] = point[i];
        k++;
    }

    return k - 1;//扣除起點
}

int main(){
    int rank, size;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Datatype mpi_point;
    MPI_Type_create_struct(3, (int[]){1, 1, 1}, (MPI_Aint[]){offsetof(Point, No), offsetof(Point, x), offsetof(Point, y)}, (MPI_Datatype[]){MPI_INT, MPI_INT, MPI_INT}, &mpi_point);
    MPI_Type_commit(&mpi_point);
    
    int n = 0;
    Point *point;
    if(rank == 0){
        char filename[50];
        scanf("%s", filename); // 不是scanf("%s", &filename)
        FILE *file = fopen(filename, "r");

        if(file == NULL){
            printf("Could not open the file %s\n", filename);
            return 1;
        }

        fscanf(file, "%d", &n);
        point = (Point *)malloc(n * sizeof(Point));
        for(int i = 0; i < n; i++){
            point[i].No = i + 1;
            fscanf(file, "%d %d", &point[i].x, &point[i].y);
        }
        fclose(file);
        qsort(point, n, sizeof(Point), cmp);
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int base_size = n / size; //每個id處理的點的數量
    int mod = n % size;
    int *send_size = (int*)malloc(size * sizeof(int));
    int *dis = (int*)malloc(size * sizeof(int));
    int chunk_size = 0;
    if(rank == 0){
        dis[0] = 0;
        for(int i = 0; i < size; i++){
            send_size[i] = base_size;
            if(i < (n % size)) send_size[i]++;
            if(i > 0) dis[i] = dis[i - 1] + send_size[i - 1];
        }
    }
    if(rank < mod) chunk_size = base_size + 1;
    else chunk_size = base_size;
    
    Point *chunk_point = (Point *)malloc(chunk_size * sizeof(Point));

    //將點分發給各個進程，每個進程接收各自的點
    MPI_Scatterv(point, send_size, dis, mpi_point, chunk_point, chunk_size, mpi_point, 0, MPI_COMM_WORLD); //每個點有2個座標值，故 * 2


    //各進程的凸包
    Point *hull_point = (Point *)malloc(chunk_size * sizeof(Point));
    //各進程的凸包點數數量
    int hull_size = convex_hull(chunk_point, chunk_size, hull_point);  
    
    //總和每個進程的凸包點
    Point *all_hull_point = NULL;
    //總和每個進程的凸包點數量
    int *recvpoint = NULL;


    if(rank == 0){
        recvpoint = (int*)malloc(size * sizeof(int));
        all_hull_point = (Point*)malloc(n * sizeof(Point));
    }


    //id 0 接收各個進程回傳的凸包點數數量
    MPI_Gather(&hull_size, 1, MPI_INT, recvpoint, 1, MPI_INT, 0, MPI_COMM_WORLD);

    
    
    int *displs;
    if(rank == 0){
      displs = (int*)malloc(size * sizeof(int));
      displs[0]= 0;//第一個進程的起始點為0
      for(int i = 1; i < size; i++){
          displs[i] = displs[i - 1] + recvpoint[i - 1]; //後續進程的起始點為前一個進程結束位置
      } 
      
    }
   
    //id 0 接收個個進程回傳的凸包點的座標
    MPI_Gatherv(hull_point, hull_size, mpi_point, all_hull_point, recvpoint, displs, mpi_point, 0, MPI_COMM_WORLD);

    
    if(rank == 0){
        int sum = 0;
        for(int i = 0; i < size; i++){
            sum += recvpoint[i];
        }
        Point *ans_hull = (Point*)malloc(n * sizeof(Point));
        int ans_count = convex_hull(all_hull_point, sum, ans_hull);
        for(int i = 0; i < ans_count; i++){
            printf("%d ", ans_hull[i].No);
        }

        free(ans_hull);
        
    }
    free(displs);
    free(recvpoint);
    free(all_hull_point);
    free(hull_point);
    free(chunk_point);
    free(point);
    
    
    MPI_Finalize();
    return 0;
}