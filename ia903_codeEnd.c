#include <stdio.h>
#include <mpi.h>  // This is the header file for MPI functions
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char* argv[]){

        int myId, numProc;
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &numProc);
        MPI_Comm_rank(MPI_COMM_WORLD, &myId);
        MPI_Status status;
        FILE* fp;
        char filename[100] = "";

        clock_t start_p1, start_p3, end_p1, end_p3;
        double local_start, local_finish, local_elapsed, elapsed;
        start_p1 = clock();
        if(argc !=5){
                printf("usage:  ./checkdiv N a b c\n");
                printf("N: the upper bound of the range [2,N]\n");
                printf("a: first divisor\n");
                printf("b: second divisor\n");
                printf("c: third divisor\n");
                exit(1);
        }

        int last_count = 0;
        int N = (unsigned int) atoi(argv[1]);
        int a = (unsigned int) atoi(argv[2]);
        int b = (unsigned int) atoi(argv[3]);
        int c = (unsigned int) atoi(argv[4]);
        int* total_values = (int*) malloc(sizeof(int)*N);

        end_p1 = clock();

        local_start = MPI_Wtime();

	if(myId==0){
		int size = N/numProc;
                if(numProc == 1)
                        size = (N/numProc) + 1;
                for(int i = 2; i < size; i++){
                        if (i % a == 0) {
                                total_values[last_count] = i;
                                last_count++;
                        } else if (i % b == 0) {
                                total_values[last_count] = i;
                                last_count++;
                        } else if (i % c == 0) {
                                total_values[last_count] = i;
                                last_count++;
                        }
                }
                int *all = 0;

                for(int i = 1; i < numProc; i++){
                        int count;
                        MPI_Recv(&count, 1, MPI_INT, i, 0, MPI_COMM_WORLD,
                                        &status);
                        all = realloc(all, sizeof(int)*count);

                        MPI_Recv(all, count, MPI_INT, i , 0,
                                        MPI_COMM_WORLD,&status);

                        for(int j = 0; j < count; j++){
                                total_values[last_count] = all[j];
                                last_count++;
                        }
                }
                free(all);
        }
        if(myId!=numProc-1){
                int partial_N = (N/numProc)*myId;
                int* all_values =  malloc(sizeof(int)*N/numProc);
                int count = 0;
                for(int i = partial_N; i < (partial_N+(N/numProc)); i++){
                        if (i % a == 0) {
                                all_values[count] = i;
                                count++;
                        } else if (i % b == 0) {
                                all_values[count] = i;
                                count++;
                        } else if (i % c == 0) {
                                all_values[count] = i;
                                count++;
                        }
                }
              MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
              MPI_Send(all_values, count, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        else{
                int partial_N = ((N/numProc)*myId);
                int count = 0;
                int* all_values=malloc(sizeof(int)*((N/numProc) +(N%numProc)));
                for(int i = partial_N; i <= N; i++){
                        if (i % a == 0) {
                                all_values[count] = i;
                                count++;
                        } else if (i % b == 0) {
                                all_values[count] = i;
                                count++;
                        } else if (i % c == 0) {
                                all_values[count] = i;
                                count++;
                        }
                }
                MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                MPI_Send(all_values, count, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        local_finish = MPI_Wtime();
        local_elapsed = local_finish - local_start;
        MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX,
                           0, MPI_COMM_WORLD);

        start_p3 = clock();
        strcpy(filename, argv[1]);
        strcat(filename, ".txt");

        if(! (fp  = fopen(filename, "w + t"))){
                 printf("Cannot create file %s\n", filename);
                 exit(1);
        }
        int prev_Num = 0;
        for(int i = 0; i < last_count; i++){
               if(total_values[i] && prev_Num!=total_values[i]){
                        fprintf(fp, "%d\n", total_values[i]);
                        prev_Num = total_values[i];
                }
        }
        fclose(fp);
        end_p3 = clock();

        MPI_Finalize();
        if(myId==0){
                double total_t1 = (double)(end_p1 - start_p1)/CLOCKS_PER_SEC;
                double total_t2 = (double)(end_p3 - start_p3)/CLOCKS_PER_SEC;
                printf("Part 1 proccesed for %1.2f seconds\n", total_t1);
                printf("Part 2 processed for %1.2f seconds\n", elapsed);
                printf("Part 3 processed for %1.2f seconds\n", total_t2);
        }
        return 0;
}
