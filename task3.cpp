#include <iostream>
#include <cmath>
#include <memory>
#include <omp.h>
#include <string>
double epsilon  = 0.00001;
#ifdef FIRSTVAR
#ifdef TYPEDYNAMIC
#define  MYVAR dynamic
#endif
#ifdef TYPESTATIC
#define MYVAR static
#endif
#ifdef TYPEGUIDED
#define MYVAR guided
#endif
int main(int argc, char const *argv[])
{
    int N = std::stoi(argv[1]);
    int thr_cnt = std::stoi(argv[2]);
    double tau = 0.0001;
    omp_set_num_threads(thr_cnt);
    // объявление векторов и матрицы
    std::unique_ptr<std::unique_ptr<double[]>[]> matrix_A(new std::unique_ptr<double[]>[N]);
    std::unique_ptr<double[]> vector_X(new double[N]); 
    std::unique_ptr<double[]> vector_b(new double[N]);
    std::unique_ptr<double[]> vector_c(new double[N]); // хранит Ax-b
    double l2normSum = 0.0;
    // double prevnormsum = 0.0;
    double l2normb = 0.0; // норма l2 для вектора правой части, посчитаем её заранее 
    // инициализация векторов и матрицы
   #pragma omp parallel 
    {
       #pragma omp for schedule(MYVAR, int(N/(thr_cnt))) 
        for(int i=0;i<N;i++){    
            matrix_A[i] = std::make_unique<double[]>(N);
            vector_X[i] = 0.0;
            vector_b[i] = (double)N+1;
            for (int j = 0; j < N; j++){
                    matrix_A[i][j] = (i==j) ? 2.0:1.0;
            }
        }
    }
    double start = omp_get_wtime();
    while (true){
    #pragma omp parallel
    {
        // проверяем, достаточно ли мы приблизили решение
        // умножаем А на вектор Х, потом отнимаем vector_b для хранения результата юзаем vector_c
        #pragma omp for schedule(MYVAR,int(N/(thr_cnt))) 
        for (int i=0 ; i < N ; i++ ){
            vector_c[i] = 0.0;
            for (size_t j = 0; j < N; j++)
            {
                vector_c[i] += matrix_A[i][j] * vector_X[j]; // c = Ax
            }
            vector_c[i] -= vector_b[i]; // c = Ax - b
        }
    }
    #pragma omp parallel
    {
        double sumloc = 0.0;
        double sumlocb = 0.0;
        #pragma omp for schedule(MYVAR,int(N/(thr_cnt))) 
        for (int i = 0;i<N;i++){
            sumloc+=vector_c[i]*vector_c[i];
            sumlocb+=vector_b[i]*vector_b[i];
        }
        #pragma omp atomic
        l2normb+=sumlocb;
        #pragma omp atomic
        l2normSum += sumloc;
    }
    //std::cout <<std::fixed << l2normSum/l2normb <<" " <<vector_X[24] <<std::endl;
    if (l2normSum/l2normb<epsilon) break; 
    l2normSum = 0 ;
    l2normb = 0;
    #pragma omp parallel
    {
        #pragma omp for schedule(MYVAR,int(N/(thr_cnt))) 
        for (int i = 0;i < N; i++)
        {
            vector_X[i] -= tau*vector_c[i]; 
        }
    }
    }
    double end = omp_get_wtime();
    std::cout<<end - start<<std::endl;
    return 0;
}
#endif

#ifdef SECONDVAR
//второй вариант
int main(int argc, char const *argv[])
{
    int N = std::stoi(argv[1]);
    int thr_cnt = std::stoi(argv[2]);
    double tau = 0.0001;

    // объявление векторов и матрицы
    std::unique_ptr<std::unique_ptr<double[]>[]> matrix_A(new std::unique_ptr<double[]>[N]);
    std::unique_ptr<double[]> vector_X(new double[N]); 
    std::unique_ptr<double[]> vector_b(new double[N]);
    std::unique_ptr<double[]> vector_c(new double[N]); // хранит Ax-b
    double l2normSum = 0.0;
    double prevnormsum = 0.0;
    double l2normb = 0.0; // норма l2 для вектора правой части, посчитаем её заранее 
    // инициализация векторов и матрицы 
    #pragma omp parallel num_threads(thr_cnt)
    {
        int nthreads = omp_get_num_threads();
        int threadid = omp_get_thread_num();
        int items_per_thread = N / nthreads;
        int lb = threadid * items_per_thread;
        int ub = (threadid == nthreads - 1) ? (N - 1) : (lb + items_per_thread - 1);
        for (int i=lb ; i < ub ; i++ ){
            matrix_A[i] = std::make_unique<double[]>(N);
            vector_X[i] = 0.0;
            vector_b[i] = (double)N+1;
            for (int j = 0; j < N; j++){
                    matrix_A[i][j] = (i==j) ? 2.0:1.0;
            }
        }
    }
    // начало алгоритма
    bool flag = true;
    double start = omp_get_wtime();
    #pragma omp parallel num_threads(thr_cnt)
    {
    while(flag){
        
        int nthreads = omp_get_num_threads();
        int threadid = omp_get_thread_num();
        int items_per_thread = N / nthreads;
        int lb = threadid * items_per_thread;
        int ub = (threadid == nthreads - 1) ? (N - 1) : (lb + items_per_thread - 1);
    // проверяем, достаточно ли мы приблизили решение
        // умножаем А на вектор Х, потом отнимаем vector_b для хранения результата юзаем vector_c
        for (int i=lb ; i < ub ; i++ ){
            vector_c[i] = 0.0;
            for (size_t j = 0; j < N; j++)
            {
                vector_c[i] += matrix_A[i][j] * vector_X[j]; // c = Ax
            }
            vector_c[i] -= vector_b[i]; // c = Ax - b
        }
        // for (size_t i = lb; i < ub; i++)
        // {
        // }
        
        // найдём л2 норму для произведения и для вектора b , сделам сумму через локальную переменную и атомарную операцию
        double sumloc = 0.0;
        double sumlocb = 0.0;
        for (int i = lb;i<ub;i++){
            sumloc+=vector_c[i]*vector_c[i];
            sumlocb+=vector_b[i]*vector_b[i];
        }
        #pragma omp atomic
        l2normb+=sumlocb;
        #pragma omp atomic
        l2normSum += sumloc;
        // #pragma omp barrier
        // #pragma omp single
        // std::cout <<std::fixed << l2normSum/l2normb <<" " <<vector_X[24] <<std::endl;
        for (int i = lb;i < ub; i++)
        {
            vector_X[i] -= tau*vector_c[i]; 
        }
        
        #pragma omp barrier
        #pragma omp single
        if (l2normSum/l2normb<epsilon) flag=false;
        l2normb = 0.0;
        l2normSum = 0.0;
    }
    }
    double end = omp_get_wtime();
    std::cout << end - start<<std::endl;
    return 0;
}
#endif