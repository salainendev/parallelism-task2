#include <iostream>
#include <cmath>
#include <memory>
#include <omp.h>
#include <string>
double epsilon  = 0.00001;

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
    if (argc<3){
        std::cout << "должно быть больше аргументов сначала размер матрицы потом колво потоков"<<std::endl;
        return 0;
    }

    int N = std::stoi(argv[1]);
    int thr_cnt = std::stoi(argv[2]);
    double tau = 0.0001;
    omp_set_num_threads(thr_cnt);

    // объявление векторов и матрицы
    std::unique_ptr<double[]> matrix_A(new double[N*N]);
    std::unique_ptr<double[]> vector_X(new double[N]); 
    std::unique_ptr<double[]> vector_b(new double[N]);
    std::unique_ptr<double[]> vector_c(new double[N]); // хранит Ax-b


    double l2normSum = 0.0;
    double prevnormsum = 0.0;
    double l2normb = 0.0; // норма l2 для вектора правой части, посчитаем её заранее 


    // инициализация векторов и матрицы
   #pragma omp parallel 
    {
       #pragma omp for schedule(MYVAR, int(N/(thr_cnt*2))) 
        for(int i=0;i<N;i++){    
            
            vector_X[i] = 0.0;
            vector_b[i] = (double)N+1;
            #pragma omp atomic
            l2normb += vector_b[i]*vector_b[i];
            
            
            for (int j = 0; j < N; j++){
                    matrix_A[i*N+j] = (i==j) ? 2.0:1.0;
            }
        }
    }

    l2normb = sqrt(l2normb); // посчитали знаменатель критерия завершения счёта 

    double start = omp_get_wtime();
    while (true)
        {
        #pragma omp parallel
        {
            // проверяем, достаточно ли мы приблизили решение
            // умножаем А на вектор Х, потом отнимаем vector_b для хранения результата юзаем vector_c
            #pragma omp for schedule(MYVAR,int(N/(thr_cnt*2))) 
            for (int i=0 ; i < N ; i++ ){
                vector_c[i] = 0.0;
                for (size_t j = 0; j < N; j++)
                {
                    vector_c[i] += matrix_A[i*N+j] * vector_X[j]; // c = Ax
                }
                vector_c[i] -= vector_b[i]; // c = Ax - b
            }
        }


        #pragma omp parallel
        {
            double sumloc = 0.0;
            
            #pragma omp for schedule(MYVAR,int(N/(thr_cnt*2))) 
            for (int i = 0;i<N;i++){
                sumloc+=vector_c[i]*vector_c[i];
            
            }

            #pragma omp atomic
            l2normSum += sumloc;
        }
    
        l2normSum = sqrt(l2normSum);

        if (l2normSum/l2normb<epsilon) break;

        if (l2normSum>prevnormsum && prevnormsum!=0){  
            std::cout << "расхождение(( "<<std::endl;
            break;
            }

    
        prevnormsum = l2normSum; 
        l2normSum = 0 ;
        
        #pragma omp parallel
        {
            #pragma omp for schedule(MYVAR,int(N/(thr_cnt*2)))
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


