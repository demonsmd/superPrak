#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <ctime>
#include <time.h>
#include <stdlib.h>
#include <mpi.h>

#define Min(A,B) ((A)<(B)?(A):(B))

using namespace std;
vector<vector<float> >uPrev;
vector<vector<float> >U;
vector<float> centers;

void generateInFile(int N){
    ofstream fin("inFile.txt");
    for (int i=0;i<N;i++){
        fin<<(int)(rand()%100000)<<" ";
    }
}

void generateMembershipDegree(const vector<float>* inVec, int C){
    uPrev.resize(inVec->size());
    U.resize(inVec->size());

    centers.resize(C);

    for (int i=0;i<inVec->size();i++){
        U[i].resize(C);
        uPrev[i].resize(C);
    }

    //генерация
    for (int i=0;i<inVec->size();i++){
        float sum = 0;
        for (int j=0;j<C;j++){
            sum+=U[i][j]=uPrev[i][j]=rand()%100;
        }
        for (int j=0;j<C;j++){
            U[i][j]/=sum;
            uPrev[i][j]/=sum;
        }
    }

}


void FCM (const vector<float>* inVec, int C, float m, float eps){
    int ProcNum, rank;
    MPI_Comm_size(MPI_COMM_WORLD,&ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    float global_numerator = 0;
    float global_denumerator = 0;
    float global_max = 0;

    float deside = 0;
    float newDeside = 0;
    generateMembershipDegree(inVec, C);            //заполняет uPrev[][] случайными нормрованными коэффициентами
    //ofstream fout("outFile.txt");

    int itaration = 0;
    while (itaration<50){
        itaration++;
        newDeside = deside;

        //вычисление центров кластеров
        for (int j=0; j<C; j++){           //для каждого кластера
            float numerator = 0;
            float denumerator = 0;
            for (int i=rank * (inVec->size() / ProcNum); i<Min((rank + 1) * (inVec->size()/ProcNum),inVec->size()); i++){        //для каждого вектора
                numerator+=powf(U[i][j], m) * (*inVec)[i];
                denumerator+=powf(U[i][j], m);
                //fout<<"c"<<numerator<<" "<<denumerator<<endl;
            }
            global_numerator = -1;
            global_denumerator = -1;
            MPI_Allreduce( &numerator, &global_numerator, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce( &denumerator, &global_denumerator, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            centers[j] = global_numerator / global_denumerator;
            //fout<<"gc"<<global_numerator<<" "<<global_denumerator<<endl;
        }

        //вычисление коэффициентов u[][]
        for (int i=rank * (inVec->size() / ProcNum); i<Min((rank + 1) * (inVec->size()/ProcNum),inVec->size()); i++)    //для каждого вектора
            for (int j=0; j<C; j++){   //для каждого кластера
                float sum = 0;
                for (int k=0; k<C; k++)
                    sum+= powf(((*inVec)[i] - centers[j])/((*inVec)[i] - centers[k]), 2/(m-1));
                U[i][j] = 1/sum;
            }

        //определения максимального различия между u[][] и uPrev[][]
        float max = 0;
        for (int i=rank * (inVec->size() / ProcNum); i<Min((rank + 1) * (inVec->size()/ProcNum),inVec->size()); i++)    //для каждого вектора
            for (int j=0; j<C; j++)   //для каждого кластера
                if (fabs(uPrev[i][j] - U[i][j]) > max)
                    max = fabs(uPrev[i][j] - U[i][j]);
        MPI_Allreduce( &max, &global_max, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

        if (eps <= global_max)
            for (int i=0; i<inVec->size(); i++)    //для каждого вектора
                for (int j=0; j<C; j++)   //для каждого кластера
                    uPrev[i][j]=U[i][j];
        else
            break;

    }

    
    
    //fout<<"====== results ======="<<endl;
    /*for (int i=rank * (inVec->size() / ProcNum); i<Min((rank + 1) * (inVec->size()/ProcNum),inVec->size()); i++){
        fout<<"vector "<<i<<": ";
        for (int j=0;j<C;j++)
            fout<<U[i][j]<<" , ";
        fout<<endl;
    }
    fout.close();*/

}


int main(int argc, char *argv[])
{
    generateInFile(10000);

    ifstream fin;
    fin.open("inFile.txt");

    if (!fin.is_open()){
        cout<<"cant open File"<<endl;
        exit(-1);
    }

    int i;
    vector<float> inVec;
    while (fin>>i){
        inVec.push_back(i);
    }

    int C = 5;
    clock_t start_time =  clock();

    MPI_Init(&argc,&argv);
    int ProcNum, rank;
    MPI_Comm_size(MPI_COMM_WORLD,&ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    FCM (&inVec, C, 2, 0.00001);

    if (rank==0){
        cout<<"time = "<<(float)(clock() - start_time)/CLOCKS_PER_SEC<<endl;
        ofstream fout("outFile.txt");
        fout<<"time = "<<(float)(clock() - start_time)/CLOCKS_PER_SEC<<endl;
        fout.close();
    }

    cout<<"done!"<<endl;

    MPI_Finalize();

    printf("\nDone\n");

    return 0;
}
