#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <ctime>
#include <time.h>
#include <stdlib.h>

using namespace std;
vector<vector<float> >uPrev;
vector<vector<float> >U;
vector<float> centers;

void generateInFile(int N){
    ofstream fin("../superPrak/inFile.txt");
    for (int i=0;i<N;i++){
        fin<<(int)(rand()%5000)<<" ";
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
    float deside = 0;
    float newDeside = 0;
    generateMembershipDegree(inVec, C);            //заполняет uPrev[][] случайными нормрованными коэффициентами
    while (true){
        newDeside = deside;

        //вычисление центров кластеров
        for (int j=0; j<C; j++){           //для каждого кластера
            float numerator = 0;
            float denumerator = 0;
            for (int i=0; i<inVec->size(); i++){        //для каждого вектора
                numerator+=powf(U[i][j], m) * (*inVec)[i];
                denumerator+=powf(U[i][j], m);
            }
            centers[j] = numerator / denumerator;
        }

        //вычисление коэффициентов u[][]
        for (int i=0; i<inVec->size(); i++)    //для каждого вектора
            for (int j=0; j<C; j++){   //для каждого кластера
                float sum = 0;
                for (int k=0; k<C; k++)
                    sum+= powf(((*inVec)[i] - centers[j])/((*inVec)[i] - centers[k]), 2/(m-1));
                U[i][j] = 1/sum;
            }

        //определения максимального различия между u[][] и uPrev[][]
        float max = 0;
        for (int i=0; i<inVec->size(); i++)    //для каждого вектора
            for (int j=0; j<C; j++)   //для каждого кластера
                if (fabs(uPrev[i][j] - U[i][j]) > max)
                    max = fabs(uPrev[i][j] - U[i][j]);

        if (eps <= max)
            for (int i=0; i<inVec->size(); i++)    //для каждого вектора
                for (int j=0; j<C; j++)   //для каждого кластера
                    uPrev[i][j]=U[i][j];
        else
            break;

    }
}



int main(int argc, char *argv[])
{
//    generateInFile(10000);

    ifstream fin;
    fin.open("../superPrak/inFile.txt");

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
    FCM (&inVec, C, 2, 0.00001);

    ofstream fout("../superPrak/outFile.txt");
    fout<<"time = "<<(float)(clock() - start_time)/CLOCKS_PER_SEC<<endl;
    cout<<"time = "<<(float)(clock() - start_time)/CLOCKS_PER_SEC<<endl;
    fout<<"====== results ======="<<endl;
    for (int i=0;i<U.size();i++){
        for (int j=0;j<C;j++)
            fout<<U[i][j]<<" , ";
        fout<<endl;
    }
    fout.close();
    cout<<"done!"<<endl;



    return 0;
}
