#include <iostream>
#include <math.h>

using namespace std;

float x[] ={0, 0, 1, 1, 0, 0,
          0, 1, 0, 0, 1, 0,
          1, 1, 1, 1, 1, 1,
          1, 0, 0, 0, 0, 1,
          1, 0, 0, 0, 0, 1,
          0, 1, 1, 1, 1, 0,
          0, 1, 0, 0, 1, 0,
          0, 1, 1, 1, 1, 0,
          0, 1, 0, 0, 1, 0,
          0, 1, 1, 1, 1, 0,
          0, 1, 1, 1, 1, 0,
          0, 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0, 0,
          0, 1, 1, 1, 1, 0};
float y[] ={1, 0, 0,
            0, 1, 0,
            0, 0, 1};

class matrix{
private:
    int m;
    int n;
    float *mat;
public:
    matrix(int m, int n):m(m),n(n){
        mat = new float[m*n];
    }
    // 새로운 allocation이 일어나지 않도록
    matrix transpose(){
        matrix ret(n, m);
        float* retf = ret.getMat();
        for(int i = 0; i < n; i++){
            for(int j = 0; j < m; j++){
                retf[i*m+j] = mat[j*n+i];
            }
        }
        return ret;
    }
    matrix dot(matrix opM){
        int n = opM.getM();
        int p = opM.getN();
        float *op = opM.getMat();
        matrix ret(n, p);
        float* retf = ret.getMat();
        for(int k = 0; k < p; k++){
            for(int i = 0; i < m; i++){
                int sum = 0;
                for(int j = 0; j < n; j++){
                    sum += mat[i*n+j] * op[k+j*p];
                }
                retf[i*p + k] = sum;
            }
        }
        return ret;
    }
    matrix multi(matrix opM){
        float *op = opM.getMat();
        matrix ret(m, n);
        float* retf = ret.getMat();
        for(int i = 0; i < m; i++){
            for(int j = 0; j < n; j++)
                retf[i*n + j] = mat[i*n + j] * op[i*n + j];
        }
        
        return ret;
    }

    void print(){
        for(int i = 0; i < m; i++){
            for(int j = 0; j < n; j++)
                cout << mat[i*n + j] << " ";
            cout << endl;
        }
    }

    float* getMat(){
        return mat;
    }
    int getM(){
        return m;
    }
    int getN(){
        return n;
    }
};

class Layer1{
private:
    float *w1;
    float *z1;
    float *a1;

public:
    Layer1(){
        w1 = new float[30*5];
        z1 = new float[5];
        a1 = new float[5];
        for(int i = 0; i < 30; i++)
            for(int j = 0; j < 5; j++)
                w1[i*5+j] = 0.5;
    }
    ~Layer1(){
        delete[] w1;
        delete[] z1;
        delete[] a1;
    }
    void calcZ1(float *x){
        int sum;
        for(int j = 0; j < 5; j++){
            sum = 0;
            for(int i = 0; i < 30; i++){
                    sum += x[i] * w1[i*5];
            }
            z1[j] = sum;
        }
    }
    void calcA1(void){
        for(int i = 0; i < 5; i++)
            a1[i] = 1/(1 + exp(-z1[i]));
    }
    float* getA1(){
        return a1;
    }
};

class Layer2{
private:
    float *w2;
    float *z2;
    float *a2;
public:
    Layer2(){
        w2 = new float[5*3];
        z2 = new float[3];
        a2 = new float[3];
        for(int i = 0; i < 5; i++)
            for(int j = 0; j < 3; j++)
                w2[i*3+j] = 0.5;
    }
    ~Layer2(){
        delete[] w2;
        delete[] z2;
        delete[] a2; 
    }
    void calcZ2(float *x){
        int sum;
        for(int j = 0; j < 3; j++){
            sum = 0;
            for(int i = 0; i < 5; i++){
                    sum += x[i] * w2[i*3];
            }
            z2[j] = sum;
        }
    }
    void calcA2(void){
        for(int i = 0; i < 3; i++)
            a2[i] = 1/(1 + exp(-z2[i]));
    }
    float* getA2(){
        return a2;
    }
};

float getLoss(float *out, float *y){
    float loss=0;
    for(int i = 0; i < 3; i++){
        loss += pow(out[i] - y[i], 2);
    }
    loss /= 3;
    return loss;
}

void forward(float *x, Layer1 *l1, Layer2 *l2){
    // hidden
    for(int i = 0; i < 3; i++){
        l1[i].calcZ1(x+i*30);
        l1[i].calcA1();
    }
    // output layer
    for(int i = 0; i < 3; i++){
        l2[i].calcZ2(l1[i].getA1());
        l2[i].calcA2();
    }
}

void backprop(float *x, Layer1 *l1, Layer2 *l2){
    matrix d1(5, 3);
    matrix d2(3, 3);

    // hidden
    for(int i = 0; i < 3; i++){
        l1[i].calcZ1(x+i*30);
        l1[i].calcA1();
    }
    // output layer
    for(int i = 0; i < 3; i++){
        l2[i].calcZ2(l1[i].getA1());
        l2[i].calcA2();
    }

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            d2.getMat()[i*3+j] = l2[i].getA2()[j] - y[i*3+j];
        }
    }
}

int main(){
    Layer1 l1[3];
    Layer2 l2[3];

    matrix m1(2, 4);
    matrix m2(4, 2);
    for(int i = 0; i < 8; i++){
        m1.getMat()[i] = i+1;
        if(i % 2 == 0)
            m2.getMat()[i] = 1;
        else m2.getMat()[i] = 2;
    }
    matrix m3 = m1.dot(m2);
    m3.print();

    return 0;
}