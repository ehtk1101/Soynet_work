#include <iostream>
#include <math.h>

using namespace std;

char data;
int size, label, num, layer;
int *outsize;
float *x, *y;

class Matrix{
private:
    int m;
    int n;
    float *mat;
public:
    Matrix(){
        m = 0;
        n = 0;
    }
    Matrix(int m, int n):m(m),n(n){
        mat = new float[m*n];
    }
    ~Matrix(){
        delete[] mat;
    }

    // 새로운 allocation이 일어나지 않도록
    Matrix* transpose(){
        Matrix* ret = new Matrix(n, m);
        float* retf = ret->getMat();
        for(int i = 0; i < n; i++){
            for(int j = 0; j < m; j++){
                retf[i*m+j] = mat[j*n+i];
            }
        }
        return ret;
    }
    void sub(Matrix *opM, Matrix *ret){
        float* op = opM->getMat();
        float* retf = ret->getMat();
        for(int i = 0; i < m; i++){
            for(int j = 0; j < n; j++){
                retf[i*n+j] = this->mat[i*n+j] - op[i*n+j];
            }
        }
    }
    void add(Matrix *opM, Matrix *ret){
        float* op = opM->getMat();
        float* retf = ret->getMat();
        for(int i = 0; i < m; i++){
            for(int j = 0; j < n; j++){
                retf[i*n+j] = this->mat[i*n+j] + op[i*n+j];
            }
        }
    }
    void dot(Matrix *opM, Matrix *ret){
        if(((ret -> getM() != this->m) && (ret->getN() != opM->getN())) && (opM->getM() != this->n)){
            cout << "this cannot be operated" << endl;
            return;
        }
        int n = opM->getM();
        int p = opM->getN();
        float *op = opM->getMat();
        float* retf = ret->getMat();
        for(int k = 0; k < p; k++){
            for(int i = 0; i < m; i++){
                int sum = 0;
                for(int j = 0; j < n; j++){
                    sum += mat[i*n+j] * op[k+j*p];
                }
                retf[i*p + k] = sum;
            }
        }
    }
    void multi(Matrix *opM, Matrix *ret){
        float *op = opM->getMat();
        float* retf = ret->getMat();
        for(int i = 0; i < m; i++){
            for(int j = 0; j < n; j++)
                retf[i*n + j] = mat[i*n + j] * op[i*n + j];
        }
    }
    void smulti(float op){
        for(int i = 0; i < n*m; i++){
            this->mat[i] *= op;
        }
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
    void setM(int m){
        this->m = m;
    }
    void setN(int n){
        this->n = n;
    }
    void setMat(float* mat){
        delete[] this->mat;
        this->mat = mat;
    }
};

class Layer{
private:
    Matrix *w;
    Matrix *z;
    Matrix *a;

public:
    Layer(int order){
        w = new Matrix(outsize[order-1], outsize[order]);
        z = new Matrix(1, outsize[order]);
        a = new Matrix(1, outsize[order]);
        float* wMat = w->getMat();
        for(int i = 0; i < outsize[order-1]; i++)
            for(int j = 0; j < outsize[order]; j++)
                wMat[i*outsize[order]+j] = 0.5;
    }
    ~Layer(){
        delete w;
        delete z;
        delete a;
    }
    void calcZ(Matrix *x){
        x->dot(w, z);
    }
    void calcA(void){
        float* aMat = a->getMat();
        float* zMat = z->getMat();
        for(int i = 0; i < 5; i++)
            aMat[i] = 1/(1 + exp(-zMat[i]));
    }
    Matrix* getW(){
        return w;
    }
    Matrix* getA(){
        return a;
    }
};

float getLoss(Matrix *out, Matrix *y){
    float loss=0;
    for(int i = 0; i < label; i++){
        loss += pow(out->getMat()[i] - y->getMat()[i], 2);
    }
    loss /= label;
    return loss;
}

Matrix* forward(Matrix *x, Layer *l[]){
    l[0]->calcZ(x);
    l[0]->calcA();
    for(int i = 1; i < layer; i++){
        l[i]->calcZ(l[i-1]->getA());
        l[i]->calcA();
    }
    return l[layer-1]->getA();
}

void backprop(Matrix *x, Matrix *y, Layer *l[], float alpha){
    l[0]->calcZ(x);
    l[0]->calcA();
    for(int i = 1; i < layer; i++){
        l[i]->calcZ(l[i-1]->getA());
        l[i]->calcA();
    }

    Matrix *d2 = new Matrix(1, l[1]->getA()->getN());
    Matrix *d1 = new Matrix(1, l[0]->getA()->getN());
    l[layer-1]->getA()->sub(y, d2);
    Matrix *tmp = new Matrix(outsize[1], 1);
    l[1]->getW()->dot(d2->transpose(), tmp);
    Matrix *tmp2 = tmp->transpose();
    delete tmp;
    l[0]->getA()->multi(l[0]->getA(), tmp); // **
    tmp->multi(tmp2, d1);

    Matrix *w1_adj;
    Matrix *w2_adj;
    x->transpose()->dot(d1, w1_adj);
    l[1]->getA()->transpose()->dot(d1, w2_adj);
    w1_adj->smulti(alpha);
    w2_adj->smulti(alpha);
    l[0]->getA()->sub(w1_adj, l[0]->getA());
    l[1]->getA()->sub(w2_adj, l[1]->getA());
}

void train(Matrix *x[], Matrix *y[], Layer *l[], float alpha, int epoch){
    Matrix *out[epoch];
    float *loss = new float[epoch];
    for(int j = 0; j < num; j++){
        for(int i = 0 ; i < epoch; i++){
            out[i] = forward(x[j*epoch+i], l);
            loss[i] = getLoss(out[i], y[i]);
            backprop(x[j*epoch+i], y[j*epoch+i], l, alpha);
        }
        int sum = 0;
        for(int i = 0; i < epoch; i++)
            sum += loss[i];
        cout << "epochs:" << j+1 << "======== acc:" << (1-(sum/epoch))*100 << endl;
    }
}

int main(){
    int xcnt = 0, ycnt = 0;
    FILE *f = fopen("data.txt", "r");
    if(f == NULL) cout << "Error in file open" << endl;
    fscanf(f, "%d\n%d\n%d\n%d\n", &size, &label, &num, &layer);
    outsize = new int[layer];

    x = new float[size*num];
    y = new float[label*num];

    outsize[0] = size;
    for(int i = 1; i < layer; i++){
        fscanf(f, "%d", outsize+i);
    }

    while(feof(f) == 0){
        fread(&data, 1, 1, f);
        if(data == '['){
            for(int i = 0; i < size;){
                fread(&data, 1, 1, f);
                if(data == '0' || data == '1'){
                    x[xcnt++] = (float)((int)(data - '0'));   
                    i++;
                }
            }

            while(data != '['){
                fread(&data, 1, 1, f);
            }

            for(int i = 0; i < label;){
                fread(&data, 1, 1, f);
                if(data == '0' || data == '1'){
                    y[ycnt++] = (float)((int)(data - '0'));   
                    i++;
                }
            }
        }
    }

    

    /*
    Matrix *input[num];
    Matrix *output[num];
    for(int i = 0; i < num; i++){
        input[i] = new Matrix(1, size);
        output[i] = new Matrix(1, label);
        float* tmp1 = input[i] -> getMat();
        float* tmp2 = output[i] -> getMat();

        for(int j = 0; j < size; j++){
            tmp1[j] = x[i*size + j];
        }
        for(int j = 0; j < label; j++){
            tmp2[j] = y[i*label + j];
        }
    }

    Layer *l[layer];
    for(int i = 0; i < layer; i++){
        l[i] = new Layer(i+1);
    }

    train(input, output, l, 0.01, 10);
    */

    return 0;
}