#include <iostream>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <thread>
#pragma warning(disable:4996)

using namespace std;

int dsize, label, num, layer, tnum;
int* outsize;
float* x, * y;
float cost;
int epoch = 10000;
float alpha = 0.001;

class Matrix {
private:
    int m;
    int n;
    float* mat;
public:
    Matrix() {
        m = 0;
        n = 0;
    }
    Matrix(int m, int n) :m(m), n(n) {
        mat = new float[m * n];
    }
    Matrix(int m, int n, float* mat) :m(m), n(n), mat(mat) {}
    ~Matrix() {
        delete[] mat;
    }

    void sub(Matrix* opM, Matrix* ret) {
        float* op = opM->getMat();
        float* retf = ret->getMat();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                retf[i * n + j] = this->mat[i * n + j] - op[i * n + j];
            }
        }
    }
    void sub(int op, Matrix* ret) {
        float* retf = ret->getMat();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                retf[i * n + j] = this->mat[i * n + j] - op;
            }
        }
    }
    void rsub(int op, Matrix* ret) {
        float* retf = ret->getMat();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                retf[i * n + j] = op - this->mat[i * n + j];
            }
        }
    }
    void add(Matrix* opM, Matrix* ret) {
        float* op = opM->getMat();
        float* retf = ret->getMat();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                retf[i * n + j] = this->mat[i * n + j] + op[i * n + j];
            }
        }
    }
    void dot(Matrix* opM, Matrix* ret) {
        if (((ret->getM() != this->m) || (ret->getN() != opM->getN())) || (opM->getM() != this->n)) {
            cout << "this cannot be operated1" << endl;
            return;
        }
        int p = opM->getN();
        float* op = opM->getMat();
        float* retf = ret->getMat();
        for (int k = 0; k < p; k++) {
            for (int i = 0; i < m; i++) {
                float sum = 0;
                int in = i * n;
                for (int j = 0; j < n; j++) {
                    sum += mat[in + j] * op[j * p + k];
                }
                retf[i * p + k] = sum;
            }
        }
    }
    void dotT(Matrix* opM, Matrix* ret) {
        if (((ret->getM() != this->m) || (ret->getN() != opM->getM())) || (opM->getN() != this->n)) {
            cout << "this cannot be operated2" << endl;
            return;
        }
        int p = opM->getM();
        float* op = opM->getMat();
        float* retf = ret->getMat();
        for (int k = 0; k < p; k++) {
            for (int i = 0; i < m; i++) {
                float sum = 0;
                int in = i * n;
                int nk = n * k;
                for (int j = 0; j < n; j++) {
                    sum += mat[in + j] * op[nk + j];
                }
                retf[i * p + k] = sum;
            }
        }
    }
    void Tdot(Matrix* opM, Matrix* ret) {
        if (((ret->getM() != this->n) || (ret->getN() != opM->getN())) || (opM->getM() != this->m)) {
            cout << "this cannot be operated3" << endl;
            return;
        }
        int m = n;
        int n = opM->getM();
        int p = opM->getN();
        float* op = opM->getMat();
        float* retf = ret->getMat();
        for (int k = 0; k < p; k++) {
            for (int i = 0; i < m; i++) {
                float sum = 0;
                for (int j = 0; j < n; j++) {
                    sum += mat[j * m + i] * op[j * p + k];
                }
                retf[i * p + k] = sum;
            }
        }
    }
    void multi(Matrix* opM, Matrix* ret) {
        float* op = opM->getMat();
        float* retf = ret->getMat();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++)
                retf[i * n + j] = mat[i * n + j] * op[i * n + j];
        }
    }
    void Tmulti(Matrix* opM, Matrix* ret) {
        int m = opM->getM();
        int n = opM->getN();
        float* op = opM->getMat();
        float* retf = ret->getMat();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++)
                retf[i * n + j] = mat[j * m + i] * op[i * n + j];
        }
    }
    void multiT(Matrix* opM, Matrix* ret) {
        float* op = opM->getMat();
        float* retf = ret->getMat();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++)
                retf[i * n + j] = mat[i * n + j] * op[j * m + i];
        }
    }
    void smulti(float op) {
        for (int i = 0; i < n * m; i++) {
            mat[i] *= op;
        }
    }
    void sdivide(float op) {
        for (int i = 0; i < n * m; i++) {
            mat[i] /= op;
        }
    }
    void copy(Matrix* opM) {
        float* op = opM->getMat();
        for (int i = 0; i < n * m; i++) {
            mat[i] = op[i];
        }
    }

    void print() {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++)
                printf("%f ", mat[i * n + j]);
            cout << endl;
        }
    }

    float* getMat() {
        return mat;
    }
    int getM() {
        return m;
    }
    int getN() {
        return n;
    }
    void setM(int m) {
        this->m = m;
    }
    void setN(int n) {
        this->n = n;
    }
    void setMat(float* mat) {
        //delete[] this->mat;
        this->mat = mat;
    }
};

class Layer {
private:
    Matrix* w;
    Matrix* a;

public:
    Layer() {}
    Layer(int order) {
        w = new Matrix(outsize[order - 1], outsize[order]);
        a = new Matrix(1, outsize[order]);

        srand((unsigned int)time(NULL));

        float* wMat = w->getMat();
        for (int i = 0; i < outsize[order - 1]; i++)
            for (int j = 0; j < outsize[order]; j++)
                wMat[i * outsize[order] + j] = (float)(rand() % 100 - 50) / 1000;
    }
    Layer(int order, Matrix* w) {
        this->w = w;
        a = new Matrix(1, outsize[order]);
    }
    ~Layer() {
        delete a;
    }

    void calcA(Matrix* x) {
        x->dot(w, a);
        float* aMat = a->getMat();
        for (int i = 0; i < a->getN(); i++)
            aMat[i] = 1 / (1 + exp(-aMat[i]));
    }
    Matrix* getW() {
        return w;
    }
    Matrix* getA() {
        return a;
    }
};

Layer** l;
Matrix** w_adj;
Matrix*** t_adj;
Matrix** input;
Matrix** output;

float getLoss(Matrix* out, Matrix* y) {
    float loss = 0;
    for (int i = 0; i < label; i++) {
        loss += pow(out->getMat()[i] - y->getMat()[i], 2);
    }
    loss /= label;
    return loss;
}

void forward(Matrix* x, Layer* l2[]) {
    l2[0]->calcA(x);
    for (int i = 1; i < layer; i++) {
        l2[i]->calcA(l2[i - 1]->getA());
    }
}

void backprop(Matrix* x, Matrix* y, Matrix* d[], Matrix* adj[], Matrix* tmpW[], Layer* l2[], Matrix* tmp1[], Matrix* tmp2[]) {
    l2[layer - 1]->getA()->sub(y, d[layer - 1]);
    for (int i = layer - 2; i >= 0; i--) {
        l2[i + 1]->getW()->dotT(d[i + 1], tmp1[i]);
        l2[i]->getA()->rsub(1, tmp2[i]);
        tmp1[i]->multiT(tmp2[i], d[i]);
    }
    x->Tdot(d[0], tmpW[0]);
    for (int i = 1; i < layer; i++) {
        l2[i - 1]->getA()->Tdot(d[i], tmpW[i]);
        tmpW[i-1]->smulti(alpha);
        adj[i-1]->add(tmpW[i-1], adj[i-1]);
    }
    adj[layer-1]->add(tmpW[layer-1], adj[layer-1]);
}

void ttrain(Matrix* adj[], int t, Layer** l2, Matrix** d, Matrix** tmpW, Matrix** tmp1, Matrix** tmp2) {
    for(int i = 0; i < layer; i++){
        adj[i]->smulti(0);
    }
    Matrix* out = l2[layer - 1]->getA();

    int range = num / tnum;

    for (int i = t * range; i < (t + 1) * range; i++) {
        input[t]->setMat(x + i * dsize);
        output[t]->setMat(y + i * label);

        forward(input[t], l2);
        cost += getLoss(out, output[t]);
        backprop(input[t], output[t], d, adj, tmpW, l2, tmp1, tmp2);
    }
}

void train(thread** tid) {
    l = new Layer *[layer];
    w_adj = new Matrix *[layer];
    t_adj = new Matrix **[tnum];
    Layer*** t_l2 = new Layer **[tnum];
    Matrix*** t_d = new Matrix **[tnum];
    Matrix*** t_tmpW = new Matrix **[tnum];

    input = new Matrix *[tnum];
    output = new Matrix *[tnum];

    Matrix*** t_tmp1 = new Matrix**[tnum]; 
    Matrix*** t_tmp2 = new Matrix**[tnum];


    for (int i = 0; i < layer; i++) {
        l[i] = new Layer(i + 1);
        w_adj[i] = new Matrix(l[i]->getW()->getM(), l[i]->getW()->getN());
    }
    for (int t = 0; t < tnum; t++) {
        t_adj[t] = new Matrix * [layer];
        t_l2[t] = new Layer * [layer];
        t_d[t] = new Matrix * [layer];
        t_tmpW[t] = new Matrix * [layer];
        input[t] = new Matrix(1, dsize, x);
        output[t] = new Matrix(1, label, y);

        t_tmp1[t] = new Matrix*[layer-1];
        t_tmp2[t] = new Matrix*[layer-1];

        for(int i = layer - 2; i >= 0; i--){
            t_tmp1[t][i] = new Matrix(l[i]->getA()->getN(), 1);
            t_tmp2[t][i] = new Matrix(1, l[i]->getA()->getN());
        }

        for (int i = 0; i < layer; i++) {
            t_adj[t][i] = new Matrix(w_adj[i]->getM(), w_adj[i]->getN());
            t_l2[t][i] = new Layer(i + 1, l[i]->getW());
            t_d[t][i] = new Matrix(l[i]->getA()->getM(), l[i]->getA()->getN());
            t_tmpW[t][i] = new Matrix(l[i]->getW()->getM(), l[i]->getW()->getN()); 
        }
    }

    for (int i = 0; i < epoch; i++) {
        cost = 0;
        for (int t = 0; t < tnum; t++) {
            tid[t] = new thread(ttrain, t_adj[t], t, t_l2[t], t_d[t], t_tmpW[t], t_tmp1[t], t_tmp2[t]);
        }

        for (int t = 0; t < tnum; t++) {
            tid[t]->join();
            tid[t]->~thread();
        }

        for (int j = 0; j < layer; j++) {
            w_adj[j]->smulti(0);
            for (int t = 0; t < tnum; t++) {
                w_adj[j]->add(t_adj[t][j], w_adj[j]);
            }
            w_adj[j]->sdivide(num);
            l[j]->getW()->sub(w_adj[j], l[j]->getW());
        }
        if((i+1)%100 == 0) 
            printf("%5d'th iteration | cost: %f | acc: %f\n", i + 1, cost / num, (1-(cost/num))*100);
    }

    delete l[0]->getW();
    for(int t = 0; t < tnum; t++){
        for(int i = layer - 2; i >= 0; i--){
            delete t_tmp1[t][i];
            delete t_tmp2[t][i];
        }
        for(int i = 0; i < layer; i++){
            delete t_adj[t][i];
            delete t_l2[t][i];
            delete t_d[t][i];
            delete t_tmpW[t][i];
        }
        delete[] t_adj[t];
        delete[] t_l2[t];
        delete[] t_d[t];
        delete[] t_tmpW[t];
        delete[] t_tmp1[t];
        delete[] t_tmp2[t];
    }
    for(int i = 0; i < layer; i++){
        delete l[i];
        delete w_adj[i];
    }
    delete[] l;
    delete[] w_adj;
    delete[] t_adj;
    delete[] t_l2;
    delete[] t_d;
    delete[] t_tmpW;
    delete[] t_tmp1;
    delete[] t_tmp2;
}

int main() {
    thread** tid;
    char data;
    int xcnt = 0, ycnt = 0;
    FILE* f = fopen("format.txt", "r");
    if (f == NULL) cout << "Error in file open" << endl;
    fscanf(f, "%d\n%d\n%d\n%d\n", &dsize, &label, &num, &layer);

    outsize = new int[layer + 1];

    x = new float[dsize * num];
    y = new float[label * num];

    outsize[0] = dsize;
    for (int i = 1; i <= layer; i++) {
        fscanf(f, "%d", outsize + i);
    }
    fscanf(f, "%d\n", &tnum);
    fclose(f);

    f = fopen("data.txt", "r");

    while (feof(f) == 0) {
        fread(&data, 1, 1, f);
        if (data == '[') {
            for (int i = 0; i < dsize;) {
                fread(&data, 1, 1, f);
                if (data == '0' || data == '1') {
                    x[xcnt++] = (float)((int)(data - '0'));
                    i++;
                }
            }

            while (data != '[') {
                fread(&data, 1, 1, f);
            }

            for (int i = 0; i < label;) {
                fread(&data, 1, 1, f);
                if (data == '0' || data == '1') {
                    y[ycnt++] = (float)((int)(data - '0'));
                    i++;
                }
            }
        }
    }
    fclose(f);

    tid = new thread * [tnum];

    epoch = 10000;
    alpha = 0.001;
    train(tid);

    return 0;
}