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

    // Matrix* transpose(){
    //     Matrix* ret = new Matrix(n, m);
    //     float* retf = ret->getMat();
    //     for(int i = 0; i < n; i++){
    //         for(int j = 0; j < m; j++){
    //             retf[i*m+j] = mat[j*n+i];
    //         }
    //     }
    //     return ret;
    // }
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
            this->mat[i] *= op;
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
                wMat[i * outsize[order] + j] = (float)(rand() % 100) / 1000;
    }
    ~Layer() {
        delete w;
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

float getLoss(Matrix* out, Matrix* y) {
    float loss = 0;
    for (int i = 0; i < label; i++) {
        loss += pow(out->getMat()[i] - y->getMat()[i], 2);
    }
    loss /= label;
    return loss;
}

void forward(Matrix* x, Layer* l[]) {
    l[0]->calcA(x);
    for (int i = 1; i < layer; i++) {
        l[i]->calcA(l[i - 1]->getA());
    }
}

Layer*** LL;

void backprop(Matrix* x, Matrix* y, Matrix* d[], Matrix* w_adj[], Layer* l[], float alpha) {
    Matrix* tmp1, * tmp2;
    l[layer - 1]->getA()->sub(y, d[layer - 1]);
    for (int i = layer - 2; i >= 0; i--) {
        tmp1 = new Matrix(l[i]->getA()->getN(), 1);
        tmp2 = new Matrix(1, l[i]->getA()->getN());
        l[i + 1]->getW()->dotT(d[i + 1], tmp1);
        l[i]->getA()->rsub(1, tmp2);
        tmp1->multiT(tmp2, d[i]);
        delete tmp1;
        delete tmp2;
    }
    x->Tdot(d[0], w_adj[0]);
    for (int i = 1; i < layer; i++) {
        l[i - 1]->getA()->Tdot(d[i], w_adj[i]);
        w_adj[i - 1]->smulti(alpha);
        l[i - 1]->getW()->sub(w_adj[i - 1], l[i - 1]->getW());
    }
    l[layer - 1]->getW()->sub(w_adj[layer - 1], l[layer - 1]->getW());
}

void train(float* x, float* y, Layer* l[], float alpha, int epoch, int t) {
    Matrix** d = new Matrix * [layer];
    Matrix** w_adj = new Matrix * [layer];
    for (int i = 0; i < layer; i++) {
        d[i] = new Matrix(l[i]->getA()->getM(), l[i]->getA()->getN());
        w_adj[i] = new Matrix(l[i]->getW()->getM(), l[i]->getW()->getN());
    }

    Matrix* input = new Matrix(1, dsize, x);
    Matrix* output = new Matrix(1, label, y);

    Matrix* out = l[layer - 1]->getA();

    int range = num / tnum;

    for (int k = 0; k < epoch; k++) {
        float cost = 0;
        for (int i = t * range; i < (t + 1) * range; i++) {
            input->setMat(x + i * dsize);
            output->setMat(y + i * label);

            forward(input, l);
            cost += getLoss(out, output);
            backprop(input, output, d, w_adj, l, alpha);
        }
        printf("%d'th thread\n%d'th iteration | cost: %f\n", t + 1, k + 1, cost / num);
        //cout << t + 1 << "'th thread == " << k + 1 << "th iteration" << " | " << "cost: " << cost / num << endl;
    }
}

void ttrain(int t) {
    Layer** l = new Layer * [layer];
    for (int i = 0; i < layer; i++) {
        l[i] = new Layer(i + 1);
    }
    LL[t] = l;
    train(x, y, l, 0.001, 100, t);
}

int main() {
    char data;
    int xcnt = 0, ycnt = 0;
    thread** tid;
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

    LL = (Layer***) new Layer * *[tnum];

    tid = (thread**)malloc(sizeof(thread*) * tnum);
    for (int t = 0; t < tnum; t++) {
        tid[t] = new thread(ttrain, t);
    }

    for (int t = 0; t < tnum; t++) {
        tid[t]->join();
    }

    /*for (int t = 0; t < tnum; t++) {
        for (int i = 0; i < layer; i++)
            delete LL[t][i];
        delete[] LL[t];
    }
    delete[] LL;*/

    return 0;
}