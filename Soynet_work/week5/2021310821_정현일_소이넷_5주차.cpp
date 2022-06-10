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
        this->w = new Matrix(outsize[order - 1], outsize[order]);
        this->w->copy(w);
        a = new Matrix(1, outsize[order]);
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

Layer** l;
Matrix** w_adj;
Matrix*** t_adj;

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

void backprop(Matrix* x, Matrix* y, Matrix* d[], Matrix* adj[], Matrix* tmpW[], Layer* l2[], float alpha) {
    Matrix* tmp1, * tmp2;
    l2[layer - 1]->getA()->sub(y, d[layer - 1]);
    for (int i = layer - 2; i >= 0; i--) {
        tmp1 = new Matrix(l2[i]->getA()->getN(), 1);
        tmp2 = new Matrix(1, l2[i]->getA()->getN());
        l2[i + 1]->getW()->dotT(d[i + 1], tmp1);
        l2[i]->getA()->rsub(1, tmp2);
        tmp1->multiT(tmp2, d[i]);
        delete tmp1;
        delete tmp2;
    }
    x->Tdot(d[0], tmpW[0]);
    for (int i = 1; i < layer; i++) {
        l2[i - 1]->getA()->Tdot(d[i], tmpW[i]);
        tmpW[i-1]->smulti(alpha);
        adj[i-1]->add(tmpW[i-1], adj[i-1]);
    }
    adj[layer-1]->add(tmpW[layer-1], adj[layer-1]);
}

void ttrain(Matrix* adj[], int t, float alpha) {
    Layer** l2 = new Layer * [layer];
    Matrix** d = new Matrix * [layer];
    Matrix** tmpW = new Matrix * [layer];
    for (int i = 0; i < layer; i++) {
        adj[i]->smulti(0);
        l2[i] = new Layer(i + 1, l[i]->getW());
        d[i] = new Matrix(l2[i]->getA()->getM(), l2[i]->getA()->getN());
        tmpW[i] = new Matrix(l2[i]->getW()->getM(), l2[i]->getW()->getN());
    }

    Matrix* input = new Matrix(1, dsize, x);
    Matrix* output = new Matrix(1, label, y);

    Matrix* out = l2[layer - 1]->getA();

    int range = num / tnum;

    for (int i = t * range; i < (t + 1) * range; i++) {
        input->setMat(x + i * dsize);
        output->setMat(y + i * label);

        forward(input, l2);
        cost += getLoss(out, output);
        backprop(input, output, d, adj, tmpW, l2, alpha);
    }

    for (int i = 0; i < layer; i++) {
        delete l2[i];
        delete d[i];
        delete tmpW[i];
    }
    delete[] l2;
    delete[] d;
    delete[] tmpW;
}

void train(thread** tid, float alpha, int epoch) {
    l = new Layer * [layer];
    w_adj = new Matrix * [layer];
    t_adj = new Matrix * *[tnum];

    for (int i = 0; i < layer; i++) {
        l[i] = new Layer(i + 1);
        w_adj[i] = new Matrix(l[i]->getW()->getM(), l[i]->getW()->getN());
    }

    for (int t = 0; t < tnum; t++) {
        t_adj[t] = new Matrix * [layer];
        for (int i = 0; i < layer; i++) {
            t_adj[t][i] = new Matrix(w_adj[i]->getM(), w_adj[i]->getN());
        }
    }

    for (int i = 0; i < epoch; i++) {
        cost = 0;
        for (int t = 0; t < tnum; t++) {
            tid[t] = new thread(ttrain, t_adj[t], t, alpha);
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

    train(tid, 0.001, 1000);

    return 0;
}