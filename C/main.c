#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

const int C = 299792458; // 光速 m/s
// const int C = 3e8; // 光速 m/s
const double a = 20 * 1e-3; // 长 mm
const double b = 10 * 1e-3; // 宽 mm
const double h = 1 * 1e-3; // 间隔 mm
const double tolerance = 1e-5; // 允许误差
const double maxIter = 1e5;
const double ezInit = 1; // Ez内点初始
const double kcInit = 0.5; // kc初始

typedef struct {
    int rows, cols;
    double** matrix;
} Matrix;

Matrix* newEz(int rows, int cols, double init);
void updateEz(Matrix* Ez, double kc);
double calcKc(Matrix* Ez);
int calcTM11(Matrix* Ez, double kcInit);
void theory(double a, double b, double* kc, double* lambda, double* freq);
void compute(double kc, double* lambda, double* freq);
void printMatrix(Matrix* Ez);
void matrixToCsv(Matrix* m, FILE* file);
void writeCsv(Matrix* m, const char* name);

void main() {
    Matrix* Ez = newEz((int)(a / h) + 1, (int)(b / h) + 1, ezInit);
    int n = calcTM11(Ez, kcInit);
    printf("迭代次数: %d\n", n);

    printMatrix(Ez);
    writeCsv(Ez, "Ez");

    double kcTheory, lambdaTheory, freqTheory;
    theory(a, b, &kcTheory, &lambdaTheory, &freqTheory);
    double kc = calcKc(Ez);
    double lambda, freq;
    compute(kc, &lambda, &freq);

    printf("理论 kc: %.6f\n", kcTheory);
    printf("计算 kc: %.6f\n", kc);
    printf("理论波长: %.6f m\n", lambdaTheory);
    printf("计算波长: %.6f m\n", lambda);
    printf("理论频率: %.6f Hz\n", freqTheory);
    printf("计算频率: %.6f Hz\n", freq);
}

Matrix* newEz(int rows, int cols, double init) {
    Matrix* m = (Matrix*)malloc(sizeof(Matrix));
    m->rows = rows;
    m->cols = cols;

    m->matrix = (double**)malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        m->matrix[i] = (double*)malloc(cols * sizeof(double));
        for (int j = 0; j < cols; j++) {
            if (i > 0 && i < rows - 1 && j > 0 && j < cols - 1) {
                m->matrix[i][j] = init;
            }
            else {
                m->matrix[i][j] = 0;
            }
        }
    }

    return m;
}

void updateEz(Matrix* Ez, double kc) {
    for (int i = 1; i < Ez->rows - 1; i++) {
        for (int j = 1; j < Ez->cols - 1; j++) {
            double E2 = Ez->matrix[i + 1][j];
            double E3 = Ez->matrix[i - 1][j];
            double E4 = Ez->matrix[i][j + 1];
            double E5 = Ez->matrix[i][j - 1];
            double E1 = (E2 + E3 + E4 + E5) / (4 - kc * kc * h * h);
            Ez->matrix[i][j] = E1;
        }
    }
}

double calcKc(Matrix* Ez) {
    double kckcSum = 0;

    for (int i = 1; i < Ez->rows - 1; i++) {
        for (int j = 1; j < Ez->cols - 1; j++) {
            double E2 = Ez->matrix[i + 1][j];
            double E3 = Ez->matrix[i - 1][j];
            double E4 = Ez->matrix[i][j + 1];
            double E5 = Ez->matrix[i][j - 1];
            double E1 = Ez->matrix[i][j];
            kckcSum += -(E2 + E3 + E4 + E5 - 4 * E1) / (E1 * h * h);
        }
    }

    return sqrt(kckcSum / ((Ez->rows - 2) * (Ez->cols - 2)));
}

int calcTM11(Matrix* Ez, double kcInit) {
    double kc1 = kcInit;
    double kcOld = kc1 + 1;
    Matrix* EzOld = newEz(Ez->rows, Ez->cols, 0);

    int n = 0;
    while (1) {
        n++;
        if (n >= maxIter) {
            break;
        }

        for (int i = 0; i < Ez->rows; i++) {
            for (int j = 0; j < Ez->cols; j++) {
                EzOld->matrix[i][j] = Ez->matrix[i][j];
            }
        }
        updateEz(Ez, kc1);

        if (n % 10 == 0) {
            kcOld = kc1;
            kc1 = calcKc(Ez);
        }

        double EzDiffNorm = 0;
        for (int i = 0; i < Ez->rows; i++) {
            for (int j = 0; j < Ez->cols; j++) {
                double diff = Ez->matrix[i][j] - EzOld->matrix[i][j];
                EzDiffNorm += diff * diff;
            }
        }
        EzDiffNorm = sqrt(EzDiffNorm);
        double kcDiff = fabs(kc1 - kcOld);
        if (EzDiffNorm < tolerance || kcDiff < tolerance) {
            break;
        }
    }

    return n;
}

double calcMatrixDistance(Matrix* m1, Matrix* m2) {
    if (m1->rows != m2->rows || m1->cols != m2->cols) {
        printf("矩阵长宽不匹配\n");
        exit(1);
    }

    double sum = 0;
    for (int i = 0; i < m1->rows; i++) {
        for (int j = 0; j < m1->cols; j++) {
            double diff = m1->matrix[i][j] - m2->matrix[i][j];
            sum += diff * diff;
        }
    }

    return sqrt(sum);
}

void theory(double a, double b, double* kc, double* lambda, double* freq) {
    *kc = sqrt(pow(M_PI / a, 2) + pow(M_PI / b, 2));
    *lambda = 2 * M_PI / *kc;
    *freq = C / *lambda;
}

void compute(double kc, double* lambda, double* freq) {
    *lambda = 2 * M_PI / kc;
    *freq = C / *lambda;
}

void printMatrix(Matrix* Ez) {
    printf("Ez:\n");
    for (int i = 0; i < Ez->rows; i++) {
        for (int j = 0; j < Ez->cols; j++) {
            printf("%.6f ", Ez->matrix[i][j]);
        }
        printf("\n");
    }
}

void matrixToCsv(Matrix* m, FILE* file) {
    for (int i = 0; i < m->rows; i++) {
        for (int j = 0; j < m->cols; j++) {
            if (j > 0) {
                fprintf(file, ",");
            }
            double num = m->matrix[i][j];
            if (num == 0) {
                fprintf(file, "0");
            }
            else {
                fprintf(file, "%.18f", num);
            }
        }
        fprintf(file, "\n");
    }
}

void writeCsv(Matrix* m, const char* name) {
    char filePath[256];
    sprintf(filePath, "%s.csv", name);

    remove(filePath);
    FILE* csvFile = fopen(filePath, "w");
    if (csvFile == NULL) {
        printf("无法打开文件\n");
        exit(1);
    }
    matrixToCsv(m, csvFile);
    fclose(csvFile);
}
