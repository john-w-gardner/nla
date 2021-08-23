// header file for some basic matrix operations
// functions found in matops.c and qrGS1.c
void matmat(double A[], int am, int an,
            double B[], int bm, int bn, 
            double (*C)[], int cm, int cn);
void matvec(double A[], int am, int an,
            double x[], int xn,
            double (*b)[]);
void printmat(double A[], int m, int n);
double dotprod(double x[], double y[], int n);
double eucnorm(double x[], int n);
double eucmetric(double x[], double y[], int n);
void addVec(double (*u)[], double v[], int n);
void scaleVec(double (*x)[], double a, int n);
void updateCol(double (*A)[], int m, int n, int i, double v[]);
void transpose1(double A[], int m, int n, double (*At)[]);
void transpose2(double (*A)[], int m, int n, double (*At)[]);
double extractEntry(double A[], int m, int n, int i, int j);
void extractColumn(double A[], int m, int n, int j, double (*v)[]);
void updateEntry(double x, double (*A)[], int m, int n, int i, int j);
double printMax(double v[], int n);
void qrGS1(double A[], int m, int n, double (*Q)[], double (*R)[]);
void qrGS2(double A[], int m, int n, double (*Q)[], double (*R)[]);
void backsub(double A[], int n, double b[], double (*x)[]);
void qrSolve(double A[], int n, double (*x)[], double b[]);
double testSolve(double A[], int n, double x[], double b[]);
