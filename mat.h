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
void transpose1(double A[], int m, int n, double (*At)[]);
void transpose2(double (*A)[], int m, int n, double (*At)[]);
double extractEntry(double A[], int m, int n, int i, int j);
void extractColumn(double A[], int m, int n, int j, double (*v)[]);
void extractRow(double A[], int m, int n, int j, double (*v)[]);
void updateCol(double (*A)[], int m, int n, int i, double v[]);
void updateRow(double (*A)[], int m, int n, int j, double v[]);
void updateEntry(double x, double (*A)[], int m, int n, int i, int j);
double printMax(double v[], int n);
void qrGS1(double A[], int m, int n, double (*Q)[], double (*R)[]);
void qrGS2(double A[], int m, int n, double (*Q)[], double (*R)[]);
void cholFac(double A[], int m, double (*R)[]);
void backsub(double A[], int n, double b[], double (*x)[]);
void forwardsub(double A[], int n, double b[], double (*x)[]);
void cholSolve(double A[], int n, double (*x)[], double b[]);
void qrSolve(double A[], int n, double (*x)[], double b[]);
void gsSolve(double A[], int n, double (*x)[], double b[]);
void hhSolve(double A[], int n, double (*x)[], double b[]);
void backsub(double A[], int n, double b[], double (*x)[]);
double testSolve(double A[], int n, double x[], double b[]);
void submat(double A[], int m, int n, int i1, int i2, 
            int j1, int j2,double (*As)[]);
int sign(double x);
void insertSubmat(double (*A)[], int m, int n, int i1, int i2, 
                  int j1, int j2, double As[]);
void qrHH(double A[], int m, int n, double (*Q)[], double (*R)[]);
