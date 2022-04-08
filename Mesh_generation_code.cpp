
/*
CC-297 - Projeto 2 - Geração de Malha Computacional
Programa Desenvolvido por : Alisson Vinicius Brito Lopes
*/

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <string>
#include <sys/time.h>

# define M_PI 3.14159265358979323846  /* pi */

using namespace std;


class WriteFile{ // Classe Malreve Arquivos

public:

static void writeMatrix(double **M, int qtdlinhas, int qtdColunas, char nomeArquivo[]){

     ofstream write;
     write.open(nomeArquivo);

     for (int i = 0; i < qtdlinhas; i++ ){

        write << "\n";

        for (int j = 0; j < qtdColunas; j++){

             write.width(10);
             write << M[i][j] << "        ";

        }
    }
             write.close();
    }


static void writeVetAgrup(double *V1, double *V2, int qtdlinhas, char nomeArquivo[]){

    double MatAgrVet[qtdlinhas][2];

    for(int i=0; i<qtdlinhas;i++){
        MatAgrVet[i][0]=V1[i];
    }

    for(int i=0; i<qtdlinhas;i++){
        MatAgrVet[i][1]=V2[i];
    }

     ofstream write;
     write.open(nomeArquivo);
       //variables="x" "y"
  //zone I=          93  ,J=          15  ,K=           1  F=point

  write << " variables = ''x'' ''y'' " << "\n";
  write << "zone I=          93  ,J=          15  ,K=           1  F=point"<<endl;

     for (int i = 0; i < qtdlinhas; i++ ){

        write << "\n";

        for (int j = 0; j < 2; j++){

             write.width(11);
             write << MatAgrVet[i][j] << " ";

        }
    }
             write.close();

}


static void writeVetor(double *V, int qtdColunas, char nomeArquivo[]){

     ofstream write;
     write.open(nomeArquivo);

     for (int i = 0; i < qtdColunas; i++ ){

            write << V[i] <<"\n";

        }
             write.close();
}

};

class Vetor{ // Classe gera Vetores

    public:

    static void imprimir(double *V, int qtdColunas){

    for (int i = 0; i < qtdColunas; i++ ){

        cout << V[i] << "\t";
    }
        cout << endl;
    }


};

class Matrix{ // Classe gera e imprime

public:

   static void  imprimir(double **M,int qtdLinhas, int qtdColunas){ // Nome de objeto minusculo

        for ( int i = 0; i < qtdLinhas; i++){

                Vetor::imprimir(M[i],qtdColunas);
        }
    }

    static void imprimirTecPlot(double **M,int qtdLinhas, int qtdColunas){

        for (int j = 0; j < qtdColunas; j++ ){
            for (int i = 0; i < qtdLinhas; i++){
                    cout << M[i][j] << "\n";
            }
         }
    }

};

class Malha {

    public:

    int IMAX, JMAX, ILE, ITE;
    double XSF, YSF;

// Criando ponteiro duplos e simples para Alocação dinâmica de memoria e "zeramento" de vetores e matrizes na memoria do computador

    double **Lx,**Ly, **Cc, **Phi, **f, **g,**DELTAx, **DELTAy;

    double *RHS, *AT, *BT, *CT, *DT, *ATX, *BTX, *CTX, *VTX,*ATY, *BTY, *CTY, *VTY;

    double **A, **B, **C, **D, **P, **Q, **AMX, **AMY;

    double **x,**y;

    double **xint, **yint, **xo,**yo;

    double *eref,*s;

    double *Rx, *R, *Ry, *Seta, *Dels;



void inicializaMatriz_Dels( int JMAX){
   int j;

    Dels = new double [JMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < JMAX; j++){

           Dels[j] = 0.0;
    }
}

void inicializaMatriz_Seta( int JMAX){
   int j;

    Seta = new double [JMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < JMAX; j++){

           Seta[j] = 0.0;
    }
}


void inicializaMatriz_Rx( int JMAX){
   int j;

    Rx = new double [JMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < JMAX; j++){

           Rx[j] = 0.0;
    }
}

void inicializaMatriz_Ry( int JMAX){
   int j;

    Ry = new double [JMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < JMAX; j++){

           Ry[j] = 0.0;
    }
}

void inicializaMatriz_R( int JMAX){
   int j;

    R = new double [JMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < JMAX; j++){

           R[j] = 0.0;
    }
}


void inicializaMatriz_s( int JMAX){
   int j;

    s = new double [JMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < JMAX; j++){

           s[j] = 0.0;
    }
}


void inicializaMatriz_eref( int JMAX){
   int j;

    eref = new double [JMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < JMAX; j++){

           eref[j] = 0.0;
    }
}


void inicializaMatriz_xint(int IMAX, int JMAX){ // MATRIX DE RESIDUOS

  int i, j;

  xint = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        xint[i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 xint[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_yint(int IMAX, int JMAX){ // MATRIX DE RESIDUOS

  int i, j;

  yint = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        yint[i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 yint[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_xo(int IMAX, int JMAX){ // MATRIX DE RESIDUOS

  int i, j;

  xo = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        xo[i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 xo[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_yo(int IMAX, int JMAX){ // MATRIX DE RESIDUOS

  int i, j;

  yo = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        yo[i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 yo[i][j] = 0.0;
         }
    }
}


void inicializaMatriz_x(int IMAX, int JMAX){ // MATRIX DE RESIDUOS

  int i, j;

  x = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        x[i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 x[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_y(int IMAX, int JMAX){ // MATRIX DE RESIDUOS

  int i, j;

  y = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        y[i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 y[i][j] = 0.0;
         }
    }
}


void inicializaMatriz_Lx(int IMAX, int JMAX){ // MATRIX DE RESIDUOS

  int i, j;

  Lx = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        Lx[i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 Lx[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_Ly(int IMAX, int JMAX){ // MATRIX DE RESIDUOS

  int i, j;

  Ly = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        Ly[i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 Ly[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_Phi(int IMAX, int JMAX){ // POTENCIAL Zerando e ja assumindo chute inicial

  int i, j;

  Phi = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        Phi [i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 Phi [i][j] = 0.0;
         }
    }
}

void inicializaMatriz_Cc(int IMAX, int JMAX){ // UPDATE Zerando e ja assumindo chute inicial

  int i, j;

  Cc = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        Cc [i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                Cc [i][j] = 0.0;
             }
    }
}

// Inicializando e alocando espaço para as variaveis calculadas no algoritmo de Thomas
void inicializaMatriz_AT( int IMAX){
   int j;

    AT = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < IMAX; j++){

           AT[j] = 0.0;
    }
}

void inicializaMatriz_BT( int IMAX){
   int j;

    BT = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < IMAX; j++){

           BT[j] = 0.0;
    }
}

void inicializaMatriz_CT( int IMAX){
   int j;

    CT = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < IMAX; j++){

           CT[j] = 0.0;
    }
}

void inicializaMatriz_DT( int IMAX){
   int j;

    DT = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j <IMAX; j++){

           DT[j] = 0.0;
    }
}


void inicializaMatriz_ATX( int IMAX){
   int j;

    ATX = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < IMAX; j++){

           ATX[j] = 0.0;
    }
}

void inicializaMatriz_BTX( int IMAX){
   int j;

    BTX = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < IMAX; j++){

           BTX[j] = 0.0;
    }
}

void inicializaMatriz_CTX( int IMAX){
   int j;

    CTX = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < IMAX; j++){

           CTX[j] = 0.0;
    }
}

void inicializaMatriz_VTX( int IMAX){
   int j;

    VTX = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j <IMAX; j++){

           VTX[j] = 0.0;
    }
}


void inicializaMatriz_ATY( int IMAX){
   int j;

    ATY = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < IMAX; j++){

           ATY[j] = 0.0;
    }
}

void inicializaMatriz_BTY( int IMAX){
   int j;

    BTY = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < IMAX; j++){

           BTY[j] = 0.0;
    }
}

void inicializaMatriz_CTY( int IMAX){
   int j;

    CTY = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < IMAX; j++){

           CTY[j] = 0.0;
    }
}

void inicializaMatriz_VTY( int IMAX){
   int j;

    VTY = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j <IMAX; j++){

           VTY[j] = 0.0;
    }
}


// Matriz dos Coeficientes da Equação 1

void inicializaMatriz_A(int IMAX, int JMAX){

  int i, j;

  A = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        A[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 A[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_B(int IMAX, int JMAX){

  int i, j;

  B = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        B[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 B[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_C(int IMAX, int JMAX){

  int i, j;

  C = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        C[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 C[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_D(int IMAX, int JMAX){

  int i, j;

  D = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        D[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 D[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_P(int IMAX, int JMAX){

  int i, j;

  P = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        P[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 P[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_Q(int IMAX, int JMAX){

  int i, j;

  Q = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        Q[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 Q[i][j] = 0.0;
         }
    }
}


void inicializaMatriz_AMX(int IMAX, int JMAX){ // UPDATE Zerando e ja assumindo chute inicial

  int i, j;

  AMX = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        AMX[i] = new double[IMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < IMAX; j++){
                AMX[i][j] = 0.0;
             }
    }
}


void inicializaMatriz_AMY(int IMAX, int JMAX){ // UPDATE Zerando e ja assumindo chute inicial

  int i, j;

  AMY = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        AMY[i] = new double[IMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < IMAX; j++){
                AMY[i][j] = 0.0;
             }
    }
}

void inicializaMatriz_f(int IMAX, int JMAX){ // MATRIX DE RESIDUOS

  int i, j;

  f = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        f[i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 f[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_g(int IMAX, int JMAX){ // MATRIX DE RESIDUOS

  int i, j;

  g = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        g[i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 g[i][j] = 0.0;
         }
    }
}



void inicializaMatriz_DELTAx(int IMAX, int JMAX){ // MATRIX DE RESIDUOS

  int i, j;

  DELTAx = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        DELTAx[i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 DELTAx[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_DELTAy(int IMAX, int JMAX){ // MATRIX DE RESIDUOS

  int i, j;

  DELTAy = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        DELTAy[i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 DELTAy[i][j] = 0.0;
         }
    }
}


void inicializarMatrizes(int IMAX, int JMAX,double th,double alfa, double omega, double eps, double deltaX){


    // Variaveis A, B, C, D, P, Q da Equação 1 Escritas no espaço Computacional
    inicializaMatriz_A(IMAX,JMAX);
    inicializaMatriz_B(IMAX,JMAX);
    inicializaMatriz_C(IMAX,JMAX);
    inicializaMatriz_D(IMAX,JMAX);
    inicializaMatriz_P(IMAX,JMAX);
    inicializaMatriz_Q(IMAX,JMAX);

    // Inicialização de matrizes para geração da Mal computacional
    inicializaMatriz_Phi(IMAX, JMAX);
    inicializaMatriz_Cc(IMAX, JMAX);

    // Inicialização para o Algoritmo de THOMAS
    inicializaMatriz_AT(IMAX);
    inicializaMatriz_BT(IMAX);
    inicializaMatriz_CT(IMAX);
    inicializaMatriz_DT(IMAX);


    // Inicialização para periodica em X
    inicializaMatriz_ATX(IMAX);
    inicializaMatriz_BTX(IMAX);
    inicializaMatriz_CTX(IMAX);
    inicializaMatriz_VTX(IMAX);

    // Inicialização para periodica em Y
    inicializaMatriz_ATY(IMAX);
    inicializaMatriz_BTY(IMAX);
    inicializaMatriz_CTY(IMAX);
    inicializaMatriz_VTY(IMAX);

    // Operador de Resuldo direção xi,j e yi,j
    inicializaMatriz_Lx(IMAX, JMAX);
    inicializaMatriz_Ly(IMAX, JMAX);

    inicializaMatriz_x(IMAX,JMAX);
    inicializaMatriz_y(IMAX,JMAX);

    inicializaMatriz_xint(IMAX,JMAX);
    inicializaMatriz_yint(IMAX,JMAX);

    inicializaMatriz_xo(IMAX,JMAX);
    inicializaMatriz_yo(IMAX,JMAX);

    inicializaMatriz_eref(JMAX);

    // SOMa da PG
    inicializaMatriz_s(JMAX);

    // Inicializando variaveis para Geração da Malha Local de referencia
    inicializaMatriz_Rx(JMAX);
    inicializaMatriz_Ry(JMAX);
    inicializaMatriz_R(JMAX);

    inicializaMatriz_Seta(JMAX);
    inicializaMatriz_Dels(JMAX);
    inicializaMatriz_eref(JMAX);

    inicializaMatriz_AMX(IMAX,JMAX);
    inicializaMatriz_AMY(IMAX,JMAX);

    inicializaMatriz_f(IMAX,JMAX);
    inicializaMatriz_g(IMAX,JMAX);

    inicializaMatriz_DELTAx(IMAX,JMAX);
    inicializaMatriz_DELTAy(IMAX,JMAX);

 }

};

//double MaxResiduo (Malha Mal, int i, int j,int IMAX,int JMAX,double deltaX,int kiter, double resmax);
Malha aerofolioPontosBiconvexo (Malha Mal, int i, int j, int IMAX, int JMAX, double R,double th, double XSF);
Malha malhaInicialAlg1 (Malha Mal,int IMAX, int JMAX, double R,double XSF, double YSF, int refj, int refjaux);
Malha malhaInicialParabolica (Malha Mal, int i, int j, int IMAX, int JMAX, double R,double XSF, double YSF, int refj);
void TridiagPX(double* a, double* b, double* c, double* dv, double *x, int n);
void TridiagPY(double* a, double* b, double* c, double* dv, double *y, int n);
Malha SLOR (Malha Mal,int i, int j, int IMAX, int JMAX,double deltaX, double r,double omega);
Malha AF1 (Malha Mal,int i, int j, int IMAX, int JMAX, double omega, double alfa);
Malha AF2 (Malha Mal,int i, int j, int IMAX, int JMAX, double omega, double alfa);
void TDMAX(double* a, double* b, double* c, double* d, double *x, int n);
void TDMAY(double* a, double* b, double* c, double* d, double *y, int n);
double MaxResiduoX (Malha Mal,int IMAX,int JMAX,int kiter, double resmaxX);
double MaxResiduoY (Malha Mal,int IMAX,int JMAX,int kiter, double resmaxY);
float FloatMachineEps();
double DoubleMachineEps();
long double LongDoubleMachineEps();
void AlfaSequence (int Melements, double alfaLow, double alfaHigh);
Malha aerofolioPontosNACA (Malha Mal, int i, int j, int IMAX, int JMAX,double R,double th, double XSF);
Malha AF2Dirichlet (Malha Mal,int i, int j, int IMAX, int JMAX, double omega, double alfa);
Malha AF2Atraction (Malha Mal,int i, int j, int IMAX, int JMAX, double omega, double alfa);
void TransMatrixVetor(double** Mx,double** My, int IMAX, int JMAX, double* xVet, double* yVet);
void PQfunction(Malha Mal,double* a, double* b, double* c, double* d, int i, int j, int L, int MM);
//double sgn(Malha Mal,double **M);
double valorAbsoluto(double x);



int main(int argc, char*argv[]){

// We can determine the value of machine epsilon by finding the value of 1.0 +1.0/2^p = 1.0

//float fep;
//double dep;
//long double ldep;
//
//fep = FloatMachineEps();
//dep = DoubleMachineEps();
//ldep = LongDoubleMachineEps();

//cout << "Machine epsilon for Single precision is:" << fep << endl;
//cout << "Machine epsilon for Double precision is:" << dep << endl;
//cout << "Machine epsilon for Long Double precision is:" << ldep << endl;

// Especificações Técnicas do computador utilizado:
// Intel(R) Core(TM) i5-3210M CPU @ 2.50 GHz, sistema operacional de 64 Bits

// "Machine epsilon for Single precision is = 5.96046e-008 ;
// "Machine epsilon for Double precision is = 11.1022e-015 ;
// "Machine epsilon for Long double precision is = 5.42101e-020 ;

// Contador de Tempo computacional para avaliação dos metodos

clock_t TempoInicial;
TempoInicial = clock();

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         GRID PARAMETERS                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

//double res;
double  resmaxX = 0.0, resmaxY = 0.0;
int maxit = 10000;

int IMAX = 93; // Número maximo de pontos na direção Csi >> i
int JMAX = 15; // Número maximo de pontos na direção Eta >> j

// Vetores para Transformação da Matrix para imprimir no formato Fortran (permite impressão tecplot)
double* xVet = new double[IMAX*JMAX];
double* yVet = new double[IMAX*JMAX];

double *resmaxiterx, *resmaxitery, RMAXX = 1.0,RMAXY = 1.0;
resmaxiterx = new double [maxit];
resmaxitery = new double [maxit];

int i,j,kiter; // i direção csi e j direção eta

double XSF = 1.25; // fator de estiramento ("Stretching") da Mal para a direção x
double YSF = 1.25; // fator de estiramento ("Stretching") da Mal para a direção y %%%%%%% USAR 1.19 para todas as simulações

 /*%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    DADOS DE ENTRADA    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

// Iterations parameters for Numerical Methods and Iteration Parameters for Eingevalues anihhilation
double Melements = 6;
double alfaLow = 1.0;
double alfaHigh = 154.0;


double r = 1.9; // Fator de Relaxação do SLOR

// Raio do Circulo externo
double R = 6.5;

// Airfoil thickness ratio
double th = 0.10;

// Iterations parameters for Numerical Methods AF1
//double alfa = 0.1299; // Combinação Boa para meu problema 0.1299 e omega = 1,90 discretização sala
//double omega = 1.9;

// Iterations parameters AF2
double alfa = 0.05; // Combinação Boa para meu problema 0.1299 e omega = 1,90 discretização sala
double omega = 1.9;

//Iterations Parameters for Eingenvalus anihhilation

// Convergence Criterion
double eps = 1.11022*pow(10,-12); //

// Calculo do deltaX na direção X apenas no intervalo em cima do perfil biconvexo
double deltaX;

Malha Mal;
Mal.inicializarMatrizes(IMAX,JMAX,th,alfa,omega,eps,deltaX);

Mal = aerofolioPontosBiconvexo (Mal,i,j,IMAX,JMAX,R,th,XSF);
//Mal = aerofolioPontosNACA (Mal,i,j,IMAX,JMAX,R,th,XSF);

//AlfaSequence (Melements,alfaLow,alfaHigh);


// LOOP GERAR TODAS AS LINHAS DA MALHA PARABOLICA
for (int refj = 1; refj<=JMAX-2;refj++){

if (refj >= JMAX-2){

    int refjaux = refj-1; // refjaux existe apenas para permitir a geração da penultima linha algebrica/parabolica

    Mal = malhaInicialAlg1 (Mal,IMAX,JMAX,R,XSF,YSF,refj,refjaux);
    Mal = malhaInicialParabolica (Mal,i,j,IMAX,JMAX,R,XSF,YSF,refj);

    }else{

    int refjaux = refj;

    Mal = malhaInicialAlg1 (Mal,IMAX,JMAX,R,XSF,YSF,refj,refjaux);
    Mal = malhaInicialParabolica (Mal,i,j,IMAX,JMAX,R,XSF,YSF,refj);

    }
}

//char nomeArquivo3[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto2/YintNEWNACA.txt";
//WriteFile::writeMatrix(Mal.y,IMAX,JMAX,nomeArquivo3);
//
//char nomeArquivo4[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto2/XintNEWNACA.txt";
//WriteFile::writeMatrix(Mal.x,IMAX,JMAX,nomeArquivo4);

// Imprimi Malha Parabolica
TransMatrixVetor(Mal.x,Mal.y,IMAX,JMAX,xVet,yVet);
char nomeArquivoParabolica[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto2/Graficos/MalhaParabolica.plt";
WriteFile::writeVetAgrup(xVet,yVet,(IMAX*JMAX),nomeArquivoParabolica);


// LOOP PARA CALCULO ITERATIVO DA MALHA ELIPTICA
for (kiter = 0; kiter < maxit; kiter ++){

    //Mal = AF1 (Mal,i,j,IMAX,JMAX,omega,alfa);
    //Mal = AF2 (Mal,i,j,IMAX,JMAX,omega,alfa);
    // Mal = AF2Dirichlet (Mal,i,j,IMAX,JMAX,omega,alfa);
     Mal = AF2Atraction (Mal,i,j,IMAX,JMAX,omega,alfa);
    //Mal = SLOR (Mal,i, j,IMAX,JMAX,deltaX,r,omega);


    RMAXX = MaxResiduoX (Mal,IMAX,JMAX,kiter,resmaxX); //MaxResiduo(esc, i, j, ILE, ITE, IMAX, JMAX, XSF, YSF, deltaX,kiter,resmax);
    RMAXY = MaxResiduoY (Mal,IMAX,JMAX,kiter,resmaxY); //MaxResiduo(esc, i, j, ILE, ITE, IMAX, JMAX, XSF, YSF, deltaX,kiter,resmax);

    resmaxiterx[kiter] =log10(RMAXX);
    resmaxitery[kiter] =log10(RMAXY);

    //kiter = kiter +1;

    //cout << RMAXX << endl;
    //cout << kiter << endl;
    cout << "ResX = " << resmaxiterx[kiter] << endl;
    //cout << "ResY = " << resmaxitery[kiter] << endl;


}



//TransMatrixVetor(Mal.x,Mal.y, IMAX, JMAX, xVet, yVet);
TransMatrixVetor(Mal.Lx,Mal.Ly,IMAX,JMAX,xVet,yVet);

char nomeArquivo3[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto2/YintesteNACA.plt";
WriteFile::writeVetAgrup(xVet,yVet,(IMAX*JMAX),nomeArquivo3);

//(yVet,IMAX,JMAX,nomeArquivo3);

//char nomeArquivo4[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto2/Xinteste.txt";
//WriteFile::writeMatrix(xVet,IMAX,JMAX,nomeArquivo4);

//char nomeArquivo3[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto2/Yint1iterBiconDirichlet.txt";
//WriteFile::writeMatrix(Mal.y,IMAX,JMAX,nomeArquivo3);
//
//char nomeArquivo4[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto2/Xint1iterBiconDirichlet.txt";
//WriteFile::writeMatrix(Mal.x,IMAX,JMAX,nomeArquivo4);

//char nomeArquivo3[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto2/Yint1iter.txt";
//WriteFile::writeMatrix(Mal.y,IMAX,JMAX,nomeArquivo3);

//char nomeArquivo4[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto2/Xint1iter.txt";
//WriteFile::writeMatrix(Mal.x,IMAX,JMAX,nomeArquivo4);

//char nomeArquivo5[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto2/ResiduoXa05w19.txt";
//WriteFile::writeVetor(resmaxiterx,maxit,nomeArquivo5);
//
//char nomeArquivo6[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto2/ResiduoYa05w19.txt";
//WriteFile::writeVetor(resmaxitery,maxit,nomeArquivo6);


double tempo_total_s = (clock()- TempoInicial)/(double)CLOCKS_PER_SEC;
cout << "Tempo Total de Execucao = " << tempo_total_s << endl;

   // Mal = SLOR(Mal,i,j,IMAX,JMAX,deltaX,r);

return 0;

};


Malha aerofolioPontosBiconvexo (Malha Mal, int i, int j, int IMAX, int JMAX,double R,double th, double XSF){


    // Gerando o circulo exterior
    for ( i = 0; i <(IMAX-1); i++){ //  IMAX-1 = 92  >>> vai de 0 - 91

    Mal.x[i][JMAX-1] = R*cos((2*M_PI*(double(i)/(IMAX-1)))); //(IMAX-2/IMAX-1)
    Mal.y[i][JMAX-1] = -R*sin((2*M_PI*(double(i)/(IMAX-1)))); // coloquei menos para ser sentido horario

    }

    Mal.x[IMAX-1][JMAX-1] = Mal.x[0][JMAX-1];
    Mal.y[IMAX-1][JMAX-1] = Mal.y[0][JMAX-1];



    // Fronteira Interna com função Cosseno
//    for ( i= 0; i<=((IMAX-1)/2);i++){
//
//        Mal.x[i][0] = 0.5*cos(2*M_PI*(double(i))/(IMAX-1))+ 0.5;
//        Mal.y[i][0] = -2.0*th*Mal.x[i][0]*(1.0 - Mal.x[i][0]);
//        Mal.x[i][0] = Mal.x[i][0]- 0.5;
//
//    }
//
//   int k = ((IMAX-1)/2);
//    // Gerando os valores de y para o terceiro e quarto quadrantes
//    for (i = ((IMAX-1)/2); i <IMAX-1; i++){
//
//        Mal.y[i][0] = -Mal.y[k][0];//2.0*th*Mal.x[i][0]*(1.0 - Mal.x[i][0]);//- Mal.y[k][0];
//        Mal.x[i][0] =  Mal.x[k][0];
//        k = k-1;
//    }
//
//     Mal.x[IMAX-1][0] = Mal.x[0][0];
//     Mal.y[IMAX-1][0] = Mal.y[0][0];


 //Implementação com Estiramento usando Progressão Geometrica (esta sera apresentada no relatorio)

        int numPointsPG = ((IMAX-1)/4);
        double somaPG = 0.5;
        double dx = somaPG*(XSF-1)/(pow(XSF,numPointsPG)-1.0);

        Mal.x[0][0] = 0.0;
        Mal.x[1][0] = Mal.x[0][0] + dx;

        Mal.y[0][0] = 2*th*Mal.x[0][0]*(1.0 -Mal.x[0][0]);
        Mal.y[1][0] = -2*th*Mal.x[1][0]*(1.0 -Mal.x[1][0]);

        Mal.x[(IMAX-1)/2][0] = 1.0;
        Mal.x[((IMAX-1)/2) -1][0] = Mal.x[(IMAX-1)/2][0] - dx;

        Mal.y[(IMAX-1)/2][0] = -2*th*Mal.x[(IMAX-1)/2][0]*(1.0 - Mal.x[(IMAX-1)/2][0]);
        Mal.y[(IMAX-1)/2 -1][0] = -2*th*Mal.x[((IMAX-1)/2) -1][0]*(1.0 - Mal.x[((IMAX-1)/2) -1][0]);


        //cout << Mal.x[0][0]<< endl;
        //cout << Mal.x[1][0]<< endl;

    for ( i = 2; i <=((IMAX-1)/4); i++){ // OK
    // Gerando a PG para estiramento x aerofolio

        Mal.x[i][0] = (Mal.x[i-1][0] + XSF*(Mal.x[i-1][0] - Mal.x[i-2][0]));
        Mal.y[i][0] = -2*th*Mal.x[i][0]*(1.0 - Mal.x[i][0]);

       // cout << Mal.x[i][0] <<endl;

    }

    // Gerando o segundo quadrante
    for ( i =((IMAX-1)/2 -2); i >=((IMAX-1)/4 +1) ; i-- ){
        Mal.x[i][0] = Mal.x[i+1][0] - (Mal.x[i+2][0] - Mal.x[i+1][0])*XSF;
        Mal.y[i][0] = -2*th*Mal.x[i][0]*(1.0 - Mal.x[i][0]);
        //cout << Mal.x[i][0] <<endl;
      }


    for (i = 0; i<=((IMAX-1)/2);i++){ // ESSE LAÇO GERA O X estirado certo de 0 até (imax-1)/2

        Mal.x[i][0] = Mal.x[i][0] - 0.5;
        Mal.x[i][0] = - Mal.x[i][0];
    }


    int k = ((IMAX-1)/2);
    // Gerando os valores de y para o terceiro e quarto quadrantes
    for (i = (((IMAX-1)/2)); i <=IMAX-2; i++){ // era i <=IMAX-1 e agora é : i <=IMAX-2

        Mal.y[i][0] = - Mal.y[k][0];
        Mal.x[i][0] =  Mal.x[k][0];

        //cout <<  Mal.y[i][0] << endl;
        //cout <<  Mal.x[i][0] << endl;

        k = k-1;
    }

    Mal.x[IMAX-1][0] = Mal.x[0][0];
    Mal.y[IMAX-1][0] = Mal.y[0][0];


//char nomeArquivo3[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto2/YintNEW.txt";
//WriteFile::writeMatrix(Mal.y,IMAX,JMAX,nomeArquivo3);
//
//char nomeArquivo4[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto2/XintNEW.txt";
//WriteFile::writeMatrix(Mal.x,IMAX,JMAX,nomeArquivo4);


return Mal;

}

Malha malhaInicialParabolica (Malha Mal, int i, int j, int IMAX, int JMAX, double R,double XSF, double YSF, int refj){

        double xqsi, yqsi, xeta, yeta;
        double *x; // Armazenando espaço e zerando o vetor solução da periodica
        x = new double [IMAX];

        double *y; // Armazenando espaço e zerando o vetor solução da periodica
        y = new double [IMAX];
        int M = IMAX-2; // Penultimo ponto da periodica

    int jj=refj;

    for (i = 0; i <=M; i++){

             if (i == 0){

        xqsi = 0.5*(Mal.x[i+1][jj+1] - Mal.x[IMAX-2][jj+1]);
        yqsi = 0.5*(Mal.y[i+1][jj+1] - Mal.y[IMAX-2][jj+1]);
        xeta = Mal.x[i][jj+1] - Mal.x[i][jj];
        yeta = Mal.y[i][jj+1] - Mal.y[i][jj];

        Mal.A[i][jj] = xeta*xeta + yeta*yeta;
        Mal.B[i][jj] = xqsi*xeta + yqsi*yeta;
        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;

        Mal.AT [M] = 2.0*Mal.A[i][jj];
        Mal.BT [i] = -4.0*(Mal.A[i][jj] + Mal.C[i][jj]);
        Mal.CT [i+1] = 2.0*Mal.A[i][jj];
        Mal.VTX [i] = Mal.B[i][jj]*(Mal.x[i+1][jj+1] - Mal.x[IMAX-2][jj+1] - Mal.x[i+1][jj-1] + Mal.x[IMAX-2][jj-1]) - 2.0*Mal.C[i][jj]*(Mal.x[i][jj+1] + Mal.x[i][jj-1]);
        Mal.VTY [i] = Mal.B[i][jj]*(Mal.y[i+1][jj+1] - Mal.y[IMAX-2][jj+1] - Mal.y[i+1][jj-1] + Mal.y[IMAX-2][jj-1]) - 2.0*Mal.C[i][jj]*(Mal.y[i][jj+1] + Mal.y[i][jj-1]);


            }else{

        xqsi = 0.5*(Mal.x[i+1][jj+1] - Mal.x[i-1][jj+1]);
        yqsi = 0.5*(Mal.y[i+1][jj+1] - Mal.y[i-1][jj+1]);
        xeta = Mal.x[i][jj+1] - Mal.x[i][jj];
        yeta = Mal.y[i][jj+1] - Mal.y[i][jj];

        Mal.A[i][jj] = xeta*xeta + yeta*yeta;
        Mal.B[i][jj] = xqsi*xeta + yqsi*yeta;
        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;

        Mal.AT [i-1] = 2.0*Mal.A[i][jj];
        Mal.BT [i] = -4.0*(Mal.A[i][jj] + Mal.C[i][jj]);
        Mal.CT [i+1] = 2.0*Mal.A[i][jj];

        Mal.VTX [i] = Mal.B[i][jj]*(Mal.x[i+1][jj+1] - Mal.x[i-1][jj+1] - Mal.x[i+1][jj-1] + Mal.x[i-1][jj-1]) - 2.0*Mal.C[i][jj]*(Mal.x[i][jj+1] + Mal.x[i][jj-1]); // LADO DIREITO RHS
        Mal.VTY [i] = Mal.B[i][jj]*(Mal.y[i+1][jj+1] - Mal.y[i-1][jj+1] - Mal.y[i+1][jj-1] + Mal.y[i-1][jj-1]) - 2.0*Mal.C[i][jj]*(Mal.y[i][jj+1] + Mal.y[i][jj-1]); // LADO DIREITO RHS

        }
    }

        Mal.CT [0] = 2.0*Mal.A[M][jj];

        TridiagPX(Mal.AT, Mal.BT, Mal.CT, Mal.VTX, x, M);
        TridiagPY(Mal.AT, Mal.BT, Mal.CT, Mal.VTY, y, M);

        for (i =0; i<=M; i++){

            Mal.x[i][jj] = x[i];
            Mal.y[i][jj] = y[i];
        }

        Mal.x[IMAX-1][jj] = Mal.x[0][jj];
        Mal.y[IMAX-1][jj] = Mal.y[0][jj];

        //refj = refj + 1;

return Mal;


}

Malha malhaInicialAlg1 (Malha Mal,int IMAX, int JMAX, double R,double XSF, double YSF, int refj, int refjaux){

    double eps,dy, deltaS, RR, Rx, Ry, seta,xqsi, yqsi,xeta,yeta; // Where R is the distance between the point i,j-1 and the outer boundary i, JMAX
    double xort, yort, xint, yint,ds;
     //Gerando a PG para o estiramento na direção eta

    int numPointsPGeta = JMAX;
    double somaPGeta = 1.0;
    dy = somaPGeta*(YSF-1)/(pow(YSF,numPointsPGeta-1)-1.0);

    Mal.s[0] = 0.0;
    Mal.s[1] = Mal.s[0] + dy;

     for (int j = 1; j < JMAX-1; j++ ){

        Mal.s[j+1] = Mal.s[j] + YSF*(Mal.s[j] - Mal.s[j-1]);

     }

    for (int j = 1; j < JMAX-1; j++ ){

         Mal.Dels[j] = (Mal.s[j] - Mal.s[j-1])/(Mal.s[JMAX-1] -Mal.s[j-1]);
     }


    for (int j = refj; j<=(refjaux+1);j++){ // j<=(refjaux+1)

        for (int i = 0; i <IMAX-1; i++){

        Rx = pow( (Mal.x[i][JMAX-1] - Mal.x[i][j-1]),2);
        Ry = pow( (Mal.y[i][JMAX-1] - Mal.y[i][j-1]),2);
        R  = pow((Rx+Ry),0.5);

        seta = Mal.Dels[j]*R;

        //cout << setprecision(10) << setiosflags(ios::scientific)<< endl;
        //cout << seta << endl;

        //cout << Mal.Dels[j] << endl;

        if (i==0){

        xqsi= ( Mal.x[i+1][j-1]- Mal.x[IMAX-2][j-1] )/2.0; //Mal.x[0][0] = Mal.x[IMAX-2][j]
        yqsi= ( Mal.y[i+1][j-1]- Mal.y[IMAX-2][j-1] )/2.0; //Mal.y[0][0] = Mal.y[IMAX-2][j]

        }else{

        xqsi= ( Mal.x[i+1][j-1]- Mal.x[i-1][j-1] )/2.0; //Mal.x[0][0] = Mal.x[IMAX-2][j] Ora pois : é periodico !!!
        yqsi= ( Mal.y[i+1][j-1]- Mal.y[i-1][j-1] )/2.0; //Mal.y[0][0] = Mal.y[IMAX-2][j]

        }

        xeta = -yqsi*seta/(pow((xqsi*xqsi + yqsi*yqsi),0.5));
        yeta =  xqsi*seta/(pow((xqsi*xqsi + yqsi*yqsi),0.5));

        xort = Mal.x[i][j-1] + xeta;
        yort = Mal.y[i][j-1] + yeta;

        xint = Mal.x[i][j-1] + Mal.Dels[j]*(Mal.x[i][JMAX-1]-Mal.x[i][j-1]);
        yint = Mal.y[i][j-1] + Mal.Dels[j]*(Mal.y[i][JMAX-1]-Mal.y[i][j-1]);

        eps= double(j)/(JMAX - 1);

        Mal.x[i][j] = eps*xint + (1.0 - eps)*xort;
        Mal.y[i][j] = eps*yint + (1.0 - eps)*yort;


    }

        Mal.x[IMAX-1][j] = Mal.x[0][j];
        Mal.y[IMAX-1][j] = Mal.y[0][j];

}


//Vetor::imprimir(Mal.Dels,JMAX);

//char nomeArquivo3[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto2/Yint.txt";
//WriteFile::writeMatrix(Mal.y,IMAX,JMAX,nomeArquivo3);
//
//char nomeArquivo4[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto2/Xint.txt";
//WriteFile::writeMatrix(Mal.x,IMAX,JMAX,nomeArquivo4);

//char nomeArquivo1[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto2/Yint.txt";
//WriteFile::writeMatrix(Mal.y,IMAX,JMAX,nomeArquivo1);
//
//char nomeArquivo2[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto2/Xint.txt";
//WriteFile::writeMatrix(Mal.x,IMAX,JMAX,nomeArquivo2);


return Mal;

}

Malha SLOR (Malha Mal, int i, int j,int IMAX,int JMAX,double deltaX, double r,double omega){

    double xqsi, yqsi, xeta, yeta, delqsiqsiX, delqsietaX, deletaetaX, delqsiqsiY, delqsietaY, deletaetaY;

    double *xx; // Armazenando espaço e zerando o vetor solução da periodica
    xx = new double [IMAX];

    double *yy; // Armazenando espaço e zerando o vetor solução da periodica
    yy = new double [IMAX];

    int M = IMAX-2; // Penultimo ponto da periodica


for (int jj = 1; jj<=JMAX-2; jj++){

    for (i = 0; i <=M; i++){

             if (i == 0){

        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[IMAX-2][jj]);
        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[IMAX-2][jj]);
        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);

        Mal.A[i][jj] = xeta*xeta + yeta*yeta;
        Mal.B[i][jj] = xqsi*xeta + yqsi*yeta;
        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;

        Mal.AT [M] = Mal.A[i][jj];
        Mal.BT [i] = -2.0*(1+ Mal.A[i][jj]);
        Mal.CT [i+1] = Mal.A[i][jj];

        delqsiqsiX= Mal.x[i+1][jj] -2.0*Mal.x[i][jj] + Mal.x[IMAX-2][jj];
        delqsietaX = 0.25*(Mal.x[i+1][jj+1] - Mal.x[i+1][jj-1] - Mal.x[IMAX-2][jj+1] + Mal.x[IMAX-2][jj-1]);
        deletaetaX = Mal.x[i][jj+1] -2.0*Mal.x[i][jj] + Mal.x[i][jj-1];

        delqsiqsiY= Mal.y[i+1][jj] -2.0*Mal.y[i][jj] + Mal.y[IMAX-2][jj];
        delqsietaY = 0.25*(Mal.y[i+1][jj+1] - Mal.y[i+1][jj-1] - Mal.y[IMAX-2][jj+1] + Mal.y[IMAX-2][jj-1]);
        deletaetaY = Mal.y[i][jj+1] -2.0*Mal.y[i][jj] + Mal.y[i][jj-1];

        Mal.Lx[i][jj] = (Mal.A[i][jj]*delqsiqsiX -2.0*Mal.B[i][jj]*delqsietaX + Mal.C[i][jj]*deletaetaX);
        Mal.Ly[i][jj] = (Mal.A[i][jj]*delqsiqsiY -2.0*Mal.B[i][jj]*delqsietaY + Mal.C[i][jj]*deletaetaY);

        Mal.VTX [i] = -r*omega*Mal.Lx[i][jj]  - r*Mal.DELTAx[i][jj-1];
        Mal.VTY [i] = -r*omega*Mal.Ly[i][jj]  - r*Mal.DELTAy[i][jj-1];

            }else{


        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[i-1][jj]);
        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[i-1][jj]);
        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);

        Mal.A[i][jj] = xeta*xeta + yeta*yeta;
        Mal.B[i][jj] = xqsi*xeta + yqsi*yeta;
        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;

        Mal.AT [i-1] = Mal.A[i][jj];
        Mal.BT [i] = -2.0*(1 + Mal.A[i][jj]);
        Mal.CT [i+1] = Mal.A[i][jj];

        delqsiqsiX = Mal.x[i+1][jj] -2.0*Mal.x[i][jj] + Mal.x[i-1][jj];
        delqsietaX = 0.25*(Mal.x[i+1][jj+1] - Mal.x[i+1][jj-1] - Mal.x[i-1][jj+1] + Mal.x[i-1][jj-1]);
        deletaetaX = Mal.x[i][jj+1] -2.0*Mal.x[i][jj] + Mal.x[i][jj-1];

        delqsiqsiY = Mal.y[i+1][jj] -2.0*Mal.y[i][jj] + Mal.y[i-1][jj];
        delqsietaY = 0.25*(Mal.y[i+1][jj+1] - Mal.y[i+1][jj-1] - Mal.y[i-1][jj+1] + Mal.y[i-1][jj-1]);
        deletaetaY = Mal.y[i][jj+1] -2.0*Mal.y[i][jj] + Mal.y[i][jj-1];

        Mal.Lx[i][jj] = (Mal.A[i][jj]*delqsiqsiX -2.0*Mal.B[i][jj]*delqsietaX + Mal.C[i][jj]*deletaetaX);
        Mal.Ly[i][jj] = (Mal.A[i][jj]*delqsiqsiY -2.0*Mal.B[i][jj]*delqsietaY + Mal.C[i][jj]*deletaetaY);

        Mal.VTX [i] = -r*omega*Mal.Lx[i][jj]  - r*Mal.DELTAx[i][jj-1];
        Mal.VTY [i] = -r*omega*Mal.Ly[i][jj]  - r*Mal.DELTAy[i][jj-1];

        }
    }

        Mal.CT [0] = Mal.A[M][jj];

//        Mal.CT [0] =  1;

        TridiagPX(Mal.AT, Mal.BT, Mal.CT, Mal.VTX, xx, M);
        TridiagPY(Mal.AT, Mal.BT, Mal.CT, Mal.VTY, yy, M);

        for (i =0; i<=M; i++){

            Mal.DELTAx[i][jj] = xx[i];
            Mal.DELTAy[i][jj] = yy[i];

            //cout << "DELTAx = " << Mal.DELTAx[i][jj] << endl;
            //cout << "DELTAy = " << Mal.DELTAy[i][jj] << endl;
        }

           //for ( int i = 1; i <= JMAX-2; i++){
           for ( int i = 0; i <= M; i++){
            //cout << "Xold = " << Mal.x[i][jj] << endl;
            //cout << "Yold = " << Mal.y[i][jj] << endl;

            //cout << "DELTAx = " << Mal.DELTAx[i][jj] << endl;
            //cout << "DELTAy = " << Mal.DELTAy[i][jj] << endl;

            Mal.x[i][jj] = Mal.DELTAx[i][jj] + Mal.x[i][jj];// + Mal.x[i][jj] ; // Xi,j(n+1) = Xi,j(n) + DeltaXi,j(n)
            Mal.y[i][jj] = Mal.DELTAy[i][jj] + Mal.y[i][jj];// + Mal.y[i][jj]; // Yi,j(n+1) = Yi,j(n) + DeltaYi,j(n)


            //cout << "X = " << Mal.x[i][jj] << endl;
            //cout << "Y = " << Mal.y[i][jj] << endl;

           // Mal.x[IMAX-1][jj] = Mal.x[0][jj];
            //Mal.y[IMAX-1][jj] = Mal.y[0][jj];

          }

            Mal.x[IMAX-1][jj] = Mal.x[0][jj];
            Mal.y[IMAX-1][jj] = Mal.y[0][jj];
}


    return Mal;

}

Malha AF1 (Malha Mal,int i, int j, int IMAX, int JMAX, double omega, double alfa){

    double xqsi, yqsi, xeta, yeta, delqsiqsiX, delqsietaX, deletaetaX, delqsiqsiY, delqsietaY, deletaetaY;

    double *x; // Armazenando espaço e zerando o vetor solução da periodica
    x = new double [IMAX];

    double *y; // Armazenando espaço e zerando o vetor solução da periodica
    y = new double [IMAX];

    int M = IMAX-2; // Penultimo ponto da periodica
    int n = JMAX-2; // Penultimo ponto da periodica

    int jj;


// Cria e Resolve a tridiagonal Periodica ( 1 Parte do AF1)

for (int jj = 1; jj<=JMAX-2; jj++){

    for (i = 0; i <=M; i++){

             if (i == 0){

//         Esquema centrado
//        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[IMAX-2][jj]);
//        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[IMAX-2][jj]);
//        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
//        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);

        // Esquema do carderno
        xqsi = 0.5*(Mal.x[i+1][jj+1] - Mal.x[IMAX-2][jj+1]);
        yqsi = 0.5*(Mal.y[i+1][jj+1] - Mal.y[IMAX-2][jj+1]);
        xeta = Mal.x[i][jj+1] - Mal.x[i][jj];
        yeta = Mal.y[i][jj+1] - Mal.y[i][jj];

         //cout << xeta << endl;
         //cout << yeta << endl;

        Mal.A[i][jj] = xeta*xeta + yeta*yeta;
        Mal.B[i][jj] = xqsi*xeta + yqsi*yeta;
        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;

        Mal.AT [M] = - Mal.A[i][jj]; // LOW DIAGONAL
        Mal.BT [i] = alfa + 2.0*Mal.A[i][jj]; // MAIN DIAGONAL
        Mal.CT [i+1] = - Mal.A[i][jj]; // UPPER DIAGONAL


        delqsiqsiX= Mal.x[i+1][jj] -2.0*Mal.x[i][jj] + Mal.x[IMAX-2][jj];
        delqsietaX = 0.25*(Mal.x[i+1][jj+1] - Mal.x[i+1][jj-1] - Mal.x[IMAX-2][jj+1] + Mal.x[IMAX-2][jj-1]);
        deletaetaX = Mal.x[i][jj+1] -2.0*Mal.x[i][jj] + Mal.x[i][jj-1];

        delqsiqsiY= Mal.y[i+1][jj] -2.0*Mal.y[i][jj] + Mal.y[IMAX-2][jj];
        delqsietaY = 0.25*(Mal.y[i+1][jj+1] - Mal.y[i+1][jj-1] - Mal.y[IMAX-2][jj+1] + Mal.y[IMAX-2][jj-1]);
        deletaetaY = Mal.y[i][jj+1] -2.0*Mal.y[i][jj] + Mal.y[i][jj-1];

        // Multiplicando B por 0.5
        Mal.Lx[i][jj] = (Mal.A[i][jj]*delqsiqsiX -2.0*Mal.B[i][jj]*delqsietaX + Mal.C[i][jj]*deletaetaX);
        Mal.Ly[i][jj] = (Mal.A[i][jj]*delqsiqsiY -2.0*Mal.B[i][jj]*delqsietaY + Mal.C[i][jj]*deletaetaY);


        Mal.VTX [i] = alfa*omega*Mal.Lx[i][jj];
        Mal.VTY [i] = alfa*omega*Mal.Ly[i][jj];

            }else{
////         Esquema centrado sugerido parao problema
//        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[i-1][jj]);
//        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[i-1][jj]);
//        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
//        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);
//
        // // Esquema do carderno
        xqsi = 0.5*(Mal.x[i+1][jj+1] - Mal.x[i-1][jj+1]);
        yqsi = 0.5*(Mal.y[i+1][jj+1] - Mal.y[i-1][jj+1]);
        xeta = Mal.x[i][jj+1] - Mal.x[i][jj];
        yeta = Mal.y[i][jj+1] - Mal.y[i][jj];

        Mal.A[i][jj] = xeta*xeta + yeta*yeta;
        Mal.B[i][jj] = xqsi*xeta + yqsi*yeta;
        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;

        Mal.AT [i-1] = - Mal.A[i][jj]; // LOW DIAGONAL
        Mal.BT [i] =   alfa + 2.0*Mal.A[i][jj]; // MAIN DIAGONAL
        Mal.CT [i+1] = - Mal.A[i][jj]; // UPPER DIAGONAL

        delqsiqsiX = Mal.x[i+1][jj] -2.0*Mal.x[i][jj] + Mal.x[i-1][jj];
        delqsietaX = 0.25*(Mal.x[i+1][jj+1] - Mal.x[i+1][jj-1] - Mal.x[i-1][jj+1] + Mal.x[i-1][jj-1]);
        deletaetaX = Mal.x[i][jj+1] -2.0*Mal.x[i][jj] + Mal.x[i][jj-1];

        delqsiqsiY = Mal.y[i+1][jj] -2.0*Mal.y[i][jj] + Mal.y[i-1][jj];
        delqsietaY = 0.25*(Mal.y[i+1][jj+1] - Mal.y[i+1][jj-1] - Mal.y[i-1][jj+1] + Mal.y[i-1][jj-1]);
        deletaetaY = Mal.y[i][jj+1] -2.0*Mal.y[i][jj] + Mal.y[i][jj-1];

        // Multiplicando B por 0.5
        Mal.Lx[i][jj] = (Mal.A[i][jj]*delqsiqsiX -2.0*Mal.B[i][jj]*delqsietaX + Mal.C[i][jj]*deletaetaX);
        Mal.Ly[i][jj] = (Mal.A[i][jj]*delqsiqsiY -2.0*Mal.B[i][jj]*delqsietaY + Mal.C[i][jj]*deletaetaY);

        Mal.VTX [i] = alfa*omega*Mal.Lx[i][jj];
        Mal.VTY [i] = alfa*omega*Mal.Ly[i][jj];

        }
    }

        Mal.CT [0] = - Mal.A[M][jj];

        TridiagPX(Mal.AT, Mal.BT, Mal.CT, Mal.VTX, x, M);
        TridiagPY(Mal.AT, Mal.BT, Mal.CT, Mal.VTY, y, M);

        for (i =0; i<=M; i++){

            Mal.f[i][jj] = x[i];
            Mal.g[i][jj] = y[i];

        }

        //Mal.f[IMAX-1][jj] = Mal.f[0][jj]; Garanto que é igual
        //Mal.g[IMAX-1][jj] = Mal.g[0][jj];
}

    // Resolve o Segundo Sweep do Metodo resolvendo agora uma tridiagonal simples pelo método TDMA
for (int i = 0; i<=IMAX-2; i++){

    for (int jj = 1; jj <=n; jj++){


        if (i == 0){

//        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[IMAX-2][jj]);
//        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[IMAX-2][jj]);

                // Esquema do carderno
        xqsi = 0.5*(Mal.x[i+1][jj+1] - Mal.x[IMAX-2][jj+1]);
        yqsi = 0.5*(Mal.y[i+1][jj+1] - Mal.y[IMAX-2][jj+1]);

        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;


        }else{

//        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[i-1][jj]);
//        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[i-1][jj]);

                // // Esquema do carderno
        xqsi = 0.5*(Mal.x[i+1][jj+1] - Mal.x[i-1][jj+1]);
        yqsi = 0.5*(Mal.y[i+1][jj+1] - Mal.y[i-1][jj+1]);

        }

//        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
//        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);

        // Caderno
        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj]);
        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj]);

        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;


        Mal.AT [jj-1] = -Mal.C[i][jj]; // LOW DIAGONAL
        Mal.BT [jj-1] = alfa + 2.0*Mal.C[i][jj] ; // MAIN DIAGONAL
        Mal.CT [jj-1] = -Mal.C[i][jj] ; // UPPER DIAGONAL

        Mal.VTX [jj-1] =  Mal.f[i][jj] ;
        Mal.VTY [jj-1] =  Mal.g[i][jj] ;

       // }

        if ( jj == 1){
         Mal.AT[0] = 0.0;
        }

        if (jj == n){
        Mal.CT[n] = 0.0;
        }

        //cout << "AT = " << Mal.AT[jj-1] << endl;
       // cout << "BT = " << Mal.BT[jj-1]<< endl;
       // cout << "CT = " << Mal.CT[jj-1] << endl;

    }

        TDMAX(Mal.AT, Mal.BT, Mal.CT, Mal.VTX, x, n);
        TDMAY(Mal.AT, Mal.BT, Mal.CT, Mal.VTY, y, n);


          for ( int jj = 1; jj <= n; jj++){

            Mal.DELTAx[i][jj] = x[jj-1];// + Mal.x[i][jj] ; // Xi,j(n+1) = Xi,j(n) + DeltaXi,j(n)
            Mal.DELTAy[i][jj] = y[jj-1];// + Mal.y[i][jj]; // Yi,j(n+1) = Yi,j(n) + DeltaYi,j(n)

           }

        // Correção dos Pontos para todos os pontos j da matrix

         for ( int jj = 1; jj <= JMAX-2; jj++){

            Mal.x[i][jj] = Mal.DELTAx[i][jj] + Mal.x[i][jj];// + Mal.x[i][jj] ; // Xi,j(n+1) = Xi,j(n) + DeltaXi,j(n)
            Mal.y[i][jj] = Mal.DELTAy[i][jj] + Mal.y[i][jj];// + Mal.y[i][jj]; // Yi,j(n+1) = Yi,j(n) + DeltaYi,j(n)

            Mal.x[IMAX-1][jj] = Mal.x[0][jj];
            Mal.y[IMAX-1][jj] = Mal.y[0][jj];

          }

}

return Mal;


}

Malha AF2 (Malha Mal,int i, int j, int IMAX, int JMAX, double omega, double alfa){

    double xqsi, yqsi, xeta, yeta, delqsiqsiX, delqsietaX, deletaetaX, delqsiqsiY, delqsietaY, deletaetaY;

    double *x; // Armazenando espaço e zerando o vetor solução da periodica
    x = new double [IMAX];

    double *y; // Armazenando espaço e zerando o vetor solução da periodica
    y = new double [IMAX];

    int M = IMAX-2; // Penultimo ponto da periodica
    int n = JMAX-2; // Penultimo ponto da periodica

    int jj;

// Cria e Resolve a tridiagonal Periodica ( 1 Parte do AF2)

for (int jj = 1; jj<=JMAX-2; jj++){

    for (i = 0; i <=M; i++){

             if (i == 0){

//         Esquema centrado
//        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[IMAX-2][jj]);
//        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[IMAX-2][jj]);
//        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
//        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);

        // Esquema do carderno
        xqsi = 0.5*(Mal.x[i+1][jj+1] - Mal.x[IMAX-2][jj+1]);
        yqsi = 0.5*(Mal.y[i+1][jj+1] - Mal.y[IMAX-2][jj+1]);
        xeta = Mal.x[i][jj+1] - Mal.x[i][jj];
        yeta = Mal.y[i][jj+1] - Mal.y[i][jj];

        Mal.A[i][jj] = xeta*xeta + yeta*yeta;
        Mal.B[i][jj] = xqsi*xeta + yqsi*yeta;
        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;

        Mal.AT [M] = - Mal.A[i][jj]; // LOW DIAGONAL
        Mal.BT [i] = alfa + 2.0*Mal.A[i][jj]; // MAIN DIAGONAL
        Mal.CT [i+1] = - Mal.A[i][jj]; // UPPER DIAGONAL

        delqsiqsiX= Mal.x[i+1][jj] -2.0*Mal.x[i][jj] + Mal.x[IMAX-2][jj];
        delqsietaX = 0.25*(Mal.x[i+1][jj+1] - Mal.x[i+1][jj-1] - Mal.x[IMAX-2][jj+1] + Mal.x[IMAX-2][jj-1]);
        deletaetaX = Mal.x[i][jj+1] -2.0*Mal.x[i][jj] + Mal.x[i][jj-1];

        delqsiqsiY= Mal.y[i+1][jj] -2.0*Mal.y[i][jj] + Mal.y[IMAX-2][jj];
        delqsietaY = 0.25*(Mal.y[i+1][jj+1] - Mal.y[i+1][jj-1] - Mal.y[IMAX-2][jj+1] + Mal.y[IMAX-2][jj-1]);
        deletaetaY = Mal.y[i][jj+1] -2.0*Mal.y[i][jj] + Mal.y[i][jj-1];

//        Mal.Lx[i][jj] = (Mal.A[i][jj]*delqsiqsiX -2.0*Mal.B[i][jj]*delqsietaX + Mal.C[i][jj]*deletaetaX);
//        Mal.Ly[i][jj] = (Mal.A[i][jj]*delqsiqsiY -2.0*Mal.B[i][jj]*delqsietaY + Mal.C[i][jj]*deletaetaY);

        // Multiplicando B por 0.5
        Mal.Lx[i][jj] = (Mal.A[i][jj]*delqsiqsiX -2.0*Mal.B[i][jj]*delqsietaX + Mal.C[i][jj]*deletaetaX);
        Mal.Ly[i][jj] = (Mal.A[i][jj]*delqsiqsiY -2.0*Mal.B[i][jj]*delqsietaY + Mal.C[i][jj]*deletaetaY);

        Mal.VTX [i] = alfa*omega*Mal.Lx[i][jj] + alfa*Mal.f[i][jj-1];
        Mal.VTY [i] = alfa*omega*Mal.Ly[i][jj] + alfa*Mal.g[i][jj-1];

        //cout << "VTX = " << Mal.VTX [i] << endl;
        //cout << "VTY = " << Mal.VTY [i] << endl;
        //cout << "f = " << Mal.f[i][jj-1] << endl;
        //cout << "g = " << Mal.g[i][jj-1] << endl;
            }else{
////         Esquema centrado sugerido parao problema
//        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[i-1][jj]);
//        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[i-1][jj]);
//        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
//        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);
//
        // // Esquema do carderno
        xqsi = 0.5*(Mal.x[i+1][jj+1] - Mal.x[i-1][jj+1]);
        yqsi = 0.5*(Mal.y[i+1][jj+1] - Mal.y[i-1][jj+1]);
        xeta = Mal.x[i][jj+1] - Mal.x[i][jj];
        yeta = Mal.y[i][jj+1] - Mal.y[i][jj];

        Mal.A[i][jj] = xeta*xeta + yeta*yeta;
        Mal.B[i][jj] = xqsi*xeta + yqsi*yeta;
        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;

//        Mal.AT [i-1] = alfa - Mal.A[i][jj]; // LOW DIAGONAL
//        Mal.BT [i] =   alfa + 2.0*Mal.A[i][jj]; // MAIN DIAGONAL
//        Mal.CT [i+1] = alfa - Mal.A[i][jj]; // UPPER DIAGONAL

        Mal.AT [i-1] = - Mal.A[i][jj]; // LOW DIAGONAL
        Mal.BT [i] =   alfa + 2.0*Mal.A[i][jj]; // MAIN DIAGONAL
        Mal.CT [i+1] = - Mal.A[i][jj]; // UPPER DIAGONAL

       //cout << "AT = " << Mal.AT[i-1] << endl;
        //cout << "BT = " << Mal.BT[i]<< endl;
        //cout << "CT = " << Mal.CT[i+1] << endl;

        delqsiqsiX = Mal.x[i+1][jj] -2.0*Mal.x[i][jj] + Mal.x[i-1][jj];
        delqsietaX = 0.25*(Mal.x[i+1][jj+1] - Mal.x[i+1][jj-1] - Mal.x[i-1][jj+1] + Mal.x[i-1][jj-1]);
        deletaetaX = Mal.x[i][jj+1] -2.0*Mal.x[i][jj] + Mal.x[i][jj-1];

        delqsiqsiY = Mal.y[i+1][jj] -2.0*Mal.y[i][jj] + Mal.y[i-1][jj];
        delqsietaY = 0.25*(Mal.y[i+1][jj+1] - Mal.y[i+1][jj-1] - Mal.y[i-1][jj+1] + Mal.y[i-1][jj-1]);
        deletaetaY = Mal.y[i][jj+1] -2.0*Mal.y[i][jj] + Mal.y[i][jj-1];

//        Mal.Lx[i][jj] = (Mal.A[i][jj]*delqsiqsiX -2.0*Mal.B[i][jj]*delqsietaX + Mal.C[i][jj]*deletaetaX);
//        Mal.Ly[i][jj] = (Mal.A[i][jj]*delqsiqsiY -2.0*Mal.B[i][jj]*delqsietaY + Mal.C[i][jj]*deletaetaY);

        // Multiplicando B por 0.5
        Mal.Lx[i][jj] = (Mal.A[i][jj]*delqsiqsiX -2.0*Mal.B[i][jj]*delqsietaX + Mal.C[i][jj]*deletaetaX);
        Mal.Ly[i][jj] = (Mal.A[i][jj]*delqsiqsiY -2.0*Mal.B[i][jj]*delqsietaY + Mal.C[i][jj]*deletaetaY);

        //cout << "Lx = " << Mal.Lx[i][jj]<< endl;
        //cout << "Ly = " << Mal.Ly[i][jj] << endl;

        Mal.VTX [i] = alfa*omega*Mal.Lx[i][jj] + alfa*Mal.f[i][jj-1];
        Mal.VTY [i] = alfa*omega*Mal.Ly[i][jj] + alfa*Mal.g[i][jj-1];

        //cout << "VTX = " << Mal.VTX [i] << endl;
        //cout << "VTY = " << Mal.VTY [i] << endl;

        }
    }

        //Mal.CT [0] = alfa - Mal.A[M][jj];

        Mal.CT [0] = - Mal.A[M][jj];
        //cout << Mal.CT [0] << endl;
        //Vetor::imprimir(Mal.CT,IMAX);

        TridiagPX(Mal.AT, Mal.BT, Mal.CT, Mal.VTX, x, M);
        TridiagPY(Mal.AT, Mal.BT, Mal.CT, Mal.VTY, y, M);

        for (i =0; i<=M; i++){

            Mal.f[i][jj] = x[i];
            Mal.g[i][jj] = y[i];

        }

        //Mal.f[IMAX-1][jj] = Mal.f[0][jj]; Garanto que é igual
        //Mal.g[IMAX-1][jj] = Mal.g[0][jj];
}

// Resolve o Segundo Sweep do Metodo resolvendo agora uma tridiagonal simples pelo método TDMA
for (int i = 0; i<=IMAX-2; i++){

    for (int jj = 1; jj <=n; jj++){

        if (i == 0){

//        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[IMAX-2][jj]);
//        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[IMAX-2][jj]);

                // Esquema do carderno
        xqsi = 0.5*(Mal.x[i+1][jj+1] - Mal.x[IMAX-2][jj+1]);
        yqsi = 0.5*(Mal.y[i+1][jj+1] - Mal.y[IMAX-2][jj+1]);

        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;


        }else{

//        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[i-1][jj]);
//        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[i-1][jj]);

                // // Esquema do carderno
        xqsi = 0.5*(Mal.x[i+1][jj+1] - Mal.x[i-1][jj+1]);
        yqsi = 0.5*(Mal.y[i+1][jj+1] - Mal.y[i-1][jj+1]);

        }

//        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
//        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);

        // Caderno
        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj]);
        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj]);


        //Mal.A[i][jj] = xeta*xeta + yeta*yeta;
        //Mal.B[i][jj] = xqsi*xeta + yqsi*yeta;
        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;

        //cout << Mal.C[i][jj] << endl;

//        Mal.AT [jj-1] = alfa -Mal.C[i][jj]; // LOW DIAGONAL
//        Mal.BT [jj-1] = alfa + 2.0*Mal.C[i][jj] ; // MAIN DIAGONAL
//        Mal.CT [jj-1] = alfa -Mal.C[i][jj] ; // UPPER DIAGONAL


        Mal.AT [jj-1] = 0.0; // LOW DIAGONAL
        Mal.BT [jj-1] = alfa + Mal.C[i][jj] ; // MAIN DIAGONAL
        Mal.CT [jj-1] = -Mal.C[i][jj] ; // UPPER DIAGONAL

       // cout << "AT = " << Mal.AT[jj-1] << endl;
        //cout << "BT = " << Mal.BT[jj-1]<< endl;
        //cout << "CT = " << Mal.CT[jj-1] << endl;

        Mal.VTX [jj-1] =  Mal.f[i][jj] ;
        Mal.VTY [jj-1] =  Mal.g[i][jj] ;

        //cout << "VTX = " << Mal.VTX [jj-1]<< endl;
        //cout << "VTY = " << Mal.VTY [jj-1]<< endl;

       // }

        if ( jj == 1){
         Mal.AT[0] = 0.0;
        }

        if (jj == n){
        Mal.CT[n] = 0.0;
        }

        //cout << "AT = " << Mal.AT[jj-1] << endl;
       // cout << "BT = " << Mal.BT[jj-1]<< endl;
       // cout << "CT = " << Mal.CT[jj-1] << endl;

    }

        TDMAX(Mal.AT, Mal.BT, Mal.CT, Mal.VTX, x, n);
        TDMAY(Mal.AT, Mal.BT, Mal.CT, Mal.VTY, y, n);


          for ( int jj = 1; jj <= n; jj++){

            Mal.DELTAx[i][jj] = x[jj-1];// + Mal.x[i][jj] ; // Xi,j(n+1) = Xi,j(n) + DeltaXi,j(n)
            Mal.DELTAy[i][jj] = y[jj-1];// + Mal.y[i][jj]; // Yi,j(n+1) = Yi,j(n) + DeltaYi,j(n)

           // cout << "DELTAx = " << Mal.DELTAx[i][jj] << endl;
           // cout << "DELTAy = " << Mal.DELTAy[i][jj] << endl;

            }

        // Correção dos Pontos para todos os pontos j da matrix

         for ( int jj = 1; jj <= JMAX-2; jj++){

            //cout << "Xold = " << Mal.x[i][jj] << endl;
            //cout << "Yold = " << Mal.y[i][jj] << endl;

            Mal.x[i][jj] = Mal.DELTAx[i][jj] + Mal.x[i][jj];// + Mal.x[i][jj] ; // Xi,j(n+1) = Xi,j(n) + DeltaXi,j(n)
            Mal.y[i][jj] = Mal.DELTAy[i][jj] + Mal.y[i][jj];// + Mal.y[i][jj]; // Yi,j(n+1) = Yi,j(n) + DeltaYi,j(n)

           // cout << "X = " << Mal.x[i][jj] << endl;
           // cout << "Y = " << Mal.y[i][jj] << endl;

            Mal.x[IMAX-1][jj] = Mal.x[0][jj];
            Mal.y[IMAX-1][jj] = Mal.y[0][jj];

          }

}

return Mal;


}

void TridiagPX(double* a, double* b, double* c, double* dv, double *x, int M){

   // Essa implementação foi feita baseando no trabalho de Doutorado do Prof. Azevedo no anexo sobre a solução tridiagonal Periodica

    int i, IMAX=93;

    double *l1,*l2,*l3,*u1,*u2,*u3,*f,*g;
    double soma1, soma2;

    l1 = new double [M];
    l2 = new double [M];
    l3 = new double [M];

    u1 = new double [M];
    u2 = new double [M];
    u3 = new double [M];

    f = new double [M];
    g = new double [M];

    for (i =0; i<=M;i++){ // Uma vez que foi alocado dinamicamente (via pointer) é conveniente zerar as posições alocadas, evitando pegar lixo de máquina

        l1[i] = 0.0;
        l2[i] = 0.0;
        l3[i] = 0.0;
        u1[i] = 0.0;
        u2[i] = 0.0;
        u3[i] = 0.0;
        f[i]  = 0.0;
        g[i]  = dv[i];
    }


    l1[0] = b[0]; // L1 = B1
    l2[0] = a[0]; // M1 = A1
    l3[0] = c[0]; // N1 = Cn+1
    f[0]  = g[0]/l1[0]; // x_barra1 = L1^-1*y1
    u2[1] = a[M]/l1[0];  // V1 = L1^-1*A0 igual doutorado do prof Azevedo

    soma1 = l3[0]*u2[1]; // Ta = N1V1
    soma2 = l3[0]*f[0];  //

         /////////////////////////
        // Forward Sweep //////
        ////////////////////////


    for (i = 1; i <=(M-2); i++){

        u1[i] = c[i]/l1[i-1];
        l1[i] = b[i] - u1[i]*l2[i-1];
        l2[i] = a[i];
        l3[i] = -l3[i-1]*u1[i];
        f[i]  = (g[i]-l2[i-1]*f[i-1])/l1[i];
        u2[i+1] = -l2[i-1]*u2[i]/l1[i];
        soma1 = soma1 + l3[i]*u2[i+1];
        soma2 = soma2 + l3[i]*f[i];
    }

        u1[M-1] = c[M-1]/l1[M-2];
        l1[M-1] = b[M-1] - l2[M-2]*u1[M-1];
        l2[M-1] = a[M-1] - l3[M-2]*u1[M-1];
         f[M-1] = (g[M-1]-l2[M-2]*f[M-2])/l1[M-1];

         u1[M] = (c[M]-l2[M-2]*u2[M-1])/l1[M-1];
         soma1 = soma1 + l2[M-1]*u1[M];
         l1[M] = b[M]-soma1;
         f[M] = (g[M] -l2[M-1]*f[M-1]-soma2)/l1[M];

        /////////////////////////
        // Backward Sweep //////
        ////////////////////////

        x[M] = f[M];

        x[M-1] = f[M-1] - u1[M]*x[M];

        for (i = (M-2); i>=0; i--){
            x[i] = f[i] - u1[i+1]*x[i+1] - u2[i+1]*x[M];

        }

 return;

}

void TridiagPY(double* a, double* b, double* c, double* dv, double *y, int M){

   // Essa implementação foi feita baseando no trabalho de Doutorado do Prof. Azevedo no anexo sobre a solução tridiagonal Periodica
int i, IMAX=93;

    double *l1,*l2,*l3,*u1,*u2,*u3,*f,*g;
    double soma1, soma2;

    l1 = new double [M];
    l2 = new double [M];
    l3 = new double [M];

    u1 = new double [M];
    u2 = new double [M];
    u3 = new double [M];

    f = new double [M];
    g = new double [M];

    for (i =0; i<=M;i++){ // Uma vez que foi alocado dinamicamente (via pointer) é conveniente zerar as posições alocadas, evitando pegar lixo de máquina

        l1[i] = 0.0;
        l2[i] = 0.0;
        l3[i] = 0.0;
        u1[i] = 0.0;
        u2[i] = 0.0;
        u3[i] = 0.0;
        f[i]  = 0.0;
        g[i]  = dv[i];
    }

    //cout << dv[0] << endl;
    //cout << dv[1] << endl;


    //Vetor::imprimir(g,IMAX);

    l1[0] = b[0]; // L1 = B1
    //cout << b[0] << endl;

    l2[0] = a[0]; // M1 = A1
   // cout << a[0] << endl;

    l3[0] = c[0]; // N1 = Cn+1
    //cout << c[0] << endl;

    f[0]  = g[0]/l1[0]; // x_barra1 = L1^-1*y1
    //cout << f[0] << endl;

    u2[1] = a[M]/l1[0];  // V1 = L1^-1*A0 igual doutorado do prof Azevedo
    //cout << u2[1] << endl;

    soma1 = l3[0]*u2[1]; // Ta = N1V1
    //cout << soma1 << endl;

    soma2 = l3[0]*f[0];  //
    //cout << soma2 << endl;

         /////////////////////////
        // Forward Sweep //////
        ////////////////////////

        //Vetor::imprimir(a,IMAX);
        //Vetor::imprimir(b,IMAX);
        //Vetor::imprimir(c,IMAX);
        //Vetor::imprimir(dv,IMAX);

    for (i = 1; i <=(M-2); i++){

        u1[i] = c[i]/l1[i-1];

        //cout <<c[i] << endl;
        //cout <<l1[i-1] << endl;
        //cout <<u1[i] << endl;

        l1[i] = b[i] - u1[i]*l2[i-1];
        l2[i] = a[i];
        l3[i] = -l3[i-1]*u1[i];
        f[i]  = (g[i]-l2[i-1]*f[i-1])/l1[i];
        u2[i+1] = -l2[i-1]*u2[i]/l1[i];
        soma1 = soma1 + l3[i]*u2[i+1];
        soma2 = soma2 + l3[i]*f[i];
    }

        u1[M-1] = c[M-1]/l1[M-2];
        l1[M-1] = b[M-1] - l2[M-2]*u1[M-1];
        l2[M-1] = a[M-1] - l3[M-2]*u1[M-1];
         f[M-1] = (g[M-1]-l2[M-2]*f[M-2])/l1[M-1];

         u1[M] = (c[M]-l2[M-2]*u2[M-1])/l1[M-1];
         soma1 = soma1 + l2[M-1]*u1[M];
         l1[M] = b[M]-soma1;
         f[M] = (g[M] -l2[M-1]*f[M-1]-soma2)/l1[M];

        /////////////////////////
        // Backward Sweep //////
        ////////////////////////

        y[M] = f[M];

       // cout << y[M] << endl;
        y[M-1] = f[M-1] - u1[M]*y[M];
        //cout << y[M-1] << endl;

        for (i = (M-2); i>=0; i--){
            y[i] = f[i] - u1[i+1]*y[i+1] - u2[i+1]*y[M];
            // cout << y[i] << endl;

            //cout << x[i] << endl;
        }



 return;

}

void TDMAX(double* a, double* b, double* c, double* d, double *x, int n){

    int k;
    double m;

    for ( k = 1 ; k <= n-1; k++){

    m = a[k]/b[k-1];
    b[k] = b[k] -m*c[k-1];
    d[k] = d[k] -m*d[k-1];

    }

    x[n-1] = d[n-1]/b[n-1];

    for (k = n-1; k >= 0; k--){

        x[k] = (d[k] - c[k]*x[k+1])/b[k];
    }

return;

}

void TDMAY(double* a, double* b, double* c, double* d, double *y, int n){

    int k;
    double m;

    for ( k = 1 ; k <= n-1; k++){

    m = a[k]/b[k-1];
    b[k] = b[k] -m*c[k-1];
    d[k] = d[k] -m*d[k-1];

    }

    y[n-1] = d[n-1]/b[n-1];

    for (k = n-1; k >= 0; k--){

        y[k] = (d[k] - c[k]*y[k+1])/b[k];
    }

return;

}

double MaxResiduoX (Malha Mal,int IMAX,int JMAX,int kiter, double resmaxX){

     //resmaxX = 0.0;

      for (int i = 1; i < IMAX-1; i++){

        for (int j = 1; j < JMAX-1; j++ ){

            if (fabs(Mal.Lx[i][j])> resmaxX ){
              resmaxX = fabs(Mal.Lx[i][j]);
            }
        }
      }

//cout << resmaxX << endl;


return resmaxX;

}

double MaxResiduoY (Malha Mal,int IMAX,int JMAX,int kiter,double resmaxY){

      //resmaxY = 0.0;

      for (int i = 1; i < IMAX-1; i++){

        for (int j = 1; j < JMAX-1; j++ ){

            if (fabs(Mal.Ly[i][j])> resmaxY ){
              resmaxY = fabs(Mal.Ly[i][j]);
            }
        }
      }


//cout << resmaxY << endl;

return resmaxY;

}

void AlfaSequence (int Melements, double alfaLow, double alfaHigh){

 double *alfa, expoente;
    int k, kaux;

    alfa = new double [Melements];

    for (k = 0; k<Melements;k++){

        kaux = k+1;

        if (k == Melements){
            kaux = k;
        }

        expoente = (kaux -1.0)/(Melements-1.0);
        alfa[k] = alfaHigh*pow((alfaLow/alfaHigh),expoente);

        //cout << kaux <<endl;
        cout << expoente <<endl;
        cout << alfa[k] <<endl;
    }

return ;

return;
}

Malha aerofolioPontosNACA (Malha Mal, int i, int j, int IMAX, int JMAX,double R,double th, double XSF){

    // Gerando o circulo exterior
    for ( i = 0; i <(IMAX-1); i++){ //  IMAX-1 = 92  >>> vai de 0 - 91

    Mal.x[i][JMAX-1] = R*cos((2*M_PI*(double(i)/(IMAX-1)))); //(IMAX-2/IMAX-1)
    Mal.y[i][JMAX-1] = -R*sin((2*M_PI*(double(i)/(IMAX-1)))); // coloquei menos para ser sentido horario

    }

    Mal.x[IMAX-1][JMAX-1] = Mal.x[0][JMAX-1];
    Mal.y[IMAX-1][JMAX-1] = Mal.y[0][JMAX-1];

 //Implementação com Estiramento usando Progressão Geometrica (esta sera apresentada no relatorio)

        int numPointsPG = ((IMAX-1)/4);
        double somaPG = 0.5;
        double dx = somaPG*(XSF-1)/(pow(XSF,numPointsPG)-1.0);

        Mal.x[0][0] = 0.0;
        Mal.x[1][0] = Mal.x[0][0] + dx;

        // Perfil NACA 0012
        Mal.y[0][0] =  (th/0.2)*(0.2969*pow(Mal.x[0][0],0.5) - 0.1260*Mal.x[0][0] - 0.3516*pow(Mal.x[0][0],2)+ 0.2843*pow(Mal.x[0][0],3) -0.1036*pow(Mal.x[0][0],4));
        Mal.y[1][0] = -(th/0.2)*(0.2969*pow(Mal.x[1][0],0.5) - 0.1260*Mal.x[1][0] - 0.3516*pow(Mal.x[1][0],2)+ 0.2843*pow(Mal.x[1][0],3) -0.1036*pow(Mal.x[1][0],4));;

        Mal.x[(IMAX-1)/2][0] = 1.0;
        Mal.x[((IMAX-1)/2) -1][0] = Mal.x[(IMAX-1)/2][0] - dx;

        Mal.y[(IMAX-1)/2][0] =    (th/0.2)*(0.2969*pow(Mal.x[(IMAX-1)/2][0],0.5) - 0.1260*Mal.x[(IMAX-1)/2][0] - 0.3516*pow(Mal.x[(IMAX-1)/2][0],2) + 0.2843*pow(Mal.x[(IMAX-1)/2][0],3)-0.1036*pow(Mal.x[(IMAX-1)/2][0],4));
        Mal.y[(IMAX-1)/2-1][0] = -(th/0.2)*(0.2969*pow(Mal.x[(IMAX-1)/2-1][0],0.5) - 0.1260*Mal.x[(IMAX-1)/2-1][0] - 0.3516*pow(Mal.x[(IMAX-1)/2-1][0],2) + 0.2843*pow(Mal.x[(IMAX-1)/2-1][0],3) -0.1036*pow(Mal.x[(IMAX-1)/2-1][0],4));

    // Gerando para o primeiro quadrante
    for ( i = 2; i <=((IMAX-1)/4); i++){ // OK

        Mal.x[i][0] = (Mal.x[i-1][0] + XSF*(Mal.x[i-1][0] - Mal.x[i-2][0]));
        Mal.y[i][0] = -(th/0.2)*(0.2969*pow(Mal.x[i][0],0.5) - 0.1260*Mal.x[i][0] - 0.3516*Mal.x[i][0]*Mal.x[i][0] + 0.2843*Mal.x[i][0]*Mal.x[i][0]*Mal.x[i][0] -0.1036*Mal.x[i][0]*Mal.x[i][0]*Mal.x[i][0]*Mal.x[i][0]);
    }

    // Gerando o segundo quadrante
    for ( i =((IMAX-1)/2 -2); i >=((IMAX-1)/4 +1) ; i-- ){

        Mal.x[i][0] = Mal.x[i+1][0] - (Mal.x[i+2][0] - Mal.x[i+1][0])*XSF;
        Mal.y[i][0]= -(th/0.2)*(0.2969*pow(Mal.x[i][0],0.5) - 0.1260*Mal.x[i][0] - 0.3516*pow(Mal.x[i][0],2) + 0.2843*pow(Mal.x[i][0],3) -0.1036*pow(Mal.x[i][0],4));

      }

    for (i = 0; i<=((IMAX-1)/2);i++){ // ESSE LAÇO GERA O X estirado certo de 0 até (imax-1)/2

        Mal.x[i][0] = Mal.x[i][0] - 0.5; // ERA : Mal.x[i][0] - 0.5
        Mal.x[i][0] = -Mal.x[i][0]; //era -Mal.x[i][0] $$$$$$$$$$$$$$$$$$$$$$$$$$
    }

    int k = ((IMAX-1)/2);
    // Gerando os valores de y para o terceiro e quarto quadrantes
    for (i = (((IMAX-1)/2)); i <=IMAX-2; i++){ // era i <=IMAX-1 e agora é : i <=IMAX-2

        Mal.y[i][0] = -Mal.y[k][0]; // ERA : - Mal.y[k][0]
        Mal.x[i][0] = Mal.x[k][0];
        k = k-1;
    }

    int kk = 0; // 92

    for (i = (((IMAX-1)/2)); i <=IMAX-1; i++){ // 46 até 92

        Mal.y[i][0] = Mal.y[kk][0];
        //cout << Mal.y[kk][0] << endl;
        kk = kk+1;
        //cout << i << endl;

    }

    int kkk = (IMAX-1);

  for (i = 0; i <=(((IMAX-1)/2)); i++){ // 46 até 92

        Mal.y[i][0] =Mal.y[kkk][0];
        kkk = kkk-1;
    }

     k = ((IMAX-1)/2);

    for (i =((IMAX-1)/2); i<=IMAX-1 ; i++ ){

    Mal.y[i][0] = -Mal.y[k][0];
    k = k-1;

    }

    Mal.x[IMAX-1][0] = Mal.x[0][0];
    Mal.y[IMAX-1][0] = Mal.y[0][0];

return Mal;

}

Malha AF2Dirichlet (Malha Mal,int i, int j, int IMAX, int JMAX, double omega, double alfa){

    double xqsi, yqsi, xeta, yeta, delqsiqsiX, delqsietaX, deletaetaX, delqsiqsiY, delqsietaY, deletaetaY;

    double *x; // Armazenando espaço e zerando o vetor solução da periodica
    x = new double [IMAX];

    double *y; // Armazenando espaço e zerando o vetor solução da periodica
    y = new double [IMAX];

    int M = IMAX-2; // Penultimo ponto da periodica
    int n = JMAX-2; // Penultimo ponto da periodica

    int jj;

// Cria e Resolve a tridiagonal Periodica ( 1 Parte do AF2)

for (int jj = 1; jj<=JMAX-2; jj++){

    for (i = 0; i <=M; i++){

             if (i == 0){

//         Esquema centrado
//        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[IMAX-2][jj]);
//        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[IMAX-2][jj]);
//        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
//        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);

        // Esquema do carderno
        xqsi = 0.5*(Mal.x[i+1][jj+1] - Mal.x[IMAX-2][jj+1]);
        yqsi = 0.5*(Mal.y[i+1][jj+1] - Mal.y[IMAX-2][jj+1]);
        xeta = Mal.x[i][jj+1] - Mal.x[i][jj];
        yeta = Mal.y[i][jj+1] - Mal.y[i][jj];

        Mal.A[i][jj] = xeta*xeta + yeta*yeta;
        Mal.B[i][jj] = xqsi*xeta + yqsi*yeta;
        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;

        Mal.AT [M] = 0.0; // LOW DIAGONAL
        Mal.BT [i] = alfa + 2.0*Mal.A[i][jj]; // MAIN DIAGONAL
        Mal.CT [i+1] = - Mal.A[i][jj]; // UPPER DIAGONAL

        delqsiqsiX= Mal.x[i+1][jj] -2.0*Mal.x[i][jj] + Mal.x[IMAX-2][jj];
        delqsietaX = 0.25*(Mal.x[i+1][jj+1] - Mal.x[i+1][jj-1] - Mal.x[IMAX-2][jj+1] + Mal.x[IMAX-2][jj-1]);
        deletaetaX = Mal.x[i][jj+1] -2.0*Mal.x[i][jj] + Mal.x[i][jj-1];

        delqsiqsiY= Mal.y[i+1][jj] -2.0*Mal.y[i][jj] + Mal.y[IMAX-2][jj];
        delqsietaY = 0.25*(Mal.y[i+1][jj+1] - Mal.y[i+1][jj-1] - Mal.y[IMAX-2][jj+1] + Mal.y[IMAX-2][jj-1]);
        deletaetaY = Mal.y[i][jj+1] -2.0*Mal.y[i][jj] + Mal.y[i][jj-1];


        // Multiplicando B por 0.5
        Mal.Lx[i][jj] = (Mal.A[i][jj]*delqsiqsiX -2.0*Mal.B[i][jj]*delqsietaX + Mal.C[i][jj]*deletaetaX);
        Mal.Ly[i][jj] = (Mal.A[i][jj]*delqsiqsiY -2.0*Mal.B[i][jj]*delqsietaY + Mal.C[i][jj]*deletaetaY);

        Mal.VTX [i] = alfa*omega*Mal.Lx[i][jj] + alfa*Mal.f[i][jj-1];
        Mal.VTY [i] = alfa*omega*Mal.Ly[i][jj] + alfa*Mal.g[i][jj-1];

            }else{
////         Esquema centrado sugerido parao problema
//        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[i-1][jj]);
//        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[i-1][jj]);
//        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
//        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);
//
        // // Esquema do carderno
        xqsi = 0.5*(Mal.x[i+1][jj+1] - Mal.x[i-1][jj+1]);
        yqsi = 0.5*(Mal.y[i+1][jj+1] - Mal.y[i-1][jj+1]);
        xeta = Mal.x[i][jj+1] - Mal.x[i][jj];
        yeta = Mal.y[i][jj+1] - Mal.y[i][jj];

        Mal.A[i][jj] = xeta*xeta + yeta*yeta;
        Mal.B[i][jj] = xqsi*xeta + yqsi*yeta;
        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;

        Mal.AT [i-1] = - Mal.A[i][jj]; // LOW DIAGONAL
        Mal.BT [i] =   alfa + 2.0*Mal.A[i][jj]; // MAIN DIAGONAL
        Mal.CT [i+1] = - Mal.A[i][jj]; // UPPER DIAGONAL

        delqsiqsiX = Mal.x[i+1][jj] -2.0*Mal.x[i][jj] + Mal.x[i-1][jj];
        delqsietaX = 0.25*(Mal.x[i+1][jj+1] - Mal.x[i+1][jj-1] - Mal.x[i-1][jj+1] + Mal.x[i-1][jj-1]);
        deletaetaX = Mal.x[i][jj+1] -2.0*Mal.x[i][jj] + Mal.x[i][jj-1];

        delqsiqsiY = Mal.y[i+1][jj] -2.0*Mal.y[i][jj] + Mal.y[i-1][jj];
        delqsietaY = 0.25*(Mal.y[i+1][jj+1] - Mal.y[i+1][jj-1] - Mal.y[i-1][jj+1] + Mal.y[i-1][jj-1]);
        deletaetaY = Mal.y[i][jj+1] -2.0*Mal.y[i][jj] + Mal.y[i][jj-1];

//        Mal.Lx[i][jj] = (Mal.A[i][jj]*delqsiqsiX -2.0*Mal.B[i][jj]*delqsietaX + Mal.C[i][jj]*deletaetaX);
//        Mal.Ly[i][jj] = (Mal.A[i][jj]*delqsiqsiY -2.0*Mal.B[i][jj]*delqsietaY + Mal.C[i][jj]*deletaetaY);

        // Multiplicando B por 0.5
        Mal.Lx[i][jj] = (Mal.A[i][jj]*delqsiqsiX -2.0*Mal.B[i][jj]*delqsietaX + Mal.C[i][jj]*deletaetaX);
        Mal.Ly[i][jj] = (Mal.A[i][jj]*delqsiqsiY -2.0*Mal.B[i][jj]*delqsietaY + Mal.C[i][jj]*deletaetaY);

        //cout << "Lx = " << Mal.Lx[i][jj]<< endl;
        //cout << "Ly = " << Mal.Ly[i][jj] << endl;

        Mal.VTX [i] = alfa*omega*Mal.Lx[i][jj] + alfa*Mal.f[i][jj-1];
        Mal.VTY [i] = alfa*omega*Mal.Ly[i][jj] + alfa*Mal.g[i][jj-1];

        //cout << "VTX = " << Mal.VTX [i] << endl;
        //cout << "VTY = " << Mal.VTY [i] << endl;

        }
    }

        //Mal.CT [0] = alfa - Mal.A[M][jj];

        Mal.CT [0] = 0.0;
        //cout << Mal.CT [0] << endl;
        //Vetor::imprimir(Mal.CT,IMAX);

        TridiagPX(Mal.AT, Mal.BT, Mal.CT, Mal.VTX, x, M);
        TridiagPY(Mal.AT, Mal.BT, Mal.CT, Mal.VTY, y, M);

        for (i =0; i<=M; i++){
            Mal.f[i][jj] = x[i];
            Mal.g[i][jj] = y[i];
        }

}

// Resolve o Segundo Sweep do Metodo resolvendo agora uma tridiagonal simples pelo método TDMA
for (int i = 0; i<=IMAX-2; i++){

    for (int jj = 1; jj <=n; jj++){

        if (i == 0){

//        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[IMAX-2][jj]);
//        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[IMAX-2][jj]);

                // Esquema do carderno
        xqsi = 0.5*(Mal.x[i+1][jj+1] - Mal.x[IMAX-2][jj+1]);
        yqsi = 0.5*(Mal.y[i+1][jj+1] - Mal.y[IMAX-2][jj+1]);

        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;


        }else{

//        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[i-1][jj]);
//        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[i-1][jj]);

                // // Esquema do carderno
        xqsi = 0.5*(Mal.x[i+1][jj+1] - Mal.x[i-1][jj+1]);
        yqsi = 0.5*(Mal.y[i+1][jj+1] - Mal.y[i-1][jj+1]);

        }

//        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
//        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);

        // Caderno
        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj]);
        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj]);


        //Mal.A[i][jj] = xeta*xeta + yeta*yeta;
        //Mal.B[i][jj] = xqsi*xeta + yqsi*yeta;
        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;

        //cout << Mal.C[i][jj] << endl;

//        Mal.AT [jj-1] = alfa -Mal.C[i][jj]; // LOW DIAGONAL
//        Mal.BT [jj-1] = alfa + 2.0*Mal.C[i][jj] ; // MAIN DIAGONAL
//        Mal.CT [jj-1] = alfa -Mal.C[i][jj] ; // UPPER DIAGONAL


        Mal.AT [jj-1] = 0.0; // LOW DIAGONAL
        Mal.BT [jj-1] = alfa + Mal.C[i][jj] ; // MAIN DIAGONAL
        Mal.CT [jj-1] = -Mal.C[i][jj] ; // UPPER DIAGONAL

        Mal.VTX [jj-1] =  Mal.f[i][jj] ;
        Mal.VTY [jj-1] =  Mal.g[i][jj] ;


       // }

        if ( jj == 1){
         Mal.AT[0] = 0.0;
        }

        if (jj == n){
        Mal.CT[n] = 0.0;
        }

    }

        TDMAX(Mal.AT, Mal.BT, Mal.CT, Mal.VTX, x, n);
        TDMAY(Mal.AT, Mal.BT, Mal.CT, Mal.VTY, y, n);


          for ( int jj = 1; jj <= n; jj++){

            Mal.DELTAx[i][jj] = x[jj-1];// + Mal.x[i][jj] ; // Xi,j(n+1) = Xi,j(n) + DeltaXi,j(n)
            Mal.DELTAy[i][jj] = y[jj-1];// + Mal.y[i][jj]; // Yi,j(n+1) = Yi,j(n) + DeltaYi,j(n)

            }

        // Correção dos Pontos para todos os pontos j da matrix

         for ( int jj = 1; jj <= JMAX-2; jj++){

            Mal.x[i][jj] = Mal.DELTAx[i][jj] + Mal.x[i][jj];// + Mal.x[i][jj] ; // Xi,j(n+1) = Xi,j(n) + DeltaXi,j(n)
            Mal.y[i][jj] = Mal.DELTAy[i][jj] + Mal.y[i][jj];// + Mal.y[i][jj]; // Yi,j(n+1) = Yi,j(n) + DeltaYi,j(n)

            Mal.x[IMAX-1][jj] = Mal.x[0][jj];
            Mal.y[IMAX-1][jj] = Mal.y[0][jj];

          }

}

return Mal;


}

Malha AF2Atraction (Malha Mal,int i, int j, int IMAX, int JMAX, double omega, double alfa){

    double xqsi, yqsi, xeta, yeta, delqsiqsiX, delqsietaX, deletaetaX, delqsiqsiY, delqsietaY, deletaetaY;

    double *x; // Armazenando espaço e zerando o vetor solução da periodica
    x = new double [IMAX];

    double *y; // Armazenando espaço e zerando o vetor solução da periodica
    y = new double [IMAX];


    // Coeficientes da Função de Atração do Tompson


    int L  = 1;
    int MM = 1;

    double *a;
    a = new double [L];
    double *b;
    b = new double [L];
    double *c;
    c = new double [L];
    double *d;
    d = new double [L];

    a[0] = 100;
    b[0] = 5;
    c[0] = 5;
    d[0] = 1;


    int M = IMAX-2; // Penultimo ponto da periodica
    int n = JMAX-2; // Penultimo ponto da periodica

    int jj;

// Cria e Resolve a tridiagonal Periodica ( 1 Parte do AF2)

for (int jj = 1; jj<=JMAX-2; jj++){

    for (i = 0; i <=M; i++){

             if (i == 0){

//         Esquema centrado
//        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[IMAX-2][jj]);
//        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[IMAX-2][jj]);
//        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
//        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);

        // Esquema do carderno
        xqsi = 0.5*(Mal.x[i+1][jj+1] - Mal.x[IMAX-2][jj+1]);
        yqsi = 0.5*(Mal.y[i+1][jj+1] - Mal.y[IMAX-2][jj+1]);
        xeta = Mal.x[i][jj+1] - Mal.x[i][jj];
        yeta = Mal.y[i][jj+1] - Mal.y[i][jj];

        Mal.A[i][jj] = xeta*xeta + yeta*yeta;
        Mal.B[i][jj] = xqsi*xeta + yqsi*yeta;
        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;
        Mal.D[i][jj] = pow((xqsi*yeta - xeta*yqsi),2);


        Mal.AT [M] = - Mal.A[i][jj]; // LOW DIAGONAL
        Mal.BT [i] = alfa + 2.0*Mal.A[i][jj]; // MAIN DIAGONAL
        Mal.CT [i+1] = - Mal.A[i][jj]; // UPPER DIAGONAL

        delqsiqsiX= Mal.x[i+1][jj] -2.0*Mal.x[i][jj] + Mal.x[IMAX-2][jj];
        delqsietaX = 0.25*(Mal.x[i+1][jj+1] - Mal.x[i+1][jj-1] - Mal.x[IMAX-2][jj+1] + Mal.x[IMAX-2][jj-1]);
        deletaetaX = Mal.x[i][jj+1] -2.0*Mal.x[i][jj] + Mal.x[i][jj-1];

        delqsiqsiY= Mal.y[i+1][jj] -2.0*Mal.y[i][jj] + Mal.y[IMAX-2][jj];
        delqsietaY = 0.25*(Mal.y[i+1][jj+1] - Mal.y[i+1][jj-1] - Mal.y[IMAX-2][jj+1] + Mal.y[IMAX-2][jj-1]);
        deletaetaY = Mal.y[i][jj+1] -2.0*Mal.y[i][jj] + Mal.y[i][jj-1];

//        Mal.Lx[i][jj] = (Mal.A[i][jj]*delqsiqsiX -2.0*Mal.B[i][jj]*delqsietaX + Mal.C[i][jj]*deletaetaX);
//        Mal.Ly[i][jj] = (Mal.A[i][jj]*delqsiqsiY -2.0*Mal.B[i][jj]*delqsietaY + Mal.C[i][jj]*deletaetaY);

        PQfunction(Mal,a,b,c,d,i,j,L,MM);

        // Multiplicando B por 0.5
       // Mal.Lx[i][jj] = (Mal.A[i][jj]*delqsiqsiX -2.0*Mal.B[i][jj]*delqsietaX + Mal.C[i][jj]*deletaetaX + Mal.D[i][jj]*(Pqsieta*xqsi + Qqsieta*xeta));
        //Mal.Ly[i][jj] = (Mal.A[i][jj]*delqsiqsiY -2.0*Mal.B[i][jj]*delqsietaY + Mal.C[i][jj]*deletaetaY + Mal.D[i][jj]*(Pqsieta*yqsi + Qqsieta*yeta));

        Mal.VTX [i] = alfa*omega*Mal.Lx[i][jj] + alfa*Mal.f[i][jj-1];
        Mal.VTY [i] = alfa*omega*Mal.Ly[i][jj] + alfa*Mal.g[i][jj-1];

        //cout << "VTX = " << Mal.VTX [i] << endl;
        //cout << "VTY = " << Mal.VTY [i] << endl;
        //cout << "f = " << Mal.f[i][jj-1] << endl;
        //cout << "g = " << Mal.g[i][jj-1] << endl;
            }else{
////         Esquema centrado sugerido parao problema
//        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[i-1][jj]);
//        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[i-1][jj]);
//        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
//        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);
//
        // // Esquema do carderno
        xqsi = 0.5*(Mal.x[i+1][jj+1] - Mal.x[i-1][jj+1]);
        yqsi = 0.5*(Mal.y[i+1][jj+1] - Mal.y[i-1][jj+1]);
        xeta = Mal.x[i][jj+1] - Mal.x[i][jj];
        yeta = Mal.y[i][jj+1] - Mal.y[i][jj];

        Mal.A[i][jj] = xeta*xeta + yeta*yeta;
        Mal.B[i][jj] = xqsi*xeta + yqsi*yeta;
        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;
        Mal.D[i][jj] = pow((xqsi*yeta - xeta*yqsi),2);

        Mal.AT [i-1] = - Mal.A[i][jj]; // LOW DIAGONAL
        Mal.BT [i] =   alfa + 2.0*Mal.A[i][jj]; // MAIN DIAGONAL
        Mal.CT [i+1] = - Mal.A[i][jj]; // UPPER DIAGONAL

        delqsiqsiX = Mal.x[i+1][jj] -2.0*Mal.x[i][jj] + Mal.x[i-1][jj];
        delqsietaX = 0.25*(Mal.x[i+1][jj+1] - Mal.x[i+1][jj-1] - Mal.x[i-1][jj+1] + Mal.x[i-1][jj-1]);
        deletaetaX = Mal.x[i][jj+1] -2.0*Mal.x[i][jj] + Mal.x[i][jj-1];

        delqsiqsiY = Mal.y[i+1][jj] -2.0*Mal.y[i][jj] + Mal.y[i-1][jj];
        delqsietaY = 0.25*(Mal.y[i+1][jj+1] - Mal.y[i+1][jj-1] - Mal.y[i-1][jj+1] + Mal.y[i-1][jj-1]);
        deletaetaY = Mal.y[i][jj+1] -2.0*Mal.y[i][jj] + Mal.y[i][jj-1];

//        Mal.Lx[i][jj] = (Mal.A[i][jj]*delqsiqsiX -2.0*Mal.B[i][jj]*delqsietaX + Mal.C[i][jj]*deletaetaX);
//        Mal.Ly[i][jj] = (Mal.A[i][jj]*delqsiqsiY -2.0*Mal.B[i][jj]*delqsietaY + Mal.C[i][jj]*deletaetaY);

        // Multiplicando B por 0.5
        //Mal.Lx[i][jj] = (Mal.A[i][jj]*delqsiqsiX -2.0*Mal.B[i][jj]*delqsietaX + Mal.C[i][jj]*deletaetaX + Mal.D[i][jj]*(Pqsieta*xqsi + Qqsieta*xeta));
        //Mal.Ly[i][jj] = (Mal.A[i][jj]*delqsiqsiY -2.0*Mal.B[i][jj]*delqsietaY + Mal.C[i][jj]*deletaetaY + Mal.D[i][jj]*(Pqsieta*yqsi + Qqsieta*yeta));

        //cout << "Lx = " << Mal.Lx[i][jj]<< endl;
        //cout << "Ly = " << Mal.Ly[i][jj] << endl;

        Mal.VTX [i] = alfa*omega*Mal.Lx[i][jj] + alfa*Mal.f[i][jj-1];
        Mal.VTY [i] = alfa*omega*Mal.Ly[i][jj] + alfa*Mal.g[i][jj-1];

        //cout << "VTX = " << Mal.VTX [i] << endl;
        //cout << "VTY = " << Mal.VTY [i] << endl;

        }
    }

        //Mal.CT [0] = alfa - Mal.A[M][jj];

        Mal.CT [0] = - Mal.A[M][jj];
        //cout << Mal.CT [0] << endl;
        //Vetor::imprimir(Mal.CT,IMAX);

        TridiagPX(Mal.AT, Mal.BT, Mal.CT, Mal.VTX, x, M);
        TridiagPY(Mal.AT, Mal.BT, Mal.CT, Mal.VTY, y, M);

        for (i =0; i<=M; i++){

            Mal.f[i][jj] = x[i];
            Mal.g[i][jj] = y[i];

        }

        //Mal.f[IMAX-1][jj] = Mal.f[0][jj]; Garanto que é igual
        //Mal.g[IMAX-1][jj] = Mal.g[0][jj];
}

// Resolve o Segundo Sweep do Metodo resolvendo agora uma tridiagonal simples pelo método TDMA
for (int i = 0; i<=IMAX-2; i++){

    for (int jj = 1; jj <=n; jj++){

        if (i == 0){

//        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[IMAX-2][jj]);
//        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[IMAX-2][jj]);

                // Esquema do carderno
        xqsi = 0.5*(Mal.x[i+1][jj+1] - Mal.x[IMAX-2][jj+1]);
        yqsi = 0.5*(Mal.y[i+1][jj+1] - Mal.y[IMAX-2][jj+1]);

        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;


        }else{

//        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[i-1][jj]);
//        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[i-1][jj]);

                // // Esquema do carderno
        xqsi = 0.5*(Mal.x[i+1][jj+1] - Mal.x[i-1][jj+1]);
        yqsi = 0.5*(Mal.y[i+1][jj+1] - Mal.y[i-1][jj+1]);

        }

//        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
//        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);

        // Caderno
        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj]);
        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj]);


        //Mal.A[i][jj] = xeta*xeta + yeta*yeta;
        //Mal.B[i][jj] = xqsi*xeta + yqsi*yeta;
        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;

        //cout << Mal.C[i][jj] << endl;

//        Mal.AT [jj-1] = alfa -Mal.C[i][jj]; // LOW DIAGONAL
//        Mal.BT [jj-1] = alfa + 2.0*Mal.C[i][jj] ; // MAIN DIAGONAL
//        Mal.CT [jj-1] = alfa -Mal.C[i][jj] ; // UPPER DIAGONAL


        Mal.AT [jj-1] = 0.0; // LOW DIAGONAL
        Mal.BT [jj-1] = alfa + Mal.C[i][jj] ; // MAIN DIAGONAL
        Mal.CT [jj-1] = -Mal.C[i][jj] ; // UPPER DIAGONAL

       // cout << "AT = " << Mal.AT[jj-1] << endl;
        //cout << "BT = " << Mal.BT[jj-1]<< endl;
        //cout << "CT = " << Mal.CT[jj-1] << endl;

        Mal.VTX [jj-1] =  Mal.f[i][jj] ;
        Mal.VTY [jj-1] =  Mal.g[i][jj] ;

        //cout << "VTX = " << Mal.VTX [jj-1]<< endl;
        //cout << "VTY = " << Mal.VTY [jj-1]<< endl;

       // }

        if ( jj == 1){
         Mal.AT[0] = 0.0;
        }

        if (jj == n){
        Mal.CT[n] = 0.0;
        }

        //cout << "AT = " << Mal.AT[jj-1] << endl;
       // cout << "BT = " << Mal.BT[jj-1]<< endl;
       // cout << "CT = " << Mal.CT[jj-1] << endl;

    }

        TDMAX(Mal.AT, Mal.BT, Mal.CT, Mal.VTX, x, n);
        TDMAY(Mal.AT, Mal.BT, Mal.CT, Mal.VTY, y, n);


          for ( int jj = 1; jj <= n; jj++){

            Mal.DELTAx[i][jj] = x[jj-1];// + Mal.x[i][jj] ; // Xi,j(n+1) = Xi,j(n) + DeltaXi,j(n)
            Mal.DELTAy[i][jj] = y[jj-1];// + Mal.y[i][jj]; // Yi,j(n+1) = Yi,j(n) + DeltaYi,j(n)

           // cout << "DELTAx = " << Mal.DELTAx[i][jj] << endl;
           // cout << "DELTAy = " << Mal.DELTAy[i][jj] << endl;

            }

        // Correção dos Pontos para todos os pontos j da matrix

         for ( int jj = 1; jj <= JMAX-2; jj++){

            //cout << "Xold = " << Mal.x[i][jj] << endl;
            //cout << "Yold = " << Mal.y[i][jj] << endl;

            Mal.x[i][jj] = Mal.DELTAx[i][jj] + Mal.x[i][jj];// + Mal.x[i][jj] ; // Xi,j(n+1) = Xi,j(n) + DeltaXi,j(n)
            Mal.y[i][jj] = Mal.DELTAy[i][jj] + Mal.y[i][jj];// + Mal.y[i][jj]; // Yi,j(n+1) = Yi,j(n) + DeltaYi,j(n)

           // cout << "X = " << Mal.x[i][jj] << endl;
           // cout << "Y = " << Mal.y[i][jj] << endl;

            Mal.x[IMAX-1][jj] = Mal.x[0][jj];
            Mal.y[IMAX-1][jj] = Mal.y[0][jj];

          }

}

return Mal;


}

void PQfunction(Malha Mal,double* a, double* b, double* c, double* d, int i, int jj, int L, int MM){

    double *PqsietaL, *QqsietaL, *PqsietaM, *QqsietaM;


    PqsietaL[0] = 0.0;
    QqsietaL[0] = 0.0;
    PqsietaM[0] = 0.0;
    QqsietaM[0] = 0.0;

    for (int l = 0; l <L; l++){

        double sinal;

        if (Mal.x[i][jj] ==0){
        sinal = 0.0;
        }
        if (Mal.x[i][jj] > 0.0){
            sinal =1.0;
        }else {
        sinal = -1.0;
        }

    PqsietaL[l+1] = Pqsieta[l] + a[l]*sinal*(Mal.x[i][jj] - Mal.x[l][jj])*exp(-c[l]*valorAbsoluto(Mal.x[i][jj] - Mal.x[l][jj]));
    Qqsieta[l+1] =2;


    }



 return;

}

double valorAbsoluto(double x){
      if (x >= 0.0) return x;
      else return -x;
    }
//double sgn(Malha Mal,double **M){
//
//double sinal;
//
//if (Mal.x[i][jj] ==0){
//
//    sinal = 0.0;
//}
//
//if (Mal.x[i][jj] > 0.0){
//    sinal =1.0;
//
//}else {
//sinal = -1.0;
//
//}
//
//return sinal;
//}

float FloatMachineEps(){

    float fmachine_e, ftest;
    fmachine_e = 1.0;

    ftest = 1.0 + fmachine_e;

    while (1.0 != ftest){

        fmachine_e = fmachine_e/2.0;
        ftest = 1.0 + fmachine_e;
    }

    return fmachine_e;
}

double DoubleMachineEps(){
    double dmachine_e, dtest;
    dmachine_e = 1.0;

    dtest = 1.0 +dmachine_e;

 while (1.0 != dtest){

        dmachine_e = dmachine_e/2.0;
        dtest = 1.0 + dmachine_e;
    }

    return dmachine_e;
}

long double LongDoubleMachineEps(){
    long double dmachine_e, dtest;
    dmachine_e = 1.0;

    dtest = 1.0 +dmachine_e;

 while (1.0 != dtest){

        dmachine_e = dmachine_e/2.0;
        dtest = 1.0 + dmachine_e;
    }

    return dmachine_e;
}

void TransMatrixVetor(double** Mx, double** My, int IMAX, int JMAX, double* xVet, double* yVet){

    int t1 = 0, t2 = 0;

    for(int j = 0; j < JMAX; j++){

        for(int i = 0; i < IMAX ; i++){

            xVet[t1++] = Mx[i][j];

            yVet[t2++] = My[i][j];

        }
    }


}
