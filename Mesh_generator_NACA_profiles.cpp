/*
Project 2 - Computational Mesh Generation on NACA Profiles
Computer code written by Alisson Vinicius Brito Lopes
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

    double MatAgrVet[qtdlinhas][2]; // Imprimir no formato Tecplot

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

// Criando ponteiros duplos e simples para Alocação dinâmica de memoria e "zeramento" de vetores e matrizes na memoria do computador
    double **Lx,**Ly, **Cc,**f, **g;
    double *RHS, *AT, *BT, *CT, *DT,*VTX,*ATY, *VTY;
    double **A, **B, **C, **D, **P, **Q, **AMX, **AMY;
    double **x,**y;
    double **xint, **yint, **xo,**yo;
    double *eref,*s;
    double *Rx, *R, *Ry, *Seta, *Dels;

void inicializaMatriz_Dels( int JMAX){
   int j;

    Dels = new double [JMAX];

    for (j = 0; j < JMAX; j++){

           Dels[j] = 0.0;
    }
}

void inicializaMatriz_Seta( int JMAX){
   int j;

    Seta = new double [JMAX];

    for (j = 0; j < JMAX; j++){

           Seta[j] = 0.0;
    }
}

void inicializaMatriz_Rx( int JMAX){
   int j;

    Rx = new double [JMAX];

    for (j = 0; j < JMAX; j++){

           Rx[j] = 0.0;
    }
}

void inicializaMatriz_Ry( int JMAX){
   int j;

    Ry = new double [JMAX];

    for (j = 0; j < JMAX; j++){

           Ry[j] = 0.0;
    }
}

void inicializaMatriz_R( int JMAX){
   int j;

    R = new double [JMAX];

    for (j = 0; j < JMAX; j++){

           R[j] = 0.0;
    }
}

void inicializaMatriz_s( int JMAX){
   int j;

    s = new double [JMAX];

    for (j = 0; j < JMAX; j++){

           s[j] = 0.0;
    }
}

void inicializaMatriz_eref( int JMAX){
   int j;

    eref = new double [JMAX];

    for (j = 0; j < JMAX; j++){

           eref[j] = 0.0;
    }
}

void inicializaMatriz_xint(int IMAX, int JMAX){

  int i, j;

  xint = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        xint[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 xint[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_yint(int IMAX, int JMAX){

  int i, j;

  yint = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        yint[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 yint[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_xo(int IMAX, int JMAX){

  int i, j;

  xo = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        xo[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 xo[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_yo(int IMAX, int JMAX){

  int i, j;

  yo = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        yo[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 yo[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_x(int IMAX, int JMAX){

  int i, j;

  x = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        x[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 x[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_y(int IMAX, int JMAX){

  int i, j;

  y = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        y[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 y[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_Lx(int IMAX, int JMAX){

  int i, j;

  Lx = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        Lx[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 Lx[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_Ly(int IMAX, int JMAX){

  int i, j;

  Ly = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        Ly[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 Ly[i][j] = 0.0;
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

    CT = new double [IMAX];
    for (j = 0; j < IMAX; j++){

           CT[j] = 0.0;
    }
}

void inicializaMatriz_DT( int IMAX){
   int j;

    DT = new double [IMAX];

    for (j = 0; j <IMAX; j++){

           DT[j] = 0.0;
    }
}

void inicializaMatriz_VTX( int IMAX){
   int j;

    VTX = new double [IMAX];

    for (j = 0; j <IMAX; j++){

           VTX[j] = 0.0;
    }
}

void inicializaMatriz_ATY( int IMAX){
   int j;

    ATY = new double [IMAX];

    for (j = 0; j < IMAX; j++){

           ATY[j] = 0.0;
    }
}

void inicializaMatriz_VTY( int IMAX){
   int j;

    VTY = new double [IMAX];

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

void inicializaMatriz_AMX(int IMAX, int JMAX){

  int i, j;

  AMX = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        AMX[i] = new double[IMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < IMAX; j++){
                AMX[i][j] = 0.0;
             }
    }
}

void inicializaMatriz_AMY(int IMAX, int JMAX){

  int i, j;

  AMY = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        AMY[i] = new double[IMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < IMAX; j++){
                AMY[i][j] = 0.0;
             }
    }
}

void inicializaMatriz_f(int IMAX, int JMAX){

  int i, j;

  f = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        f[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 f[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_g(int IMAX, int JMAX){

  int i, j;

  g = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        g[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 g[i][j] = 0.0;
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

    // Inicialização para o Algoritmo de THOMAS
    inicializaMatriz_AT(IMAX);
    inicializaMatriz_BT(IMAX);
    inicializaMatriz_CT(IMAX);
    inicializaMatriz_DT(IMAX);

    inicializaMatriz_VTX(IMAX);
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

 }

};

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
Malha aerofolioPontosNACA (Malha Mal, int i, int j, int IMAX, int JMAX,double R,double th, double XSF);
Malha AF2Atraction (Malha Mal,int i, int j, int IMAX, int JMAX, double omega, double alfa);
void TransMatrixVetor(double** Mx,double** My, int IMAX, int JMAX, double* xVet, double* yVet);
double PQfunction(Malha Mal,int i, int j);
double valorAbsoluto(double x);
Malha AF1alfa (Malha Mal,int i, int j, int IMAX, int JMAX, double omega, double alfa,int Melements,double alfaLow,double alfaHigh);
void AlfaSequence (int Melements, double alfaLow, double alfaHigh,double *alfaa);
float FloatMachineEps();
double DoubleMachineEps();
long double LongDoubleMachineEps();

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
// Intel(R) Core(TM) i5-3210M CPU @ 2.50 GHz, sistema operacional Windows 7 de 64 Bits

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

int IMAX = 93; // Número maximo de pontos na direção qsi >> i
int JMAX = 15; // Número maximo de pontos na direção Eta >> j

// Vetores para Transformação da Matrix para imprimir no formato Tecplot/Fortran (permite impressão tecplot)
double* xVet = new double[IMAX*JMAX];
double* yVet = new double[IMAX*JMAX];

double *resmaxiterx, *resmaxitery, RMAXX = 1.0,RMAXY = 1.0;
resmaxiterx = new double [maxit];
resmaxitery = new double [maxit];

int i,j,kiter=1; // i direção qsie j direção eta

double XSF = 1.18; // fator de estiramento ("Stretching") da Mal para a direção x
double YSF = 1.25; // fator de estiramento ("Stretching") da Mal para a direção y

 /*%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    DADOS DE ENTRADA    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

double r = 1.9; // Fator de Relaxação do SLOR
// Raio do Circulo externo
double R = 6.5;
// Airfoil thickness ratio
double th = 0.10;
// Iterations parameters AF2
double alfa = 0.05; // Combinação Boa para meu problema 0.1299 e omega = 1,90 discretização sala
double omega = 1.88;
// Iterations parameters for Numerical Methods and Iteration Parameters for Eingevalues anihhilation
int Melements = 5;
double alfaLow = 0.001;
double alfaHigh = 0.35; // Utilizando o Menor intervalo conforme recomendando por Jameson
// Convergence Criterion
double eps = 1.0*pow(10,-12); //
// Calculo do deltaX na direção X apenas no intervalo em cima do perfil biconvexo
double deltaX;

Malha Mal;
Mal.inicializarMatrizes(IMAX,JMAX,th,alfa,omega,eps,deltaX);

Mal = aerofolioPontosBiconvexo (Mal,i,j,IMAX,JMAX,R,th,XSF);
//Mal = aerofolioPontosNACA (Mal,i,j,IMAX,JMAX,R,th,XSF);

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

// Avalia o Parametro B que esta relacionado com ortogonalidade da malha
char nomeArquivo4[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto2/OrtoBparabolica.txt";
WriteFile::writeMatrix(Mal.B,IMAX,JMAX,nomeArquivo4);

//Imprimi Malha Parabolica
TransMatrixVetor(Mal.x,Mal.y,IMAX,JMAX,xVet,yVet);
char nomeArquivoParabolicaBiconvexo[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto2/Graficos/MalhaParabolicaNACA.plt";
WriteFile::writeVetAgrup(xVet,yVet,(IMAX*JMAX),nomeArquivoParabolicaBiconvexo);

    double *alfaa;
    alfaa = new double [Melements];
    int kk =0;
    AlfaSequence(Melements,alfaLow,alfaHigh,alfaa);

// LOOP PARA CALCULO ITERATIVO DA MALHA ELIPTICA

//while (eps <= RMAXX){
for (kiter = 0; kiter < maxit; kiter ++){

    // Mal = AF1alfa (Mal,i,j,IMAX,JMAX,omega,alfaa[kk],Melements,alfaLow,alfaHigh);
    // Mal = AF1 (Mal,i,j,IMAX,JMAX,omega,alfa);
     Mal = AF2 (Mal,i,j,IMAX,JMAX,omega,alfa);
    // Mal = AF2 (Mal,i,j,IMAX,JMAX,omega,alfaa[kk]); // With Alfa sequence test
    // Mal = AF2Atraction (Mal,i,j,IMAX,JMAX,omega,alfa);
    // Mal = SLOR (Mal,i, j,IMAX,JMAX,deltaX,r,omega);

    RMAXX = MaxResiduoX (Mal,IMAX,JMAX,kiter,resmaxX);
    RMAXY = MaxResiduoY (Mal,IMAX,JMAX,kiter,resmaxY);

    resmaxiterx[kiter] =log10(RMAXX);
    resmaxitery[kiter] =log10(RMAXY);

     cout << "ResX = " << resmaxiterx[kiter] << endl; // Imprimir residuo na tela para monitoramente

    if ( kk>Melements-2){
        kk =0;
    }else{
       kk = kk+1;
    }
}

// kiter = kiter +1;
// cout << kiter <<endl;
// cout << "Numero Total de Iterações = " <<kiter<<endl;

double tempo_total_s = (clock()- TempoInicial)/(double)CLOCKS_PER_SEC;
cout << "Tempo Total de Execucao (segundos) = " << tempo_total_s << endl;
cout << "Numero Total de Iterações = " <<kiter<<endl;
cout << " Tempo por Iteração" << tempo_total_s /kiter << endl;

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

 //Implementação com Estiramento usando Progressão Geometrica (esta sera apresentada no relatorio)
        int numPointsPG = ((IMAX-1)/4);
        double somaPG = 0.5;
        double dx = somaPG*(XSF-1)/(pow(XSF,numPointsPG)-1.0);

        cout << dx <<endl;

        Mal.x[0][0] = 0.0;
        Mal.x[1][0] = Mal.x[0][0] + dx;

        Mal.y[0][0] = 2*th*Mal.x[0][0]*(1.0 -Mal.x[0][0]);
        Mal.y[1][0] = -2*th*Mal.x[1][0]*(1.0 -Mal.x[1][0]);

        Mal.x[(IMAX-1)/2][0] = 1.0;
        Mal.x[((IMAX-1)/2) -1][0] = Mal.x[(IMAX-1)/2][0] - dx;

        Mal.y[(IMAX-1)/2][0] = -2*th*Mal.x[(IMAX-1)/2][0]*(1.0 - Mal.x[(IMAX-1)/2][0]);
        Mal.y[(IMAX-1)/2 -1][0] = -2*th*Mal.x[((IMAX-1)/2) -1][0]*(1.0 - Mal.x[((IMAX-1)/2) -1][0]);

    for ( i = 2; i <=((IMAX-1)/4); i++){ // OK
    // Gerando a PG para estiramento x aerofolio
        Mal.x[i][0] = (Mal.x[i-1][0] + XSF*(Mal.x[i-1][0] - Mal.x[i-2][0]));
        Mal.y[i][0] = -2*th*Mal.x[i][0]*(1.0 - Mal.x[i][0]);

    }

    // Gerando o segundo quadrante
    for ( i =((IMAX-1)/2 -2); i >=((IMAX-1)/4 +1) ; i-- ){
        Mal.x[i][0] = Mal.x[i+1][0] - (Mal.x[i+2][0] - Mal.x[i+1][0])*XSF;
        Mal.y[i][0] = -2*th*Mal.x[i][0]*(1.0 - Mal.x[i][0]);
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
        k = k-1;
    }

    Mal.x[IMAX-1][0] = Mal.x[0][0];
    Mal.y[IMAX-1][0] = Mal.y[0][0];

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

        if (i==0){

        xqsi= ( Mal.x[i+1][j-1]- Mal.x[IMAX-2][j-1] )/2.0;
        yqsi= ( Mal.y[i+1][j-1]- Mal.y[IMAX-2][j-1] )/2.0;

        }else{

        xqsi= ( Mal.x[i+1][j-1]- Mal.x[i-1][j-1] )/2.0;
        yqsi= ( Mal.y[i+1][j-1]- Mal.y[i-1][j-1] )/2.0;

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

        TridiagPX(Mal.AT, Mal.BT, Mal.CT, Mal.VTX, xx, M);
        TridiagPY(Mal.AT, Mal.BT, Mal.CT, Mal.VTY, yy, M);

        for (i =0; i<=M; i++){

            Mal.DELTAx[i][jj] = xx[i];
            Mal.DELTAy[i][jj] = yy[i];

        }

           for ( int i = 0; i <= M; i++){

            Mal.x[i][jj] = Mal.DELTAx[i][jj] + Mal.x[i][jj];
            Mal.y[i][jj] = Mal.DELTAy[i][jj] + Mal.y[i][jj];
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
    int n = JMAX-2; //
    int jj;

// Cria e Resolve a tridiagonal Periodica ( 1 Parte do AF1)

for (int jj = 1; jj<=JMAX-2; jj++){

    for (i = 0; i <=M; i++){

             if (i == 0){

         //Esquema centrado
        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[IMAX-2][jj]);
        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[IMAX-2][jj]);
        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);

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
         //Esquema centrado sugerido parao problema
        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[i-1][jj]);
        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[i-1][jj]);
        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);

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
}

    // Resolve o Segundo Sweep do Metodo resolvendo agora uma tridiagonal simples pelo método TDMA
for (int i = 0; i<=IMAX-2; i++){

    for (int jj = 1; jj <=n; jj++){


        if (i == 0){

        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[IMAX-2][jj]);
        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[IMAX-2][jj]);

        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;

        }else{

        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[i-1][jj]);
        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[i-1][jj]);

        }

        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);

        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;

        Mal.AT [jj-1] = -Mal.C[i][jj]; // LOW DIAGONAL
        Mal.BT [jj-1] = alfa + 2.0*Mal.C[i][jj] ; // MAIN DIAGONAL
        Mal.CT [jj-1] = -Mal.C[i][jj] ; // UPPER DIAGONAL

        Mal.VTX [jj-1] =  Mal.f[i][jj] ;
        Mal.VTY [jj-1] =  Mal.g[i][jj] ;

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

            Mal.x[i][jj] = Mal.DELTAx[i][jj] + Mal.x[i][jj];//
            Mal.y[i][jj] = Mal.DELTAy[i][jj] + Mal.y[i][jj];//

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
        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[IMAX-2][jj]);
        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[IMAX-2][jj]);
        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);

        // Esquema do carderno
//        xqsi = 0.5*(Mal.x[i+1][jj+1] - Mal.x[IMAX-2][jj+1]);
//        yqsi = 0.5*(Mal.y[i+1][jj+1] - Mal.y[IMAX-2][jj+1]);
//        xeta = Mal.x[i][jj+1] - Mal.x[i][jj];
//        yeta = Mal.y[i][jj+1] - Mal.y[i][jj];

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

            }else{
////         Esquema centrado sugerido parao problema
        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[i-1][jj]);
        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[i-1][jj]);
        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);

//        xqsi = 0.5*(Mal.x[i+1][jj+1] - Mal.x[i-1][jj+1]);
//        yqsi = 0.5*(Mal.y[i+1][jj+1] - Mal.y[i-1][jj+1]);
//        xeta = Mal.x[i][jj+1] - Mal.x[i][jj];
//        yeta = Mal.y[i][jj+1] - Mal.y[i][jj];

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

        Mal.VTX [i] = alfa*omega*Mal.Lx[i][jj] + alfa*Mal.f[i][jj-1];
        Mal.VTY [i] = alfa*omega*Mal.Ly[i][jj] + alfa*Mal.g[i][jj-1];

        }
    }

        Mal.CT [0] = - Mal.A[M][jj];

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

        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[IMAX-2][jj]);
        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[IMAX-2][jj]);

//        xqsi = 0.5*(Mal.x[i+1][jj+1] - Mal.x[IMAX-2][jj+1]);
//        yqsi = 0.5*(Mal.y[i+1][jj+1] - Mal.y[IMAX-2][jj+1]);

        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;

        }else{

        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[i-1][jj]);
        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[i-1][jj]);

                // // Esquema do carderno
//        xqsi = 0.5*(Mal.x[i+1][jj+1] - Mal.x[i-1][jj+1]);
//        yqsi = 0.5*(Mal.y[i+1][jj+1] - Mal.y[i-1][jj+1]);

        }

        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);

        // Caderno
//        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj]);
//        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj]);

        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;

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

void TridiagPX(double* a, double* b, double* c, double* dv, double *x, int M){

   // Decomposição L-U
   // Implementar a do Hirsh também e comparar o custo computacional dos dois procedimentos

    int i, IMAX=93;

    double *l1,*m2,*n3,*u1,*u2,*u3,*f,*g;
    double Ta, Tb;

    l1 = new double [M];
    m2 = new double [M];
    n3 = new double [M];
    u1 = new double [M];
    u2 = new double [M];
    u3 = new double [M];
    f = new double [M];
    g = new double [M];

    for (i =0; i<=M;i++){ // Uma vez que foi alocado dinamicamente (via pointer) é conveniente zerar as posições alocadas, evitando pegar lixo de máquina

        l1[i] = 0.0;
        m2[i] = 0.0;
        n3[i] = 0.0;
        u1[i] = 0.0;
        u2[i] = 0.0;
        u3[i] = 0.0;
        f[i]  = 0.0;
        g[i]  = dv[i];
    }

    l1[0] = b1[0]; // L1 = B1 ZERO no c++!!!
    m2[0] = a[0]; // M2 = A1
    n3[0] = c[0]; // N3[0] = C[0]
    f[0]  = g[0]/l1[0]; // f ou x_barra1 = L1^-1*y1
    u2[1] = a[M]/l1[0];  // V1 = L1^-1*A0
    Ta = n3[0]*u2[1]; // Ta = N1V1
    Tb = n3[0]*f[0];  //

/////////////////////
// Forward Sweep ////
/////////////////////

    for (i = 1; i <=(M-2); i++){

        u1[i] = c[i]/l1[i-1];
        l1[i] = b1[i] - u1[i]*m2[i-1];
        m2[i] = a[i];
        n3[i] = -n3[i-1]*u1[i];
        f[i]  = (g[i]-m2[i-1]*f[i-1])/l1[i];
        u2[i+1] = -m2[i-1]*u2[i]/l1[i];
        Ta = Ta + n3[i]*u2[i+1];
        Tb = Tb + n3[i]*f[i];
    }
        u1[M-1] = c[M-1]/l1[M-2];
        l1[M-1] = b1[M-1] - m2[M-2]*u1[M-1];
        m2[M-1] = a[M-1] - n3[M-2]*u1[M-1];
         f[M-1] = (g[M-1]-m2[M-2]*f[M-2])/l1[M-1];

         u1[M] = (c[M]-m2[M-2]*u2[M-1])/l1[M-1];
         Ta = Ta + m2[M-1]*u1[M];
         l1[M] = b1[M]-Ta;
         f[M] = (g[M] -m2[M-1]*f[M-1]-Tb)/l1[M];
        /////////////////////////
        // Backward Sweep //////
        ////////////////////////
        x[M] = f[M];
        x[M-1] = f[M-1] - u1[M]*x[M];

        for (i = (M-2); i>=0; i--){
            x[i] = f[i] - u1[i+1]*x[i+1] - u2[i+1]*x[M];
        }

    delete [] l1;
    delete [] m2;
    delete [] n3;
    delete [] u1;
    delete [] u2;
    delete [] u3;
    delete [] f;
    delete [] g;

 return;

}

void TridiagPY(double* a, double* b, double* c, double* dv, double *y, int M){

int i, IMAX=93;

    double *l1,*m2,*n3,*u1,*u2,*u3,*f,*g;
    double Ta, Tb;

    l1 = new double [M];
    m2 = new double [M];
    n3 = new double [M];

    u1 = new double [M];
    u2 = new double [M];
    u3 = new double [M];

    f = new double [M];
    g = new double [M];

    for (i =0; i<=M;i++){ // Uma vez que foi alocado dinamicamente (via pointer) é conveniente zerar as posições alocadas, evitando pegar lixo de máquina

        l1[i] = 0.0;
        m2[i] = 0.0;
        n3[i] = 0.0;
        u1[i] = 0.0;
        u2[i] = 0.0;
        u3[i] = 0.0;
        f[i]  = 0.0;
        g[i]  = dv[i];
    }

    l1[0] = b1[0]; // L1 = B1 // imprimir para achar o erro 07/05 tarde
    //cout << b[0] << endl;

    m2[0] = a[0]; // M1 = A1
   // cout << a[0] << endl;

    n3[0] = c[0]; // N1 = Cn+1
    //cout << c[0] << endl;

    f[0]  = g[0]/l1[0]; //f = L1^-1*y1
    //cout << f[0] << endl;

    u2[1] = a[M]/l1[0];  // V1 = L1^-1*A0
    //cout << u2[1] << endl;

    Ta = n3[0]*u2[1]; // Ta = N1V1
    //cout << Ta << endl;

    Tb = n3[0]*f[0];  //
    //cout << Tb << endl;

         /////////////////////////
        // Forward Sweep //////
        ////////////////////////
    for (i = 1; i <=(M-2); i++){

        u1[i] = c[i]/l1[i-1];
        l1[i] = b1[i] - u1[i]*m2[i-1];
        m2[i] = a[i];
        n3[i] = -n3[i-1]*u1[i];
        f[i]  = (g[i]-m2[i-1]*f[i-1])/l1[i];
        u2[i+1] = -m2[i-1]*u2[i]/l1[i];
        Ta = Ta + n3[i]*u2[i+1];
        Tb = Tb + n3[i]*f[i];
    }

        u1[M-1] = c[M-1]/l1[M-2];
        l1[M-1] = b1[M-1] - m2[M-2]*u1[M-1];
        m2[M-1] = a[M-1] - n3[M-2]*u1[M-1];
         f[M-1] = (g[M-1]-m2[M-2]*f[M-2])/l1[M-1];

         u1[M] = (c[M]-m2[M-2]*u2[M-1])/l1[M-1];
         Ta = Ta + m2[M-1]*u1[M];
         l1[M] = b1[M]-Ta;
         f[M] = (g[M] -m2[M-1]*f[M-1]-Tb)/l1[M];
        /////////////////////////
        // Backward Sweep //////
        ////////////////////////
        y[M] = f[M];
        y[M-1] = f[M-1] - u1[M]*y[M];

        for (i = (M-2); i>=0; i--){
           y[i] = f[i] - u1[i+1]*y[i+1] - u2[i+1]*y[M];
        }

    // Important Key Concepts : Dynamic memory allocation within functions is quite natural, and most
    // programmers have no problem allocating arrays within functions. However, many programmers
    // become negligent and do not deallocate the temporary memory that they needed.

    //BEWARE of memory leaks! For every allocate (new) there should be
    // a deallocate (delete[]) for local variables!

    delete [] l1;
    delete [] m2;
    delete [] n3;
    delete [] u1;
    delete [] u2;
    delete [] u3;
    delete [] f;
    delete [] g;

    // Common Trick : Notice that every time we call this routine we must "allocate" and "deallocate" memory
    // Suppose we are calling this routine over and over using the same size allocation each time.

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
return resmaxY;
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
        double deltax = 1.0/dx;
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

    for (i = 0; i<=((IMAX-1)/2);i++){

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
        kk = kk+1;
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

Malha AF1alfa (Malha Mal,int i, int j, int IMAX, int JMAX, double omega, double alfa,int Melements, double alfaLow, double alfaHigh){

    double xqsi, yqsi, xeta, yeta, delqsiqsiX, delqsietaX, deletaetaX, delqsiqsiY, delqsietaY, deletaetaY;

    double *x; // Armazenando espaço e zerando o vetor solução da periodica
    x = new double [IMAX];

    double *y; // Armazenando espaço e zerando o vetor solução da periodica
    y = new double [IMAX];

    int M = IMAX-2; // Penultimo ponto da periodica
    int n = JMAX-2; // Penultimo ponto da periodica

    int jj;
    int k;

// Cria e Resolve a tridiagonal Periodica ( 1 Parte do AF1)

for (int jj = 1; jj<=JMAX-2; jj++){

    for (i = 0; i <=M; i++){

             if (i == 0){

//         Esquema centrado
        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[IMAX-2][jj]);
        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[IMAX-2][jj]);
        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);

        Mal.A[i][jj] = xeta*xeta + yeta*yeta;
        Mal.B[i][jj] = xqsi*xeta + yqsi*yeta;
        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;

        Mal.AT [M] = - Mal.A[i][jj]; // LOW DIAGONAL
        Mal.BT [i] = alfa+ 2.0*Mal.A[i][jj]; // MAIN DIAGONAL
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

//        Esquema centrado sugerido para o problema
        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[i-1][jj]);
        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[i-1][jj]);
        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);
//
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
}

    // Resolve o Segundo Sweep do Metodo resolvendo agora uma tridiagonal simples pelo método TDMA
for (int i = 0; i<=IMAX-2; i++){

    for (int jj = 1; jj <=n; jj++){

        if (i == 0){

        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[IMAX-2][jj]);
        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[IMAX-2][jj]);

        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;

        }else{

        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[i-1][jj]);
        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[i-1][jj]);

        }

        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);

        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;

        Mal.AT [jj-1] = -Mal.C[i][jj]; // LOW DIAGONAL
        Mal.BT [jj-1] = alfa + 2.0*Mal.C[i][jj] ; // MAIN DIAGONAL
        Mal.CT [jj-1] = -Mal.C[i][jj] ; // UPPER DIAGONAL

        Mal.VTX [jj-1] =  Mal.f[i][jj] ;
        Mal.VTY [jj-1] =  Mal.g[i][jj] ;

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

    double Pqsieta, Qqsieta ;

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

    a[0] = 10;
    b[0] = 5;
    c[0] = 5;
    d[0] = 1;

    int M = IMAX-2; // Penultimo ponto da periodica
    int n = JMAX-2; // Penultimo ponto da periodica

    int jj;
    int l = 0;

// Cria e Resolve a tridiagonal Periodica ( 1 Parte do AF2)

for (int jj = 1; jj<=JMAX-2; jj++){

    for (i = 0; i <=M; i++){

   if (jj <=JMAX-6){

   Pqsieta = 0.2; //PQfunction(Mal,i,jj); testar valores pequenos de P inicialmente
   Qqsieta = 0.02;


   }else{

   Pqsieta = 0.0; //PQfunction(Mal,i,jj); testar valores pequenos de P inicialmente
   Qqsieta = 0.0;

   }
         if (i == 0){

//         Esquema centrado
        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[IMAX-2][jj]);
        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[IMAX-2][jj]);
        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);

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

        PQfunction(Mal,a,b,c,d,i,jj,L,MM); // Chama function P e Q

        // Multiplicando B por 0.5
        Mal.Lx[i][jj] = (Mal.A[i][jj]*delqsiqsiX -2.0*Mal.B[i][jj]*delqsietaX + Mal.C[i][jj]*deletaetaX + Mal.D[i][jj]*(Pqsieta*xqsi + Qqsieta*xeta));
        Mal.Ly[i][jj] = (Mal.A[i][jj]*delqsiqsiY -2.0*Mal.B[i][jj]*delqsietaY + Mal.C[i][jj]*deletaetaY + Mal.D[i][jj]*(Pqsieta*yqsi + Qqsieta*yeta));

        Mal.VTX [i] = alfa*omega*Mal.Lx[i][jj] + alfa*Mal.f[i][jj-1];
        Mal.VTY [i] = alfa*omega*Mal.Ly[i][jj] + alfa*Mal.g[i][jj-1];

            }else{
////         Esquema centrado sugerido parao problema
        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[i-1][jj]);
        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[i-1][jj]);
        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);

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

        // Multiplicando B por 0.5
        Mal.Lx[i][jj] = (Mal.A[i][jj]*delqsiqsiX -2.0*Mal.B[i][jj]*delqsietaX + Mal.C[i][jj]*deletaetaX + Mal.D[i][jj]*(Pqsieta*xqsi + Qqsieta*xeta));
        Mal.Ly[i][jj] = (Mal.A[i][jj]*delqsiqsiY -2.0*Mal.B[i][jj]*delqsietaY + Mal.C[i][jj]*deletaetaY + Mal.D[i][jj]*(Pqsieta*yqsi + Qqsieta*yeta));

        Mal.VTX [i] = alfa*omega*Mal.Lx[i][jj] + alfa*Mal.f[i][jj-1];
        Mal.VTY [i] = alfa*omega*Mal.Ly[i][jj] + alfa*Mal.g[i][jj-1];

        }
    }

        Mal.CT [0] = - Mal.A[M][jj];

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

        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[IMAX-2][jj]);
        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[IMAX-2][jj]);

        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;

        }else{

        xqsi = 0.5*(Mal.x[i+1][jj] - Mal.x[i-1][jj]);
        yqsi = 0.5*(Mal.y[i+1][jj] - Mal.y[i-1][jj]);

        }

        xeta = 0.5*(Mal.x[i][jj+1] - Mal.x[i][jj-1]);
        yeta = 0.5*(Mal.y[i][jj+1] - Mal.y[i][jj-1]);

        Mal.C[i][jj] = xqsi*xqsi + yqsi*yqsi;

        Mal.AT [jj-1] = 0.0; // LOW DIAGONAL
        Mal.BT [jj-1] = alfa + Mal.C[i][jj] ; // MAIN DIAGONAL
        Mal.CT [jj-1] = -Mal.C[i][jj] ; // UPPER DIAGONAL

        Mal.VTX [jj-1] =  Mal.f[i][jj] ;
        Mal.VTY [jj-1] =  Mal.g[i][jj] ;

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

            Mal.DELTAx[i][jj] = x[jj-1];
            Mal.DELTAy[i][jj] = y[jj-1];

            }

        // Correção dos Pontos para todos os pontos j da matrix
         for ( int jj = 1; jj <= JMAX-2; jj++){

            Mal.x[i][jj] = Mal.DELTAx[i][jj] + Mal.x[i][jj];
            Mal.y[i][jj] = Mal.DELTAy[i][jj] + Mal.y[i][jj];

            Mal.x[IMAX-1][jj] = Mal.x[0][jj];
            Mal.y[IMAX-1][jj] = Mal.y[0][jj];

          }
}
return Mal;

}

double PQfunction(Malha Mal,int i, int jj){

// Será apenas um teste
    double Pqsieta, Qqsieta;
    double bb;
    double dd;
    bb = 10.0;
    dd = 1.0;

    //Pqsieta = 0.0;
    int refpointi = 2;
    int refpointjj = 2;

    double sinal1;
    double sinal2;

        if ( (i - refpointi) == 0){
        sinal1 = 0.0;
        }

        if ((i-refpointi) > 0.0){
            sinal1 =1.0;

        }else {

        sinal1 = -1.0;
        }

    Pqsieta = - bb*sinal1*exp(-dd*pow((pow((i - refpointi),2)+pow((jj - refpointjj),2)),0.5));

// Calculo do Q
        if ((jj - refpointjj) == 0){
        sinal2 = 0.0;
        }

        if ((jj-refpointjj) > 0.0){
            sinal2 =1.0;

        }else {

        sinal2 = -1.0;
        }

    Qqsieta = - bb*sinal2*exp(-dd*pow((pow((i - refpointjj),2)+pow((jj - refpointjj),2)),0.5));

    //cout << sinal << endl;
    //cout << Pqsieta << endl;

    //}

  // Laço para M

//    for (int m = 0; m <MM; m++){
//
//        double sinal;
//
//        if (() == 0){
//        sinal = 0.0;
//        }
//        if (Mal.x[i][jj] > 0.0){
//            sinal =1.0;
//        }else {
//        sinal = -1.0;
//        }
//
//    PqsietaM[m+1] = -(PqsietaM[m] + b[m]*sinal*(Mal.x[i][jj] - Mal.x[m][jj])*exp(-d[m]*(pow((Mal.x[i][jj] - Mal.x[m][jj]),2))+ )  );
//    QqsietaM[m+1] = -(QqsietaM[m] + b[m]*sinal*(Mal.x[i][jj] - Mal.x[m][jj])*exp(-c[m]*valorAbsoluto(Mal.x[i][jj] - Mal.x[m][jj]))  );
//
//    cout << PqsietaL[m+1] << endl;
//    cout << QqsietaL[m+1] << endl;
//
//    }

 return Pqsieta, Qqsieta;

}

double valorAbsoluto(double x){
      if (x >= 0.0) return x;
      else return -x;
    }

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

void AlfaSequence (int Melements, double alfaLow, double alfaHigh,double *alfaa){

 double expoente;
 int k, kaux;

    for (k = 0; k<Melements;k++){

        kaux = k+1;

        if (k == Melements){
            kaux = k;
        }

        expoente = (kaux -1.0)/(Melements-1.0);
        alfaa[k] = alfaHigh*pow((alfaLow/alfaHigh),expoente);

    }

return;
}
