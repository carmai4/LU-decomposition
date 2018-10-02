/*
 set the diagonals of U to 1.0
 copy the first column of A into the first column of L
 divide the first row of A by A[0][0] and copy to the first row of U
 you will have to complete the following n-1 times (since first column of L and first row of U are done)
     1) sequentially compute all the rows of a column of L by doing the previous multiplys on known values of x
     2) sequentially compute all the columns of a row of U by doing the previous multiplys on known values of x

 Once L and U are computed
 find the vector y in L * y = b by doing a forward solve
 then find the vector x in U * x = y by doing a back solve
 check the two norm of the error vector by ||b - Ax||
*/

//reference for matrix operation: http://paulbourke.net/miscellaneous/determinant/


#include <iostream>
#include <random>
#include <ctime>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


using namespace std;


    //for printing (n x n) matrix
    void printMat(double**mat, int n)
    {
        if(n == 1)
        {
            cout << mat[0][0] << "\n" << endl;
            return;
        }

        for(int i=0; i<n; ++i)
        {
            cout << " " << endl;
            for(int j=0;j<n;++j)
            {
                cout << " " << mat[i][j] << "  " ;
            }
        }
        cout << "\n" <<endl;
    }

    //for printing (n x 1) matrix, or Vector
    void printVec(double* vec, int n)
    {
        for(int i=0; i<n; ++i)
        {
            cout << " " << vec[i] << endl;
        }
        cout << "\n" << endl;
    }



    //allocates memory
    double** createMat( int n )
    {
        double **a;
        a = (double**)malloc((n)*sizeof(double *));

        for(int i=0; i < n; ++i)
        {
            a[i] = (double*)malloc((n)*sizeof(double));
        }
        return a;
    }



    //get determinant of matrix
    double getDet( double **a,int n )
    {
    //cout << "entering getDet \n" <<endl;
       int i,j,j1,j2;
       double det = 0.0;
       double **m = NULL;

        if (n < 1)
        { /* Error */
        }
        else if (n == 1)
        {
          det = a[0][0];
        }
        else if (n == 2)
        {
          det = (a[0][0] * a[1][1]) - (a[1][0] * a[0][1]);
        }
        else
        {
          det = 0.0;
          for (j1=0; j1<n; j1++)
          {
             m = (double**)malloc((n-1)*sizeof(double *));

             for (i=0; i<n-1; i++)
                m[i] = (double*)malloc((n-1)*sizeof(double));

             for (i=1; i<n; i++)
             {
                j2 = 0;

                for (j=0; j<n; j++)
                {
                   if (j == j1)
                      continue;
                   m[i-1][j2] = a[i][j];
                   j2++;
                }
             }
             det += pow(-1.0,j1+2.0) * a[0][j1] * getDet(m,n-1);
             for (i=0;i<n-1;i++)
                free(m[i]);
             free(m);
          }
        }
        return(det);
    }



    void coFactor(double **a, int n, double **b)
    {
cout <<"entering coFactor " <<endl;
        int i,j,ii,jj,i1,j1;
        double det = 0.0;
        double **c;

        c = (double**)malloc((n-1)*sizeof(double *));
        for (i=0; i<n-1; i++)
        {
            c[i] = (double*)malloc((n-1)*sizeof(double));
        }

        for (j=0; j < n; j++)
        {
            for (i=0; i < n; i++)
            {
                i1 = 0;
                for (ii=0; ii < n; ii++)
                {
                    if (ii == i)
                       continue;

                    j1 = 0;

                    for (jj=0; jj < n; jj++)
                    {
                       if (jj == j)
                          continue;

                       c[i1][j1] = a[ii][jj];
                       j1++;
                    }
                    i1++;
                }

cout << " \nc = "<<endl;
printMat(c , n-1);

             det = getDet( c , n-1 );

    cout << " in for loop, det = " << det << " , j=" << j << " , i=" << i << endl;

             b[i][j] = pow(-1.0 , i+j) * det;
            }
        }

        for (i=0; i<n-1; i++)
        {
            free(c[i]);
        }   free(c);
    }

    //transpose matrix
    void transpose(double **a, int n)
    {
        cout<< "\n before transpose: \n" << endl;
        printMat(a,n);

       int i,j;
       double tmp;

        for (i=1; i<n; i++)
        {
            for (j=0; j<i; j++)
            {
                 tmp = a[i][j];
                 a[i][j] = a[j][i];
                 a[j][i] = tmp;
            }
       }
    }


    void getInv(double** a, double **ainv, int n)
    {
        double det = getDet(a , n);
cout <<" det = " << det << "\n"<<endl;

        double** cofac;
        cofac = createMat(n);
//cout << " created cofac[] " << endl;

        coFactor(a, n, cofac);//get cofactor matrix
//cout << " called coFactor()" << endl;

        transpose(cofac, n); //next, transpose cofactor matrix


        for(int i=0; i < n; ++i) //divide each element by det
        {
            for(int j=0; j < n; ++j)
            {
                ainv[i][j]  = cofac[i][j] / det;
            }
        }
    }


int main()
{
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();

    cout << "\nSpring 2016 CS417 Morris\n Project 2\n\n****************************************************************************\n" << endl;

    int n;
    cout << "size of matrix? ";
    cin >> n;

    myclock::duration d = myclock::now() - beginning;
    unsigned seed3 =  d.count();
    cout << "( seed = " << seed3 << " )\n" << endl;
    minstd_rand0 generator ( seed3 );
    uniform_int_distribution<int> distribution( -654321 , 9876543 );

    double **A;
    A = createMat(n);
    double **L;
    L = createMat(n);
    double **U;
    U = createMat(n);


    cout << "\nA = " << endl;
    for(int r=0; r<n; r++)
    {
        for(int c=0; c<n; c++)
        {
            A[r][c] = double( distribution(generator) ) / 56789;

            if( c == 0 )//copy first col of A into first col of L
                L[r][0] = A[r][0];
        }

        if( r == 0 ) //divide first ro of A by A[0][0] and copy to the first ro of U
        {
            for(int k = 0; k < n; ++k)
            {
                U[0][k] = ( A[0][k] / A[0][0] );
            }
        }
    }
    printMat(A,n);


    //L to have upper triangular entries zero
    for(int r = 0; r < n; ++r)
    {
        for(int c = 0; c < n; ++c)
        {
            if( c > r )
                L[r][c] = 0.0;
        }
    }


    //U to have lower triangular entries zero
    for(int r = 0; r < n; ++r)
    {
        for(int c = 0; c < n; ++c)
        {
            if(c == r) //set the diagonal entries to 1.0
                U[r][c] = 1.0;

            else if(c < r)
                U[r][c] = 0.0;
        }
    }

    //need (n-1)^2 steps to compute all unknown entries of L and U
    double sumL = 0.0;
    double sumU = 0.0;

    for( int k = 1; k < n; ++k )
    {
        //compute entries of L and U
        for( int ro = k; ro < n; ++ro ) //starting at ro = 1
        {
            sumL = 0.0;
            for( int co = 0; co < ro; ++co )
            {
                sumL += ( L[ro][co] * U[co][k] );
            }
            L[ro][k] = ( A[ro][k] - sumL );

            if(ro > k)
            {
                for(int col = k; col < ro; ++col)
                {
                    sumU = 0.0;
                    for(int row = 0; row < col; ++row)
                    {
                        sumU += ( L[col][row] * U[row][ro] );
                    }
                    U[k][ro] = ( A[col][ro] - sumU ) / L[k][k];
                }
            }

        }
    }

//PRINT
    cout << "\nL = " << endl;
    for(int r=0; r<n; r++)
    {
        for(int c=0; c<n; c++)
        {
            cout << "L[" << r << "][" << c << "]=" <<L[r][c] << "     ";
        }
        cout << "\n"<< endl;
    }

    cout << "\nU = " << endl;
    for(int r=0; r<n; r++)
    {
        for(int c=0; c<n; c++)
        {
            cout << "U[" << r << "][" << c << "]="<< U[r][c] << "     ";
        }
        cout << "\n"<< endl;
    }


//is L*U = A ???
/* not working
    double **LxU;
    LxU = createMat(n);

    double sum = 0.0;
    for(int k=0; k < n-1; ++k)
    {
        for(int r=0; r < n-1; ++r)
        {
            for(int col = 0; col < n-1; ++col)
            {
                sum = 0.0;
                for(int c=0; c < n-1; ++c)
                {
                    sum += L[r][c]*U[c][col];
                }
                LxU[r][col] = sum;
            }
        }
    }

    cout << "\n L * U = " << endl;
    printMat(LxU,n);
    */

    // generate random solution b
    double *b;
    b = new double [n];

    cout << "\nGenerated random solution b : " << endl;

    for(int i=0; i < n; i++)
    {
        b[i] = double ( distribution(generator) ) / 56789;

    }
    printVec(b,n);



    //compute (L^-1)
    double **Linv; //inverse of L (L^-1)
    Linv = createMat(n);
    getInv(L, Linv, n);

//PRINT
cout << "\nL^-1 = " << endl;
    for(int r=0; r<n; r++)
    {
        for(int c=0; c<n; c++)
        {
            cout << " Linv[" << r << "][" << c << "]=" << Linv[r][c] << "     ";
        }
        cout << "\n"<< endl;
    }


    // ( Linv * b ) = y
    // forward solve for y in L * y = b
    // y = (L^-1) * b
    double *y;
    y = new double[n];

    for(int r=0; r<n; ++r)
    {
        double sum = 0.0;
        for(int c=0; c<n; ++c)
        {
            sum += ( Linv[r][c] * b[r] );
        }
        y[r] = sum;
    }
    cout << "\ny = " << endl;
    printVec(y,n);



    //backsolve to find x in U * x = y
    double *x;
    x = new double[n];

    for(int i = n-1; i >= 0; --i) //going from bottom to top row
    {
        double sum = 0.0;
        for(int j = n-1; j > i+1; --j) //lower triangle 0's
        {
            sum += ( U[i][j] * x[j] );
        }

        x[i] = ( y[i] - sum ) / U[i][i];
    }

    cout << "\nx = " << endl;
    printVec(x,n);



    //compute b' = A * x  in order to compare it with b
    double *bp;
    bp = new double[n];
    for( int r=0; r < n; ++r )
    {
        double sum = 0.0;

        for(int c=0; c < n; ++c)
        {
            sum += ( A[r][c] * x[c] );
        }
        bp[r] = sum;
    }
    cout <<"\nb' = " << endl;
    printVec( bp , n );

    //compute two-norm error vector
    double *err;
    err = new double [n];

    for(int r = 0; r < n; ++r)
    {
        err[r] = ( bp[r] - b[r] );
    }

    cout << "\nerror vector: e = b' - b " << endl;
    printVec(err , n);

    //STEP 5 e) compute two-norm error vector
    double sigma = 0.0;
    double evec = 0.0;
    for(int i=0; i < n; ++i)
    {
        sigma += ( abs(err[i])*abs(err[i]) );
    }
    evec = sqrt(sigma);
    cout << "\n two-norm error vector = " << evec << endl;

    return 0;

}
