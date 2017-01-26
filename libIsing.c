#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

/*typedef struct LatticeData
{
    int* lattice;
    double h;
    double J;
    double beta;
    int L;
    int R;
} LatticeData;
typedef double (*funcPtr_t)(LatticeData);
*/
//=========================================================
// Slave functions for "settle"
//=========================================================
void fixIndex(int* i, int L)
{
    if (*i < 0)
    {
        *i += L;
    }
    else if (*i >= L)
    {
        *i -= L;
    }
    return;
}

int idx(int i, int j, int L)
{
    fixIndex(&i, L);
    fixIndex(&j, L);
    return (i*2*L + 2*j);
}

// This is entirely for debugging purposes
void printLattice(int* lattice, int L)
{
    printf("\n----------------------------------");
    for (int i = 0; i < L; i++)
    {
        printf("\n");
        for (int j = 0; j < L; j++)
        {
            printf("%2d ", lattice[idx(i,j,L)]);
        }
    }
    printf("\n----------------------------------\n");
    
    return;
}

extern "C" int getM(int* lattice, int L)
{   
    int M = 0;
    
    //printLattice(lattice, L);
    
    for(int i=0;i<L;i++)
    {
        for(int j=0; j<L; j++)
        {
            M += lattice[idx(i,j,L)];
        }
    }
    return M;
}

extern "C" int getMForSpin(int* lattice, int L, int R, int i, int j)
{
    int M = 0;
    
    for (int di=-R; di<R; di++)
    {
        for (int dj=-R; dj<R; dj++)
        {
            M += lattice[idx(i+di,j+dj,L)];
        }
    }
    
    return M;
}

extern "C" void getMLattice(int* lattice, int* mLatt, int L, int R)
{
    for (int i = 0; i<L; i++)
    {
        for (int j=0; j<L; j++)
        {
            mLatt[idx(i,j,L)] = getMForSpin(lattice, L, R, i, j);
        }
    }
}

extern "C" int getInterSpin(int* lattice, int L, int R, int i, int j)
{
    int E_ij = 0;
    
    for (int di=-R; di<R; di++)
    {
        for (int dj=-R; dj<R; dj++)
        {
            if (di!=0 || dj != 0) // if they're not both zero
            {
                E_ij += lattice[idx(i,j,L)]*lattice[idx(i+di, j+dj,L)];
            }
        }
    }
    
    return E_ij;
}

extern "C" int getInterVal(int* lattice, int L, int R)
{
    int E = 0;
    
    for (int i=0; i<L; i++)
    {
        for (int j=0; j<L; j++)
        {
            E += getInterSpin(lattice, L, R, i, j);
        }
    }
    
    return E;
}

extern "C" int getInterLattice(int* lattice, int* iLattice, int L, int R, double h)
{
    int E_ij = 0;
    int h_sgn = 1;
    
    if signbit(h)
    {
        h_sgn = -1;
    }
    
    for (int i=0; i<L; i++)
    {
        for (int j=0; j<L; j++)
        {
            E_ij = getInterSpin(lattice,L,R,i,j);
            iLattice[idx(i,j,L)] = E_ij;
        }
    }
}

/*double getm(int* lattice, int L)
{
    printf("In the getm function\n");
    int M = getM(lattice, L);
    double A = pow(L,2);
    double m = M/A;
    
    return m;
}*/
/*
void parseFuncNames(char* trackFuncNames, int tfnLen, funcPtr_t* funcList)//(**func_list)(double*, int))
{
    char letter = ' ';
    char funcName[tfnLen];
    int i = 0;
    int j = 0;
    while (i < tfnLen)
    {
        while ((letter != ';') && (letter != '\0'))
        {
            letter = trackFuncNames[i];
            funcName[i-j] = letter;
            i ++;
        }
        funcName[i-j] = '\0';
        if (strcmp(funcName, "getM")==0)
        {
            funcList[j] = &getM;
        }
        else if (strcmp(funcName,"getm")==0)
        {
            funcList[j] = &getm;
        }
        j ++;
    }
}
*/
/* Returns an integer in the range [0, n).
 *
 * Uses rand(), and so is affected-by/affects the same seed.
 */
int randint(int n) 
{
    if ((n - 1) == RAND_MAX) 
    {
        return rand();
    } 
    else 
    {
        // Chop off all of the values that would cause skew...
        long end = RAND_MAX / n; // truncate skew
        assert (end > 0L);
        end *= n;
        
        // ... and ignore results from rand() that fall above that limit.
        // (Worst case the loop condition should succeed 50% of the time,
        // so we can expect to bail out of this loop pretty quickly.)
        int r;
        while ((r = rand()) >= end);

        return r % n;
    }
}

double calcDE(int* lattice, double h, double J, int L, int R, int i, int j)
{
    int q = pow(2*R + 1, 2) - 1;
    double E_by_s = -h;
    for (int di = -R; di < R+1; di++)
    {
        for (int dj = -R; (dj < R+1); dj++)
        {
            if (di != 0 || dj != 0)
            {
                //printf("%d %d\n", is, js);
                E_by_s -= J*lattice[idx(i+di,j+dj,L)]/q;
            }
        }
    }
    // dE = E_flip - E_same = -2*E_same = -2*s[i,j]*E_by_s
    return -2*lattice[idx(i,j,L)]*E_by_s;
}


//=========================================================
// One of the primary functions: settle
//=========================================================
extern "C" void settle(int* lattice, int L, int R, double h, double J, 
                       double beta, int randSeed, int N_iter)
{
    //printLattice(lattice, L);
    //int* lattice = (int*) latticev;
    //double* trackVals = (double*) trackValsv;
    int i,j;
    double dE, P_flip;
    int Np = 100000;
    
    //printf("%f\n", h);
    
    //LatticeData ld = {lattice, h, J, beta, L, R};
    
    // set the random seed
    srand(randSeed);
    
    // Begin the main loop
    for (int n = 0; n < N_iter; n++)
    {
        i = randint(L);
        j = randint(L);
        /*if (i > ld.L/2 && j > ld.L/2)
        {
            printf("(i,j) = (%d, %d) %d\n", i, j, lattice[idx(i,j,L)]);
        }*/
        
        // Get dE = E_flip - E_same = -2*E_same
        dE = calcDE(lattice, h, J, L, R, i, j);
        P_flip = 1/(1 + exp(beta*dE));
        if ((dE < 0) || (randint(Np) < Np*P_flip))
        {
            lattice[idx(i,j,L)] = -lattice[idx(i,j,L)];
            //printf("%2d, %-2.2f, %2.2f\n", -lattice[idx(i,j,L)], dE, P_flip);
        }
        
    }
    //printLattice(lattice, L);
}

