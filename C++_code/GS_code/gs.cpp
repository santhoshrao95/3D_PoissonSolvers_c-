#include<iostream>
#include<math.h>
#include<fstream>
using namespace std;
const int N = 34;
struct ArrClass
{
   double fun[N][N][N], p[N][N][N];
   double x=0.0,y=0.0,z=0.0;
   double h=1.0/(N-1), rcd2= -1.0/6.0;
};
void printAll(ArrClass *var)
{
    for (int j = 0; j < N; j++)
    {
        for (int k = 0; k < N; k++)
        {
            for (int i = 0; i < N; i++)
            {
                cout << "p[" << j<<"]"<<"["<<k<<"]"<<"["<<i<<"]: " << var->p[j][k][i] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
}
void printLevel(ArrClass *var, int level)
{
    ofstream gs("gs_33_5000_j=16.txt",ios::trunc);
    for (int k = 0; k < N; k++)
        {
            for (int i = 0; i < N; i++)
            {
                //cout << var->p[i][0][k] << endl;
                //cout << "p[" << level<<"]"<<"["<<k<<"]"<<"["<<i<<"]: " << var->p[level][k][i] << " ";
                gs << var->p[level][k][i] <<",";

            }
            gs << endl;
        }
        //cout << endl;
    gs.close();
}
void initializeArr(ArrClass *var)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < N; k++)
            {
                var->p[i][j][k] = 0.0;
                var->fun[i][j][k] = 0.0;
            }         
        }       
    }
}
double rhs(double x, double y, double z)
{
    double f_x = 50000.0 * exp(-50.0*(((1-x)*(1-x))+(z*z))) * (100*(((1-x)*(1-x))+(z*z))-2.0);
    return f_x;
}

void calculateFun(ArrClass *var)
{
    var->x=0.0, var->y = 0.0, var->z = 0.0;
    for(int j = 0; j < N; j++)
    {
        var->z=0.0;
        for (int k = 0; k < N; k++)
        {
            var->x = 0.0;
            for (int i = 0; i < N; i++)
            {
                var->fun[j][k][i] = rhs(var->x,var->y,var->z);
                var->x = var->x+var->h;
            }
            var->z = var->z+var->h;   
        }
        var->y = var->y+var->h;
    }    
}
void boundaryValues(ArrClass *var)
{
    var->x = 0.0, var->y = 0.0, var->z = 0.0;
    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < N; i++)
        {
            var->p[j][0][i] = (100.0*var->x)+(500.0*exp(-50.0*((1.0-var->x)*(1.0-var->x))));
            var->x = var->x+var->h;
        }
        var->x = 0.0;    
    }
    var->x = 0.0, var->y = 0.0, var->z = 0.0;
    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < N; i++)
        {
            var->p[j][N-1][i] = 500.0*exp(-50.0*(((1.0-var->x)*(1-var->x))+1.0));
            var->x = var->x+var->h;
        }
        var->x = 0.0;    
    }
    var->x = 0.0, var->y = 0.0, var->z = 0.0;
    for (int j = 0; j < N; j++)
    {
        for (int k = 0; k < N; k++)
        {
            var->p[j][k][0] = 500.0*exp(-50.0*(1.0+(var->z*var->z)));
            var->z = var->z+var->h;
        }
        var->z = 0.0;    
    }
    var->x = 0.0, var->y = 0.0, var->z = 0.0;
    for (int j = 0; j < N; j++)
    {
        for (int k = 0; k < N; k++)
        {
            var->p[j][k][N-1] = (100.0*(1.0-var->z))+(500.0*exp(-50.0*var->z*var->z));
            var->z = var->z+var->h;
        }
        var->z = 0.0;    
    }
}
void gaussSidel(ArrClass *var)
{
    for (int iter = 0; iter< 10000; iter++)
    {
        for (int j = 1; j < N-1; j++)
        {
            for (int k = 1; k < N-1; k++)
            {
                for (int i = 1; i < N-1; i++)
                {
                    var->p[j][k][i] = (var->rcd2*var->h*var->h*var->fun[j][k][i]) - (var->rcd2*(\
                    var->p[j-1][k][i] + var->p[j+1][k][i] + \
                    var->p[j][k-1][i] + var->p[j][k+1][i] + \
                    var->p[j][k][i-1] + var->p[j][k][i+1] ));
                }
          
            }
       
        }   

        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < N; k++)
            {
                for (int i = 0; i < N; i++)
                {
                    if(j==0)
                    {
                        var->p[j][k][i] = var->p[N-2][k][i];
                    }
                    if(j==N-1)
                    {
                        var->p[j][k][i] = var->p[1][k][i];
                    }

                }
          
            }
       
        }         
    }
}
int main()
{
    ArrClass arr;
    initializeArr(&arr);
    calculateFun(&arr);
    boundaryValues(&arr);
    gaussSidel(&arr);
    printLevel(&arr, 16);
    return 0;
}
