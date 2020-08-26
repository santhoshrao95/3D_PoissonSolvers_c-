#include<iostream>
#include<math.h>
#include<fstream>
using namespace std;
const int N = 33;

struct ArrClass
{
   double fun[N][N][N], p[N][N][N], r[N][N][N];
   double x=0.0,y=0.0,z=0.0;
   double h=1.0/(N-1);
};

double rhs(double x, double y, double z)
{
    double f_x = 50000*exp(-50*(((1-x)*(1-x))+(z*z))) * (100*((1-x)*(1-x))+(z*z)-2);
    return f_x;
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
                var->fun[i][j][k]=0.0;
                var->r[i][j][k]=0.0;
                //var->rtemp[i][j][k]=0.0;
               // var->d[i][j][k]=0.0;
                //var->ad[i][j][k]=0.0;
            }
          
        }
       
    }
    
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
                var->fun[j][k][i] = rhs(var->y,var->z,var->x);
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



void calResidual(ArrClass *var)
{
    for (int j = 1; j < N-1; j++)
    {
        for (int k = 1; k < N-1; k++)
        {
            for (int i = 1; i < N-1; i++)
            {
                var->r[j][k][i] = (var->h*var->h*var->fun[j][k][i]) - (\
                    var->p[j-1][k][i] + var->p[j+1][k][i] + \
                    var->p[j][k-1][i] + var->p[j][k+1][i] + \
                    var->p[j][k][i-1] + var->p[j][k][i+1] - \
                    (6.0*var->p[j][k][i]) );
            }
            
        }
        
    }
}

void printLevel(ArrClass *var, int level)
{
    ofstream cg("cg_63_5000_j=31.txt",ios::trunc);
    for (int k = 0; k < N; k++)
        {
            for (int i = 0; i < N; i++)
            {
                //cout << var.p[i][0][k] << endl;
                //cout << "p[" << level<<"]"<<"["<<k<<"]"<<"["<<i<<"]: " << var.p[level][k][i] << " ";
                cg << var->p[level][k][i] <<",";

            }
            cg << endl;
        }
        //cout << endl;
    cg.close();

}

void CGMethod()
{

    double d[N][N][N], rtemp[N][N][N], ad[N][N][N], rold[N][N][N], \
    dnom=0.0, alpha = 0.0, beta = 0.0, l2norm = 0.0; 
    int iter = 0.0;
    double num2=0.0, num1=0.0;
    ArrClass arr;
    
    //initialzie local arrays
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < N; k++)
            {
                d[i][j][k] = 0.0;
                rtemp[i][j][k] =0.0;
                rold[i][j][k] = 0.0;
                ad[i][j][k] = 0.0;
            }
            
        }
        
    }


    initializeArr(&arr);
    calculateFun(&arr);
    boundaryValues(&arr);
    calResidual(&arr); 
    
    for (int j = 1; j < N-1; j++) //copy r to d and rtemp
    {
        for (int k = 1; k  < N-1; k++)
        {
            for (int i = 1; i < N-1; i++)
            {
                d[j][k][i] = arr.r[j][k][i];
                rtemp[j][k][i] = arr.r[j][k][i];
            }
            
        }
        
    }

    
    for (int j = 1; j < N-1; j++) // calculate num2 using rtemp
    {
        for (int k = 1; k < N-1; k++)
        {
            for (int i = 1; i < N-1; i++)
            {
                num2 = num2 + (rtemp[j][k][i]*rtemp[j][k][i]);
            }
            
        }
        
    }
    
    num1=num2;


    do
    {
        for (int j = 1; j < N-1; j++) //calculate ad
        {
            for (int k = 1; k < N-1; k++)
            {
                for (int i = 1; i < N-1; i++)
                {
                    ad[j][k][i] =  d[j-1][k][i] + d[j+1][k][i] + \
                    d[j][k-1][i] + d[j][k+1][i] + \
                    d[j][k][i-1] + d[j][k][i+1] - \
                    (6.0*d[j][k][i]) ;
                }
                
            }
            
        }

        dnom = 0.0;

        for (int j = 1; j < N-1; j++) //calculate dnom=d*ad
        {
            for (int k = 1; k < N-1; k++)
            {
                for (int i = 1; i < N-1; i++)
                {
                    dnom = dnom+(d[j][k][i]*ad[j][k][i]);
                }
                
            }
            
        }

        alpha = num2/dnom; //rSquare divided by (d*ad)

        for (int j = 1; j < N-1; j++) //update p and r
        {
            for (int k = 1; k < N-1; k++)
            {
                for (int i = 1; i < N-1; i++)
                {
                    arr.p[j][k][i] = arr.p[j][k][i] + (alpha*d[j][k][i]);
                    arr.r[j][k][i] = arr.r[j][k][i] - (alpha*ad[j][k][i]);
                }
                
            }
            
        }

    
 




        num2 = 0.0;
        for (int j = 1; j < N-1; j++) // update num2
        {
            for (int k = 1; k < N-1; k++)
            {
                for (int i = 1; i < N-1; i++)
                {
                    num2 = num2+(arr.r[j][k][i]*arr.r[j][k][i]);
                }
                
            }
            
        }

        beta = num2/num1;




        for (int j = 1; j < N-1; j++) // update d
        {
            for (int k = 1; k < N-1; k++)
            {
                for (int i = 1; i < N-1; i++)
                {
                    d[j][k][i] = arr.r[j][k][i] + (beta*d[j][k][i]);
                }
                
            }
            
        }


        for (int j = 0; j <= N-1; j++) //applying periodic boundary condition on d
        {
            for (int k = 0; k <= N-1; k++)
            {
                for (int i = 1; i <= N-1; i++)
                {
                    if(j==0)
                    {
                        arr.p[j][k][i] = arr.p[N-2][k][i];
                    }
                    if(j==N-1)
                    {
                        arr.p[j][k][i] = arr.p[1][k][i];
                    }

                }
          
            }
       
        }

    for (int j = 1; j < N-1; j++) // update rtemp
    {
        for (int k = 1; k  < N-1; k++)
        {
            for (int i = 1; i < N-1; i++)
            {
                rtemp[j][k][i] = arr.r[j][k][i];
            }
            
        }
        
    }

    num1=num2;


    iter++;
    cout << iter << endl;

    for (int j = 0; j < N; j++)
    {
        for (int k = 0; k < N; k++)
        {
            for (int i = 0; i < N; i++)
            {
                rold[j][k][i] = arr.r[j][k][i];
            }
            
        }
        
    }


    double err = 0.0;
    for (int j = 1; j < N-1; j++) // calculate err
    {
        for (int k = 1; k < N-1; k++)
        {
            for (int i = 1; i < N-1; i++)
            {
                err = err+((arr.r[j][k][i]*arr.r[j][k][i]) - (rold[j][k][i]*rold[j][k][i]));
            }
                
        }
            
    }    

    l2norm = sqrt(err);
    


    } while (l2norm>0.0000001);
    

printLevel(&arr, 16);

    



}

int main()
{
    CGMethod();

}