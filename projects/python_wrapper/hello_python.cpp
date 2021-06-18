#include <stdio.h>
// #include <conio.h>
#include <pyhelper.hpp>

int main()
{
    double L = 10.0;
    double nnP[3] ={0.,0.,0.};
    double Estart = 3.0033;
    double Eend = 3.0040;
    int steps = 40;
    
    CPyInstance pyInstance;
    
    PyRun_SimpleString("print('Hello World from Embedded Python!!!')");
    double r=python_detQC(Estart, Eend, steps,  L,  nnP);
    printf("solution is: %g\n",r);
//     printf("\nPress any key to exit...\n");
//     if(!_getch()) _getch();
    return 0;
}
