#include <stdio.h>
#include <QC3_interface.hpp>
#include <stdlib.h>

#include "fileutils.h"
#include "pylifecycle.h"
#include "pymem.h"

int main(int argc, char **argv)
{
    double L = 10.0;
    double nnP[3] ={0.,0.,0.};
    double Estart = 3.0033;
    double Eend = 3.0040;
    int steps = 4;
   
    
    wchar_t *program = Py_DecodeLocale(argv[0], NULL);
    if (program == NULL) {
        fprintf(stderr, "Fatal error: cannot decode argv[0]\n");
        exit(1);
    }  
    Py_SetProgramName(program);  /* optional but recommended */
        
    Py_Initialize();
//     PyRun_SimpleString("print('Hello World from Embedded Python!!!')");
    
    
    double Pkcot[2]={1,0};
    double Pkiso[2]={1,0};
    int Nkcot=2;
    int Nkiso=2;
    
    
    for(int i=0;i<10;i++){
        printf("starting iteration %d\n",i);
//         double r=test_python(i);
         double r=python_detQC(Estart, Eend, steps,  L,  nnP, Nkcot,Pkcot,Nkiso, Pkiso);
         printf("solution is: %g\n",r);
    }
    
    
    if (Py_FinalizeEx() < 0) {
        exit(120);
    }
    PyMem_RawFree(program);
//     printf("\nPress any key to exit...\n");
//     if(!_getch()) _getch();
    return 0;
}
