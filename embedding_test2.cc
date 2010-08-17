// Compile with 
// g++ -lpython2.6 -o embedding_test2 embedding_test2.cc

#include "python_options.h"
#include <stdio.h>


int main() 
{
    Py_Initialize();
    PyRun_SimpleString("import sys\nsys.path.append('/research/src/Multistrand-Python/')\n");
    
    PyObject *options = newObject(options_test, Options);
    
    class stopcomplexes *sc = getStopcomplexList(options);
    
    Py_Finalize();
    return 0;
}

