// Compile with 
// g++ -lpython2.6 -o embedding_test embedding_test.cc

#include "options-python.h"
#include <stdio.h>

int main() 
{
    Py_Initialize();
    
    PyRun_SimpleString("import random\nwho = random.choice(['world', 'y\\'all', 'fellas'])\nprint 'hello %s' % who\n");
    
    Py_Finalize();
    return 0;
}

