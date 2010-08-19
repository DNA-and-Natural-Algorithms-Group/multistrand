// Compile with 
// g++ -g -lpython2.6 -o embedding_test2 optionlists.o embedding_test2.cc python_options.cc

#include "python_options.h"
#include <stdio.h>


int main() 
{
    Py_Initialize();
    PyRun_SimpleString("import sys\nsys.path.append('/research/src/Multistrand-Python/')\n");
    
    PyObject *options = newObject(options_test, Options);
    
    int a, b;
    getLongAttr( options, integer, &a);
    getLongAttr( options, neg_integer, &b);
    printf("The positive integer was %d.\n", a);
    printf("The negative integer was %d.\n", b);
    
    if (testLongAttr( options, integer, >, 0))
        printf("The positive integer is greater than 0.\n");
    else
        printf("The positive integer is not greater than 0.\n");
    
    if (testLongAttr( options, neg_integer, >, 0))
        printf("The negative integer is greater than 0.\n");
    else
        printf("The negative integer is not greater than 0.\n");
    
    class stopcomplexes *sc = getStopComplexList(options, 0);
    class identlist *il = getID_list(options, 1);
    
    Py_Finalize();
    return 0;
}

