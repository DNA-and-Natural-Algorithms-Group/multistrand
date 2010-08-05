// Compile with 
// g++ -lpython2.6 -o embedding_test embedding_test.cc

#include "options-python.h"
#include <stdio.h>


int main() 
{
    Py_Initialize();
    PyRun_SimpleString("import sys\nsys.path.append('/research/src/Multistrand-Python/')\n");
    
    PyObject *options = newObject(options_test, Options);
    
    char *str1 = getStringAttr(options, name);
    printf("String attribute is '%s'\n", str1);
    
    int num1 = getLongAttr(options, integer);
    double num2 = getDoubleAttr(options, decimal);
    printf("Number attributes are %d and %f\n", num1, num2);
    
    setDoubleAttr(options, decimal, 4.5);
    num2 = getDoubleAttr(options, decimal);
    printf("Number attributes are %d and %f\n", num1, num2);
        
    char *str2 = getStringItem(getListAttr(options, list_of_strings), 1);
    printf("String from list is '%s'\n", str2);
    
    callFunc_NoArgsToNone(options, no_args_no_return);
    callFunc_DoubleToNone(options, one_arg_no_return, 2.3);
    
    Py_Finalize();
    return 0;
}

