// Compile with 
// g++ -lpython2.6 -o embedding_test embedding_test.cc

#include "options-python.h"
#include <stdio.h>


int main() 
{
    Py_Initialize();
    
    
    // This stuff works:
    // PyRun_SimpleString("import random\n");
    // PyRun_SimpleString("who = random.choice(['world', 'y\\'all', 'fellas'])\nprint 'hello %s' % who\n");
    
    
    // But this stuff doesn't work:
    PyObject *module, *type, *options;
    
    PyRun_SimpleString("import sys\nsys.path.append('/research/src/Multistrand-Python/')\n");
    //PyRun_SimpleString("import options_test");
    
    //module = PyImport_ImportModule("options_test");
    //type = PyObject_GetAttrString(module, "Options");
    //options = PyObject_CallObject(type, NULL);
    
    
    // This is what I would like to use:
    options = PyObject_CallObject(PyObject_GetAttrString(PyImport_ImportModule("options_test"), "Options"), NULL);
    
    char *str1 = getStringAttr(options, name);
    printf("String attribute is '%s'\n", str1);
    
    int num1 = getLongAttr(options, integer);
    double num2 = getDoubleAttr(options, decimal);
    printf("Number attributes are %d and %f\n", num1, num2);
    
    setDoubleAttr(options, decimal, 4.5);
    printf("Number attributes are %d and %f\n", num1, num2);
        
    char *str2 = getStringItem(getListAttr(options, list_of_strings), 0);
    printf("String from list is %s\n", str1);
    
    //callFunc_NoArgsToNone(options, no_args_no_return);
    //callFunc_DoubleToNone(options, one_arg_no_return, 2.3);
    
    //Py_DECREF(module_name);
    //Py_DECREF(module);
    //Py_DECREF(type);
    //Py_DECREF(options);
    Py_Finalize();
    return 0;
}

/*
int main(int argc, char *argv[])
{
    PyObject *pName, *pModule, *pDict, *pFunc;
    PyObject *pArgs, *pValue;
    int i;

    if (argc < 3) {
        fprintf(stderr,"Usage: call pythonfile funcname [args]\n");
        return 1;
    }

    Py_Initialize();
    pName = PyString_FromString(argv[1]);
    // Error checking of pName left out

    pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule != NULL) {
        pFunc = PyObject_GetAttrString(pModule, argv[2]);
        // pFunc is a new reference

        if (pFunc && PyCallable_Check(pFunc)) {
            pArgs = PyTuple_New(argc - 3);
            for (i = 0; i < argc - 3; ++i) {
                pValue = PyInt_FromLong(atoi(argv[i + 3]));
                if (!pValue) {
                    Py_DECREF(pArgs);
                    Py_DECREF(pModule);
                    fprintf(stderr, "Cannot convert argument\n");
                    return 1;
                }
                // pValue reference stolen here: 
                PyTuple_SetItem(pArgs, i, pValue);
            }
            pValue = PyObject_CallObject(pFunc, pArgs);
            Py_DECREF(pArgs);
            if (pValue != NULL) {
                printf("Result of call: %ld\n", PyInt_AsLong(pValue));
                Py_DECREF(pValue);
            }
            else {
                Py_DECREF(pFunc);
                Py_DECREF(pModule);
                PyErr_Print();
                fprintf(stderr,"Call failed\n");
                return 1;
            }
        }
        else {
            if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "Cannot find function \"%s\"\n", argv[2]);
        }
        Py_XDECREF(pFunc);
        Py_DECREF(pModule);
    }
    else {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"%s\"\n", argv[1]);
        return 1;
    }
    Py_Finalize();
    return 0;
}
*/
