// Compile with 
// g++ -lpython2.6 -o embedding_test embedding_test.cc

#include "options-python.h"
#include <stdio.h>

int main() 
{
    Py_Initialize();
    
    
    /* This stuff works: */
    // PyRun_SimpleString("import random\n");
    // PyRun_SimpleString("who = random.choice(['world', 'y\\'all', 'fellas'])\nprint 'hello %s' % who\n");
    
    
    /* But this stuff doesn't work: */
    PyObject *module, *type, *options;
    
    PyRun_SimpleString("import sys\nsys.path.append('/research/src/Mulistrand-Python')\n");
    PyRun_SimpleString("import options_test");
    
    //module = PyImport_ImportModule("options_test");
    //type = PyObject_GetAttrString(module, "Options");
    //options = PyObject_CallObject(type, NULL);
    
    
    /* This is what I would like to use: */
    //options = PyObject_CallObject(PyObject_GetAttrString(PyImport_Import(PyString_FromString("options_test.py")), "Options"), NULL);
    
    
    //Py_DECREF(module_name);
    //Py_DECREF(module);
    //Py_DECREF(type);
    //Py_DECREF(options);
    Py_Finalize();
    return 0;
}

