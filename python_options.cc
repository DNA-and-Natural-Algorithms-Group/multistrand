#include "python_options.h"
#include "optionlists.h"


class stopcomplexes *convertStopcomplexList(PyObject *stop_conditions)
{
    int i, n;
    PyObject *stop_condition
    
    n = PyList_Size(stop_conditions);
    for (i = 0; i < n; i++)
    {
        stop_condition = PyList_GetItem(stop_conditions, i);
        
    }
}
