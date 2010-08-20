#include "python_options.h"
#include "optionlists.h"

// TODO: check deallocation of new objects





class stopcomplexes *getStopComplexList(PyObject *options, int index)
{
    class identlist *id_list;
    class complex_item *complexes;
    class stopcomplexes *return_list;
    
    PyObject *stop_conditions = PyObject_GetAttrString(options, "stop_conditions");
    
    int n = PyList_GET_SIZE(stop_conditions);
    for (int i = n - 1; i >= index; i--)
    {
        PyObject *stop_condition = PyList_GET_ITEM(stop_conditions, i);
        char *tag = getStringAttr(stop_condition, tag);
        PyObject *tuple_list = getListAttr(stop_condition, complex_items);
        
        int m = PyList_GET_SIZE(tuple_list);
        for (int j = m - 1; j >= 0; j--)
        {
            PyObject *tuple = PyList_GET_ITEM(tuple_list, j);
            PyObject *cmplx = PyTuple_GET_ITEM(tuple, 0);
            int stoptype = getLongItemFromTuple(tuple, 1);
            int count = getLongItemFromTuple(tuple, 2);
            char *structure = getStringAttr(cmplx, structure);
            PyObject *strand_list = getListAttr(cmplx, strand_list);
            
            id_list = makeID_list(strand_list);
            
            if (j == m - 1)
                complexes = new complex_item(structure, id_list, NULL, stoptype, count);
            else
                complexes = new complex_item(structure, id_list, complexes, stoptype, count);
        }
        
        if (i == n - 1)
            return_list = new stopcomplexes(tag, complexes, NULL);
        else
            return_list = new stopcomplexes(tag, complexes, return_list);
    }
    
    return return_list;
}


class identlist *getID_list(PyObject *options, int index)
{
    return makeID_list(PyObject_GetAttrString(PyList_GET_ITEM(PyObject_GetAttrString(options, "start_complexes"), index), "strand_list"));
}


class identlist *makeID_list(PyObject *strand_list)
{
    class identlist *id_list;
    
    int n = PyList_GET_SIZE(strand_list);
    
    PyObject* py_strand = PyList_GetItem(strand_list, n - 1);
    id_list = new identlist( getLongAttr(py_strand, id), getStringAttr(py_strand, name), NULL );
    
    for (int i = n - 2; i >= 0; i--)
    {
        py_strand = PyList_GetItem(strand_list, i);
        id_list = new identlist( getLongAttr(py_strand, id), getStringAttr(py_strand, name), id_list );
    }
    
    return id_list;
}








