#include "python_options.h"
#include "optionlists.h"

// TODO: check deallocation of new objects


class stopcomplexes *convertStopcomplexList(PyObject *stop_conditions)
{
    int i, j, k, n, m, l, stoptype, count;
    PyObject *stop_condition, *tuple_list, *tuple, *complex, *strand;
    char *tag, *id;
    class identlist *id_list;
    class complex_item *complexes;
    class stopcomplexes *return_list;
    
    n = PyList_GET_SIZE(stop_conditions);
    for (i = 0; i < n; i++)
    {
        stop_condition = PyList_GET_ITEM(stop_conditions, i);
        tag = getStringAttr(stop_condition, tag);
        tuple_list = getListAttr(stop_condition, complex_items);
        
        m = PyList_GET_SIZE(tuple_list);
        for (j = 0; j < m; j++)
        {
            tuple = PyList_GET_ITEM(tuple_list, j);
            complex = PyTuple_GET_ITEM(tuple, 0);
            stoptype = getLongItemFromTupletuple, 1);
            count = getLongItemFromTuple(tuple, 2);
            structure = getStringAttr(complex, structure);
            strand_list = getListAttr(complex, strand_list);
            
            l = PyList_GET_SIZE(strand_list);
            id_list = new identlist(getStringAttr(PyList_GetItem(strand_list, l-1), id), NULL);
            for (k = l-2; k >= 0; k--)
            {
                id_list = new identlist(getStringAttr(PyList_GetItem(strand_list, k), id), id_list);
            }
            
            if (j == 0)
                complexes = new complex_item(structure, strand_list, NULL, stoptype, count);
            else
                complexes = new complex_item(structure, strand_list, complexes, stoptype, count);
        }
        
        if (j == 0)
            return_list = new complex_item(tag, complexes, NULL);
        else
            return_list = new complex_item(tag, complexes, return_list);
    }
    
    return return_list
}
