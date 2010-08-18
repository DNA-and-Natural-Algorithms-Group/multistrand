#include "python_options.h"
#include "optionlists.h"

// TODO: check deallocation of new objects


class stopcomplexes *getStopComplexList(PyObject *options, int index)
{
    int i, j, k, n, m, l, stoptype, count;
    PyObject *stop_conditions, *stop_condition, *tuple_list, *tuple, *cmplx, *strand_list;
    char *tag, *id, *structure;
    class identlist *id_list;
    class complex_item *complexes;
    class stopcomplexes *return_list;
    
    stop_conditions = PyObject_GetAttrString(options, "stopcomplexes");
    
    n = PyList_GET_SIZE(stop_conditions);
    for (i = n - 1; i >= index; i--)
    {
        stop_condition = PyList_GET_ITEM(stop_conditions, i);
        tag = getStringAttr(stop_condition, tag);
        tuple_list = getListAttr(stop_condition, complex_items);
        
        m = PyList_GET_SIZE(tuple_list);
        for (j = m - 1; j >= 0; j--)
        {
            tuple = PyList_GET_ITEM(tuple_list, j);
            cmplx = PyTuple_GET_ITEM(tuple, 0);
            stoptype = getLongItemFromTuple(tuple, 1);
            count = getLongItemFromTuple(tuple, 2);
            structure = getStringAttr(cmplx, structure);
            strand_list = getListAttr(cmplx, strand_list);
            
            l = PyList_GET_SIZE(strand_list);
            id_list = new identlist(getStringAttr(PyList_GetItem(strand_list, l - 1), id), NULL);
            for (k = l - 2; k >= 0; k--)
            {
                id_list = new identlist(getStringAttr(PyList_GetItem(strand_list, k), id), id_list);
            }
            
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









//class stopcomplexes *convertStopcomplexList(PyObject *stop_conditions)
//{
//    int i, j, k, n, m, l, stoptype, count;
//    PyObject *stop_condition, *tuple_list, *tuple, *cmplx, *strand_list;
//    char *tag, *id, *structure;
//    class identlist *id_list;
//    class complex_item *complexes;
//    class stopcomplexes *return_list;
//    
//    n = PyList_GET_SIZE(stop_conditions);
//    for (i = 0; i < n; i++)
//    {
//        stop_condition = PyList_GET_ITEM(stop_conditions, i);
//        tag = getStringAttr(stop_condition, tag);
//        tuple_list = getListAttr(stop_condition, complex_items);
//        
//        m = PyList_GET_SIZE(tuple_list);
//        for (j = 0; j < m; j++)
//        {
//            tuple = PyList_GET_ITEM(tuple_list, j);
//            cmplx = PyTuple_GET_ITEM(tuple, 0);
//            stoptype = getLongItemFromTuple(tuple, 1);
//            count = getLongItemFromTuple(tuple, 2);
//            structure = getStringAttr(cmplx, structure);
//            strand_list = getListAttr(cmplx, strand_list);
//            
//            l = PyList_GET_SIZE(strand_list);
//            id_list = new identlist(getStringAttr(PyList_GetItem(strand_list, l-1), id), NULL);
//            for (k = l-2; k >= 0; k--)
//            {
//                id_list = new identlist(getStringAttr(PyList_GetItem(strand_list, k), id), id_list);
//            }
//            
//            if (j == 0)
//                complexes = new complex_item(structure, id_list, NULL, stoptype, count);
//            else
//                complexes = new complex_item(structure, id_list, complexes, stoptype, count);
//        }
//        
//        if (j == 0)
//            return_list = new stopcomplexes(tag, complexes, NULL);
//        else
//            return_list = new stopcomplexes(tag, complexes, return_list);
//    }
//    
//    return return_list;
//}












