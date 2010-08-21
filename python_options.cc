#include "python_options.h"
#include "optionlists.h"

// TODO: check deallocation of new objects





class stopcomplexes *getStopComplexList(PyObject *options, int index)
{
    class identlist *id_list;
    class complex_item *complexes;
    class stopcomplexes *return_list;
    
    PyObject *strand_list;
    char *structure;
    PyObject *pyo_structure;
    int count;
    int stoptype;
    PyObject *cmplx;
    PyObject *tuple;
    PyObject *stop_condition;
    PyObject *pyo_tag;
    char *tag;
    PyObject *tuple_list;
    PyObject *stop_conditions;

    stop_conditions= PyObject_GetAttrString(options, "stop_conditions");
    // new reference
    
    int n = PyList_GET_SIZE(stop_conditions);
    for (int i = n - 1; i >= index; i--)
    {
        stop_condition = PyList_GET_ITEM(stop_conditions, i);
        pyo_tag = NULL;

        tag = getStringAttr(stop_condition, tag, pyo_tag);
        // new reference is stored in tmp_tag.

        tuple_list = getListAttr(stop_condition, complex_items);
        // new reference here
        
        int m = PyList_GET_SIZE(tuple_list);
        for (int j = m - 1; j >= 0; j--)
        {
            tuple = PyList_GET_ITEM(tuple_list, j);
            cmplx = PyTuple_GET_ITEM(tuple, 0);
            stoptype = getLongItemFromTuple(tuple, 1);
            count = getLongItemFromTuple(tuple, 2);
            // none of these steal the reference, only borrow.

            structure = getStringAttr(cmplx, structure, pyo_structure);

            strand_list = getListAttr(cmplx, strand_list);
            // new reference

            id_list = makeID_list(strand_list);
            // does not steal the reference to strand_list.

            Py_DECREF( strand_list );
            // strand_list reference is now clean.
            
            if (j == m - 1)
                complexes = new complex_item(structure, id_list, NULL, stoptype, count);
            else
                complexes = new complex_item(structure, id_list, complexes, stoptype, count);

            Py_DECREF( pyo_structure );
            structure = NULL;
            // pyo_structure is now clean again.
        }
        
        Py_DECREF( tuple_list );
        // tuple_list is now clean.

        if (i == n - 1)
            return_list = new stopcomplexes(tag, complexes, NULL);
        else
            return_list = new stopcomplexes(tag, complexes, return_list);
        
        Py_DECREF( pyo_tag );
        tag = NULL;
    }

    Py_DECREF( stop_condition );
    return return_list;
}


class identlist *getID_list(PyObject *options, int index)
{
  PyObject *pyo_strand_list = NULL;
  PyObject *pyo_start_complex_list = NULL;
  pyo_start_complex_list = PyObject_GetAttrString( options, "start_state");
  // new reference

  pyo_strand_list = PyObject_GetAttrString(PyList_GET_ITEM( pyo_start_complex_list, index ));
  // new reference

  Py_DECREF( pyo_start_complex_list );
  // pyo_start_complex_list is now clean, pyo_strand_list has a new reference.

  class identlist *tmp = NULL;

  tmp = makeID_list( pyo_strand_list );
  // makeID_list does not steal the reference. 

  Py_DECREF( pyo_strand_list );
  // pyo_strand_list is now clean
}

// TODO still.
class identlist *makeID_list(PyObject *strand_list)
{
    class identlist *id_list;
    Py_INCREF( strand_list );

    int n = PyList_GET_SIZE(strand_list);
    
    PyObject* py_strand = PyList_GetItem(strand_list, n - 1);
    // borrowed reference, strand_list is passed to us by 
    long id = getLongAttr(py_strand, id);
    char *name;
    PyObject *pyo_name;

    name = getStringAttr(py_strand, name, pyo_name, NULL);
    // new reference
    
    id_list = new identlist( id, name, NULL );
    // id_list should make a copy of name, it is now ok to free it.

    Py_DECREF( pyo_name );
    name = NULL;
    // pyo_name is now clear.

    for (int i = n - 2; i >= 0; i--)
    {
        py_strand = PyList_GetItem(strand_list, i);
        name = getStringAttr(py_strand, name, pyo_name);
        // new reference
        
        id_list = new identlist( id, name, pyo_name, id_list );
        // id_list should make a copy of name, it is now ok to free it.
        
        Py_DECREF( pyo_name );
        name = NULL;
        // pyo_name is now clear.
    }
    
    return id_list;
}








