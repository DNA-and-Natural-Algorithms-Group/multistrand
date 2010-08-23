// Compile with 
// g++ -lpython2.6 -o embedding_test embedding_test.cc

#include "../include/python_options.h"
#include <stdio.h>
#include <string.h>

int main(int argc, char **argv) 
{
  char path[300];
  strncpy(path,argv[0],300);
  int i = strlen(argv[0])-2;
  while(i>= 0)
	{
	  if( argv[0][i] == '/' )
		{
		  argv[0][i+1] = '\0';
		  i = -1;
		}
	  i--;
	}
  snprintf( path, 300, "import sys\nsys.path.append(\"%s\")\n", argv[0] );
  //  Py_SetProgramName(argv[0]);
  Py_Initialize();
  //  PySys_SetArgv(argc, argv);
  PyRun_SimpleString(path);
  
  PyObject *options = newObject(options_test, Options);

  PyObject *pyo_str1;
  char *str1 = getStringAttr(options, name, pyo_str1);
  // str1 is new reference, stored in pyo_str1
  printf("String attribute is '%s'\n", str1);
  Py_DECREF( pyo_str1 );
  str1= NULL;
  // str1 is clear

  long num1; 
#ifdef DEBUG_MACROS
  getLongAttr(options, integerddd, &num1);
  // testing debug macros
#endif DEBUG_MACROS
  getLongAttr(options, integer, &num1);
  double num2;
  getDoubleAttr(options, decimal, &num2);
  printf("Number attributes are %d and %f\n", num1, num2);
  
  setDoubleAttr(options, decimal, 4.5);
  getDoubleAttr(options, decimal, &num2);
  printf("Number attributes are %d and %f\n", num1, num2);
  
  char *str2 = getStringItem(getListAttr(options, list_of_strings), 1);
  printf("String from list is '%s'\n", str2);

  // commented out til I work out the deprecated warnings issues.
  //callFunc_NoArgsToNone(options, no_args_no_return);
  //  callFunc_DoubleToNone(options, one_arg_no_return, 2.3);
  
  Py_Finalize();
  return 0;
}

