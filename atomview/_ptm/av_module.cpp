#define PY_SSIZE_T_CLEAN
#include <Python.h>
// Empty module so the PTM library builds as an importable extension; functions are called through ctypes.
static struct PyModuleDef m = {PyModuleDef_HEAD_INIT, "_ptmlib", NULL, -1, NULL};
PyMODINIT_FUNC PyInit__ptmlib(void) { return PyModule_Create(&m); }
