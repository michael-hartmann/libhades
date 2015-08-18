#ifndef NPY_DICT
#define NPY_DICT

const char *npy_dict_get_value(const char *dict, const char *key);

int npy_dict_get_fortran_order(const char *dict, int *fortran_order);
int npy_dict_get_shape(const char *dict, int *rows, int *columns);
int npy_dict_get_descr(const char *dict, int *is_complex);

#endif
