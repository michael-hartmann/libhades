#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define WHITESPACE " \r\n\t"

#include <libhades/parse_npy_dict.h>

const char *npy_dict_get_value(const char *dict, const char *key)
{
    int delta;
    const char *p;

    if((p = strstr(dict, key)) == NULL)
        return NULL;

    /* format: KEY\s*:\s* */
    p += strlen(key);
    p += strspn(p, WHITESPACE);
    delta = strspn(p, ":");
    if(delta == 0)
        return NULL;
    p += delta;
    p += strspn(p, WHITESPACE);

    return p;
}

int npy_dict_get_fortran_order(const char *dict, int *fortran_order)
{
    const char *p;

    p = npy_dict_get_value(dict, "\"fortran_order\"");
    if(p == NULL)
        p = npy_dict_get_value(dict, "'fortran_order'");
    if(p == NULL)
        return -1;

    if(strstr(p, "True"))
    {
        *fortran_order = 1;
        return 0;
    }
    if(strstr(p, "False"))
    {
        *fortran_order = 0;
        return 0;
    }

    return -2;
}

int npy_dict_get_shape(const char *dict, int *rows, int *columns)
{
    const char *p;

    p = npy_dict_get_value(dict, "\"shape\"");
    if(p == NULL)
        p = npy_dict_get_value(dict, "'shape'");
    if(p == NULL)
        return -1;

    if(*p != '(')
        return -2;
    
    *columns = 1;

    p++;
    p += strspn(p, WHITESPACE);

    if(!isdigit(*p))
        return -3;

    *rows = atoi(p);

    p = strchr(p, ',');
    if(p == NULL)
        return -4;

    p++;
    p += strspn(p, WHITESPACE);

    if(!isdigit(*p))
        *columns = 1;
    else
        *columns = atoi(p);

    return 0;
}

int npy_dict_get_descr(const char *dict, int *is_complex)
{
    const char *p;

    p = npy_dict_get_value(dict, "\"descr\"");
    if(p == NULL)
        p = npy_dict_get_value(dict, "'descr'");
    if(p == NULL)
        return -1;

    if(strncmp(p+1, "<c16", 4) == 0)
    {
        *is_complex = 1;
        return 0;
    }

    if(strncmp(p+1, "<d8", 3) == 0)
    {
        *is_complex = 0;
        return 0;
    }

    return -2;
}
