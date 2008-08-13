#ifndef _GAP_CLI_ARG_H_
#define _GAP_CLI_ARG_H_

#include <stddef.h>

#define ARG_INT    1
#define ARG_STR    2
#define ARG_IO     3
#define ARG_ARR    4
#define ARG_FLOAT  5
#define ARG_DOUBLE 6
#define ARG_DBL ARG_DOUBLE
#define ARG_OBJ    7

#define offsetofa(type, field) ((int) ((char *) ((type *)0)->field))

typedef struct {
    char *command;	/* What to recognise, including the '-' symbol */
    int type;		/* ARG_??? */
    int value;		/* Set if this argument takes an argument */
    char *def;		/* NULL if non optional argument */
    int offset;		/* Offset into the 'store' address */
} cli_args;

/*
 * Parses command line arguments.
 * 'args' specifies an array cli_args. Each cli_arg contains a default. Each
 * argument with a default of NULL are non optional. This function will return
 * -1 if any of these arguments aren't specified, or if unknown arguments
 * are specified. Returns 0 otherwise.
 *
 * The _args functions set the default values into the store and then run
 * the _config functions, which do the main work of parsing the
 * arguments and filling out the relevant store bits.
 * Ie _args is primarily for argument parsing to functions and _config 
 * is for subsequent configuration of data objects. The _config
 * function ignores the .def section of cli_args.
 */
int gap_parse_args(cli_args *args, void *store, int argc, char **argv);
int gap_parse_config(cli_args *args, void *store, int argc, char **argv);
int gap_parse_obj_args(cli_args *args, void *store, int objc,
		       Tcl_Obj * const objv[]);
int gap_parse_obj_config(cli_args *args, void *store, int objc,
			 Tcl_Obj * const objv[]);

#endif
