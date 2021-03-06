#ifndef _TG_IO_LOW_H_
#define _TG_IO_LOW_H_

#include "tg_cache_item.h"

/*
 * The iface concept is to provide a standardised interface for accessing
 * (both read and write) the primary database objects; contigs, sequences,
 * annotations, etc.
 *
 * The purpose of this is two fold.
 *
 * 1) It allows us to replace the underlying database format with something
 *    different (eg g-library, hdf5, SQL, etc) or in-memory faked systems
 *    such as the contig editor.
 *
 * 2) It allows us to implement I/O agnostic algorithms. An example is the
 *    consensus algorithm which can operate directly on the disk structures,
 *    from within a portion of the contig editor's in memory structures,
 *    or on a hybrid "overlay" of fake sequences and real sequences as used
 *    in the prefinish code.
 *
 * Although this is initially based on the g-library, we try to acknowledge
 * how other systems may work (such as multiple tables in Oracle vs one
 * big shared table in dbm or g-lib).
 *
 * Hence every primary object has its own set of allocate, read and write
 * methods. A "contig" number therefore refers to the number returned from
 * a contig allocate() call. It should no longer be assumed that all numbers
 * count from 1 to N inclusive with no gaps. 
 */


typedef GCardinal GRec;

/*
 * Common interfaces to be used in every primary gap5 object.
 */
#define STANDARD_IFACE \
    /* Allocate and deallocate. Init is allocate() + lock + initialise */ \
    GRec (*create)(void *dbh, void *from);				  \
    int (*destroy)(void *dbh, GRec rec);				  \
									  \
    /* Locking and unlocking */						  \
    GView (*lock)(void *dbh, GRec rec, int mode);			  \
    int (*unlock)(void *dbh, GView view);				  \
    int (*upgrade)(void *dbh, GView view, int mode);			  \
    int (*abandon)(void *dbh, GView view);				  \
									  \
    /* Read/Write */				                          \
    cached_item *(*read)(void *dbh, GRec rec);	                          \
    int (*write)(void *dbh, cached_item *ci);                             \
									  \
    /* Queries on size */						  \
    int (*info)(void *dbh, GView view, GViewInfo *vi);


/*
 * Primary Gap5 data types.
 * Each of these has the common set of interface pointers, in addition
 * to local methods.
 */
typedef struct {
    STANDARD_IFACE
} io_array;

typedef struct {
    STANDARD_IFACE
} io_database;

typedef struct {
    STANDARD_IFACE
    GRec (*index_query)(void *dbh, char *name);
    int  (*index_add)(void *dbh, char *name, GRec rec);
} io_contig;

typedef struct {
    STANDARD_IFACE
    GRec (*index_query)(void *dbh, char *name);
    int  (*index_add)(void *dbh, char *name, GRec rec);
} io_seq;

typedef struct {
    STANDARD_IFACE
} io_anno;

typedef struct {
    STANDARD_IFACE
} io_anno_ele;

typedef struct {
    STANDARD_IFACE
} io_dnasrc;

typedef struct {
    STANDARD_IFACE
} io_vector;

typedef struct {
    STANDARD_IFACE
} io_bin;

typedef struct {
    STANDARD_IFACE
} io_track;

typedef struct {
    /* Higher level database-level functions */
    int (*create)(char *dbname);
    void *(*connect)(char *dbname, int ro);
    int (*disconnect)(void *dbh);
    int (*commit)(void *dbh);
    int (*lock)(void *dbh);
    int (*unlock)(void *dbh);

    /* The objects themselves */
    io_array        array; /* generic array */
    io_database     database;
    io_contig       contig;
    io_bin          bin;
    io_track        track;
    io_seq          seq;
    io_anno_ele     anno_ele;
    io_anno         anno;
    io_dnasrc       dnasrc;
    io_vector       vector;
} iface;

#endif /* _TG_IO_LOW_H_ */
