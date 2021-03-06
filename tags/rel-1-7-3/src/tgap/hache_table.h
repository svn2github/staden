#ifndef _HACHE_TABLE_H_
#define _HACHE_TABLE_H_

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

/* The data referenced by the hash table */
typedef union {
    uint64_t i;
    void *p;
} HacheData;

/* A hash item with "next" pointer to use in a linked list */
typedef struct HacheItemStruct {
    struct HacheItemStruct *next;
    HacheData data;        /* user defined data attached to this key */
    char    *key;         /* key we hashed on */
    int      key_len;     /* and its length */
    int      order;	  /* Index into HacheTable->ordering[] */
    int	     ref_count;	  /* Starts as 1 */
} HacheItem;

typedef struct HacheOrderStruct {
    HacheItem *hi;
    int next;
    int prev;
} HacheOrder;

/* The main hash table structure itself */
typedef struct {
    int	      cache_size; /* Max size before purging cached entries */
    int       options;  /* HASH_FUNC & HASH_OPT macros */
    uint32_t  nbuckets; /* Number of hash buckets; power of 2 */
    uint32_t  mask;	/* bit-mask equiv of nbuckets */
    int       nused;    /* How many hash entries we're storing */
    HacheItem **bucket;  /* The bucket "list heads" themselves */

    /* Cyclic cache array */
    HacheOrder *ordering; 
    int head, tail, free;

    /* Cache I/O callback functions */
    void      *clientdata;
    HacheData *(*load)(void *clientdata, char *key, int key_len, HacheItem *hi);
    void      (*del)(void *clientdata, HacheData data);

    int	      searches; /* number of total queries */
    int       hits;	/* number of cached queries */
} HacheTable;

/* An iterator on HacheTable items */
typedef struct {
    int bnum;
    HacheItem *hi;
} HacheIter;

/* Functions to to use HacheTable.options */
#define HASH_FUNC_HSIEH       0
#define HASH_FUNC_TCL         1
#define HASH_FUNC_JENKINS     2
#define HASH_FUNC_MASK        7

/* Other HacheTable.options values */
#define HASH_NONVOLATILE_KEYS (1<<3)
#define HASH_ALLOW_DUP_KEYS   (1<<4)
#define HASH_DYNAMIC_SIZE     (1<<5)
#define HASH_OWN_KEYS	      (1<<6)

/* Hacheing prototypes */
uint32_t hash(int func, uint8_t *key, int key_len);
uint32_t HacheJenkins(uint8_t *k, int length);
uint32_t HacheTcl(uint8_t *data, int len);
uint32_t HacheHsieh(uint8_t *k, int length);

/* HacheTable management prototypes */
HacheTable *HacheTableCreate(int size, int options);
void HacheTableDestroy(HacheTable *h, int deallocate_date);
int HacheTableResize(HacheTable *h, int newsize);
HacheItem *HacheTableAdd(HacheTable *h, char *key, int key_len,
			 HacheData data, int *new);
int HacheTableDel(HacheTable *h, HacheItem *hi, int deallocate_data);
int HacheTableRemove(HacheTable *h, char *key, int key_len,
		     int deallocate_data);
HacheItem *HacheTableQuery(HacheTable *h, char *key, int key_len);
HacheItem *HacheTableSearch(HacheTable *h, char *key, int key_len);
HacheItem *HacheTableNext(HacheItem *hi, char *key, int key_len);
void HacheTableIncRef(HacheTable *h, HacheItem *hi);
void HacheTableDecRef(HacheTable *h, HacheItem *hi);

void HacheTablePurge(HacheTable *h);
void HacheTableStats(HacheTable *h, FILE *fp);
void HacheTableDump(HacheTable *h, FILE *fp);

void HacheOrderAccess(HacheTable *h, HacheItem *hi);
void HacheOrderRemove(HacheTable *h, HacheItem *hi);
int HacheOrderAdd(HacheTable *h, HacheItem *hi);

/* Iterator prototypes */
HacheIter *HacheTableIterCreate(void);
void HacheTableIterDestroy(HacheIter *iter);
HacheItem *HacheTableIterNext(HacheTable *h, HacheIter *iter);
void HacheTableIterReset(HacheIter *iter);

#endif /* _HACHE_TABLE_H_ */
