#define USE_RINTERNALS

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#define R_USE_SIGNALS 1
#include <Defn.h>
#include <Internal.h>
#include <Print.h>

#define intCHARSXP 73 /* Very important that this is equivalent in memory.c */

/* Use __e__ to acess the current element */
#define TRAVERSE_ST_SLOT(x) {\
	size_t __esize__; \
	R_len_t __slotn__ = (x); \
	R_str_slot_t *__slot__ = R_StringTable->slot[__slotn__]; \
	if (__slot__) { \
	R_str_elem_t *__e__=(R_str_elem_t *)(__slot__+1); \
	R_str_elem_t *__end__=(R_str_elem_t *)((char *)__slot__ + __slot__->size); \
	while (__e__ != __end__) { \
	    __esize__ = __e__->size;

#define END_TRAVERSE_ST_SLOT __e__ = (R_str_elem_t *)((char *)__e__ + __esize__); } } }

#define ELEMENT_CHAR(e) ((char *)((R_str_elem_t *)e+1))
#define FIRST_ELEMENT(s) ((R_str_elem_t *)((R_str_slot_t *)s+1))
#define NEXT_ELEMENT(e) ((R_str_elem_t *)((char *)e+e->size))
#define SLOT_END(s) ((R_str_elem_t *)((char *)s+s->size))
static inline R_len_t ST_SLOT(const char *s,R_len_t len){
     return (R_len_t)(StringHash(s,len) & (R_StringTable->len - 1));
}

#define INSERT_SYMBOL(s) do { \
    R_sym_table_t *__sym__ = calloc(1,sizeof(R_sym_table_t)); \
    if (!__sym__) \
	error("couldn't allocate memory for symbol table"); \
    __sym__->symsxp = s; \
    if (R_SymbolTable) __sym__->next = R_SymbolTable; \
    R_SymbolTable = __sym__; \
} while(0)


#define STSIZE 65536

#define SET_CHAR_PTR(c,s) *((char **)DATAPTR(c)) = s

/* Encodings on R_str_elem_t ints */
#define ST_SET_UTF8(x) ((x)->encoding |= UTF8_MASK)
#define ST_SET_LATIN1(x) ((x)->encoding |= LATIN1_MASK)
#define ST_SET_BYTES(x) ((x)->encoding |= BYTES_MASK)
#define ST_ENC_KNOWN(x) ((x)->encoding & (LATIN1_MASK | UTF8_MASK))
#define ST_IS_BYTES(x) ((x)->encoding & BYTES_MASK)

static R_str_elem_t *ST_NA_STRING;

void attribute_hidden InitStringTable()
{
    R_StringTable = 
	(R_str_table_t *)calloc(1,sizeof(*R_StringTable));
    if (!R_StringTable)
    	R_Suicide("couldn't allocate memory for string table");
    R_StringTable->slot = (R_str_slot_t **)calloc(STSIZE,sizeof(R_str_slot_t *));
    if (!R_StringTable->slot)
    	R_Suicide("couldn't allocate memory for string table");
    R_StringTable->len = STSIZE;
    
    R_SymbolTable = NULL;

    /* Note: we don't want NA_STRING to be in the CHARSXP cache, so that
       mkChar("NA") is distinct from NA_STRING */
    /* NA_STRING */
    NA_STRING = allocVector(intCHARSXP,sizeof(char *));
    ST_NA_STRING = calloc(1,sizeof(R_str_elem_t) + 3);
    ST_NA_STRING->symsxp = NULL;
    ST_NA_STRING->charsxp = NA_STRING;
    ST_NA_STRING->len = 2;
    memcpy(ELEMENT_CHAR(ST_NA_STRING),"NA",2);
    ELEMENT_CHAR(ST_NA_STRING)[2] = 0;
    SET_CHAR_PTR(NA_STRING, ELEMENT_CHAR(ST_NA_STRING));
    SETLENGTH(NA_STRING, 2);
    SET_CACHED(NA_STRING);
    R_print.na_string = NA_STRING;
}

#include <stdlib.h>
static size_t st_multiple(size_t requested_size, size_t multiple){
    return (
	( 
	    requested_size/multiple + 
	    (((requested_size % multiple)>0)? 1: 0)
	) * multiple
    );
}
static void *st_malloc(size_t size){
    void *memptr;
    if (posix_memalign(&memptr, 8, size)==0) return memptr;
    return NULL;
}

static void *st_realloc(void *oldmem, size_t oldsize, size_t newsize){
    void *newmem;

    if (posix_memalign(&newmem, 8, newsize)!=0) return NULL;

    memcpy(newmem, oldmem, oldsize);

    free(oldmem);

    return newmem;
}

static R_str_elem_t *SlotInsElem(R_len_t slotn, SEXP charSXP, const char *name, R_len_t len)
{
    R_str_slot_t *slot;
    R_str_elem_t *e;
    size_t esize, size;

    slot = R_StringTable->slot[slotn];
    /* TODO: determine best word alignment: 8, 16, ? */
    esize = st_multiple(sizeof(R_str_elem_t) + BYTE2VEC(len + 1)*sizeof(VECREC), 8);

    if (slot){
	R_str_elem_t *elem, *last_elem;
	size_t oldsize = slot->size;
	size = st_multiple(oldsize + esize, 8);
	slot = (R_str_slot_t *)st_realloc(slot,oldsize, size);
	if (!slot)
	    error("couldn't allocate memory for string table");
	slot->size = size;
	last_elem = (R_str_elem_t *)((char *)slot + oldsize);

	/* Fix up charsxp pointers to str for all but last element*/
	elem = FIRST_ELEMENT(slot);
	while (elem != last_elem){
	    if (elem->charsxp){
		SET_CHAR_PTR(elem->charsxp,ELEMENT_CHAR(elem));
	    }
	    elem = NEXT_ELEMENT(elem);
	    
	}
	e = last_elem;
    } else {
	size = st_multiple(sizeof(R_str_slot_t) + esize, 8);
	esize = size - sizeof(R_str_slot_t);
	slot = (R_str_slot_t *)st_malloc(size);
	if (!slot)
	    error("couldn't allocate memory for string table");
	slot->size = size;
	e = FIRST_ELEMENT(slot);
    }
    R_StringTable->slot[slotn] = slot;

    e->symsxp = NULL;
    e->charsxp = charSXP;
    e->size = esize;
    e->len = len;
    e->encoding = 0;
    memcpy(ELEMENT_CHAR(e),name,len);
    ELEMENT_CHAR(e)[len] = 0;
    SET_CHAR_PTR(charSXP, ELEMENT_CHAR(e));
    SETLENGTH(charSXP, len);

    return e;
}

/* Insert CHARSXP using String. Called from allocCharsxp in memory.c */
SEXP R_STInsChrStr(SEXP charSXP, const char *name, R_len_t len)
{
    SlotInsElem(ST_SLOT(name,len),charSXP, name, len);
    return charSXP;
}

void R_STDelCHAR(const char *str){
    R_str_elem_t *e;
    
    e = (R_str_elem_t *)str - 1;
    if (!e->charsxp){
	error("CHARSXP(%s) IS NULL!\n", str);
    }
    if (e->symsxp){
	error("Cannot Delete CHARSXP(%s) because of SYMSXP!\n",str);
    }
    e->charsxp = NULL;
}

const char *R_STCHAR(SEXP charsxp){
    if (TYPEOF(charsxp) != CHARSXP){
	error("Not a CHARSXP but a %d!\n",TYPEOF(charsxp));
    }
    const char *str = *(char **)DATAPTR(charsxp);
    R_str_elem_t *e;
    
    e = (R_str_elem_t *)str - 1;
    if (e->charsxp != charsxp){
	const char *str2 = *(char **)DATAPTR(e->charsxp);
	error("CHARSXP MISMATCH! '%s' '%s'", str, str2);
    }
    
    return str;
}

/* Call this from gdb with

       call do_show_cache(10)

   for the first 10 cache chains in use. */
void do_show_cache()
{
    for (R_len_t i=0; i < R_StringTable->len; i++){
	size_t sum=0;
	TRAVERSE_ST_SLOT(i)
	    sum++;
	END_TRAVERSE_ST_SLOT
	printf("%lu, ", sum);
    }
    printf("\n");
}
//    int i, j;
//    Rprintf("Cache size: %d\n", LENGTH(R_StringHash));
//    Rprintf("Cache pri:  %d\n", HASHPRI(R_StringHash));
//    for (i = 0, j = 0; j < n && i < LENGTH(R_StringHash); i++) {
//	SEXP chain = VECTOR_ELT(R_StringHash, i);
//	if (! ISNULL(chain)) {
//	    Rprintf("Line %d: ", i);
//	    do {
//		if (IS_UTF8(CXHEAD(chain)))
//		    Rprintf("U");
//		else if (IS_LATIN1(CXHEAD(chain)))
//		    Rprintf("L");
//		else if (IS_BYTES(CXHEAD(chain)))
//		    Rprintf("B");
//		Rprintf("|%s| ", CHAR(CXHEAD(chain)));
//		chain = CXTAIL(chain);
//	    } while(! ISNULL(chain));
//	    Rprintf("\n");
//	    j++;
//	}
//    }
//}

void do_write_cache()
{
//    int i;
//    FILE *f = fopen("/tmp/CACHE", "w");
//    if (f != NULL) {
//	fprintf(f, "Cache size: %d\n", LENGTH(R_StringHash));
//	fprintf(f, "Cache pri:  %d\n", HASHPRI(R_StringHash));
//	for (i = 0; i < LENGTH(R_StringHash); i++) {
//	    SEXP chain = VECTOR_ELT(R_StringHash, i);
//	    if (! ISNULL(chain)) {
//		fprintf(f, "Line %d: ", i);
//		do {
//		    if (IS_UTF8(CXHEAD(chain)))
//			fprintf(f, "U");
//		    else if (IS_LATIN1(CXHEAD(chain)))
//			fprintf(f, "L");
//		    else if (IS_BYTES(CXHEAD(chain)))
//			fprintf(f, "B");
//		    fprintf(f, "|%s| ", CHAR(CXHEAD(chain)));
//		    chain = CXTAIL(chain);
//		} while(! ISNULL(chain));
//		fprintf(f, "\n");
//	    }
//	}
//	fclose(f);
//    }
}


/*  install - probe the symbol table */
/*  If "name" is not found, it is installed in the symbol table.
    The symbol corresponding to the string "name" is returned. */

SEXP install(const char *name)
{
    SEXP sym;
    size_t len;
    R_str_elem_t *e;
    R_len_t slotn;
    const char *estr;

    /* Quick sanity check */
    if (*name == '\0')
	error(_("attempt to use zero-length variable name"));

    /* Variable name length check */
    len = strlen(name);
    if (len > MAXIDSIZE)
	error(_("variable names are limited to %d bytes"), MAXIDSIZE);

    slotn = ST_SLOT(name,len);
    TRAVERSE_ST_SLOT(slotn)
    {
	if (__e__->charsxp && __e__->len == len){
	    estr = ELEMENT_CHAR(__e__);
	    sym = __e__->symsxp;
	    if ((estr == name) || memcmp(ELEMENT_CHAR(__e__),name,len)==0){
		/* Symbol already exists */
		if (sym) return sym;

		/* Charsxp exists, need to create sym. */
		sym = mkSYMSXP(__e__->charsxp, R_UnboundValue);
		INSERT_SYMBOL(sym);

		/* Add to cache */
		__e__->symsxp = sym;

		return sym;
	    }
	}
    } END_TRAVERSE_ST_SLOT;

    e = SlotInsElem(slotn, allocVector(intCHARSXP,sizeof(char *)),name, len);
    sym = mkSYMSXP(e->charsxp,R_UnboundValue);
    INSERT_SYMBOL(sym);
    e->symsxp = sym;

    return sym;
}

SEXP installChar(SEXP charSXP)
{
    SEXP sym, csxp;
    const char *str = *(char **)DATAPTR(charSXP);
    R_str_elem_t *e = (R_str_elem_t *)str - 1;
    R_len_t len = e->len;

    if (len == 0)
        error(_("attempt to use zero-length variable name"));
    if (len > MAXIDSIZE)
        error(_("variable names are limited to %d bytes"), MAXIDSIZE);

    if (IS_ASCII(charSXP) || (IS_UTF8(charSXP) && utf8locale) ||
                                        (IS_LATIN1(charSXP) && latin1locale) ){
	csxp = charSXP;
    } else {
        /* This branch is to match behaviour of install (which is older):
           symbol C-string names are always interpreted as if
           in the native locale, even when they are not in the native locale */
	csxp = mkCharLen(CHAR(charSXP), len);
	str = *(char **)DATAPTR(csxp);
	e = (R_str_elem_t *)str - 1;
    }

    if (e->charsxp == csxp){
	/* Symbol already exists */
	if (e->symsxp) return e->symsxp;

	sym = mkSYMSXP(e->charsxp, R_UnboundValue);
	INSERT_SYMBOL(sym);
	/* Add to cache */
	e->symsxp = sym;
	return sym;
    }
    /* Should never happen */
    error(_("non-cached CHARSXP in the wild!"));
}

SEXP mkCharCE(const char *name, cetype_t enc)
{
    size_t len =  strlen(name);
    if (len > INT_MAX)
	error("R character strings are limited to 2^31-1 bytes");
   return mkCharLenCE(name, (int) len, enc);
}

/* no longer used in R but docuented in 2.7.x */
SEXP mkCharLen(const char *name, int len)
{
    return mkCharLenCE(name, len, CE_NATIVE);
}

SEXP mkChar(const char *name)
{
    size_t len =  strlen(name);
    if (len > INT_MAX)
	error("R character strings are limited to 2^31-1 bytes");
    return mkCharLenCE(name, (int) len, CE_NATIVE);
}

static inline SEXP traverse_st_slot(R_len_t slotn, const char *name, int len, int need_enc){
    const char *estr;
    size_t esize;
    R_str_slot_t *__slot__ = R_StringTable->slot[slotn];
    if (!__slot__) return R_NilValue;
    R_str_elem_t *__e__=(R_str_elem_t *)(__slot__+1);
    R_str_elem_t *__end__=(R_str_elem_t *)((char *)__slot__ + __slot__->size);
    while (__e__ != __end__) {
	esize = __e__->size;
	if (
	    __e__->charsxp && 
	    __e__->len == len &&
	    (need_enc == (ST_ENC_KNOWN(__e__) | ST_IS_BYTES(__e__)))
	    ) {
	    estr = ELEMENT_CHAR(__e__);
	    if ( (estr == name) || (memcmp(estr, name, len) == 0) )
		return __e__->charsxp;
	}

	__e__ = (R_str_elem_t *)((char *)__e__ + esize);
    }

    return R_NilValue;
}

/* mkCharCE - make a character (CHARSXP) variable and set its
   encoding bit.  If a CHARSXP with the same string already exists in
   the global CHARSXP cache, R_StringHash, it is returned.  Otherwise,
   a new CHARSXP is created, added to the cache and then returned. */


SEXP mkCharLenCE(const char *name, int len, cetype_t enc)
{
    SEXP cval;
    int need_enc;
    Rboolean embedNul = FALSE, is_ascii = TRUE;
    R_len_t slotn;
    R_str_elem_t *e;

    switch(enc){
    case CE_NATIVE:
    case CE_UTF8:
    case CE_LATIN1:
    case CE_BYTES:
    case CE_SYMBOL:
    case CE_ANY:
	break;
    default:
	error(_("unknown encoding: %d"), enc);
    }
    for (int slen = 0; slen < len; slen++) {
	if ((unsigned int) name[slen] > 127) is_ascii = FALSE;
	if (!name[slen]) embedNul = TRUE;
    }
    if (embedNul) {
	SEXP c;
	/* This is tricky: we want to make a reasonable job of
	   representing this string, and EncodeString() is the most
	   comprehensive */
	c = allocVector(intCHARSXP, sizeof(char *));
	SlotInsElem(ST_SLOT(name,len),c, name, len);
	switch(enc) {
	case CE_UTF8: SET_UTF8(c); break;
	case CE_LATIN1: SET_LATIN1(c); break;
	case CE_BYTES: SET_BYTES(c); break;
	default: break;
	}
	if (is_ascii) SET_ASCII(c);
	error(_("embedded nul in string: '%s'"),
	      EncodeString(c, 0, 0, Rprt_adj_none));
    }

    if (enc && is_ascii) enc = CE_NATIVE;
    switch(enc) {
    case CE_UTF8: need_enc = UTF8_MASK; break;
    case CE_LATIN1: need_enc = LATIN1_MASK; break;
    case CE_BYTES: need_enc = BYTES_MASK; break;
    default: need_enc = 0;
    }

    slotn = ST_SLOT(name,len);
    /*cval = traverse_st_slot_stats(slotn, name, len, need_enc);*/
    cval = traverse_st_slot(slotn, name, len, need_enc);
    if (cval != R_NilValue) return cval;
    
    /* no cached value; need to allocate one and add to the cache */
    PROTECT(cval = allocVector(intCHARSXP, sizeof(char *)));
    e = SlotInsElem(slotn,cval, name, len);
    switch(enc) {
    case CE_NATIVE:
	break;          /* don't set encoding */
    case CE_UTF8:
	SET_UTF8(cval); ST_SET_UTF8(e);
	break;
    case CE_LATIN1:
	SET_LATIN1(cval); ST_SET_LATIN1(e);
	break;
    case CE_BYTES:
	SET_BYTES(cval); ST_SET_BYTES(e);
	break;
    default:
	error("unknown encoding mask: %d", enc);
    }
    if (is_ascii) SET_ASCII(cval);
    SET_CACHED(cval);  /* Mark it */

    UNPROTECT(1);
    
    return cval;
}
