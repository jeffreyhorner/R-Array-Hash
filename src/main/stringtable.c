#define USE_RINTERNALS

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#define R_USE_SIGNALS 1
#include <Defn.h>
#include <Internal.h>
#include <Print.h>

#define intCHARSXP 73 /* Very important that this is equivalent in memory.c */

#define MAX_CHAR_LEN 128

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
static unsigned int char_hash(const char *s, int len)
{
    /* djb2 as from http://www.cse.yorku.ca/~oz/hash.html */
    char *p;
    int i;
    unsigned int h = 5381;
    for (p = (char *) s, i = 0; i < len; p++, i++)
	h = ((h << 5) + h) + (*p);
    return h;
}
static inline R_len_t ST_SLOT(const char *s,R_len_t len){
     return (R_len_t)(char_hash(s,len) & (R_StringTable->len - 1));
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

void attribute_hidden InitStringTable()
{
    int memuse=0;
    R_StringTable = 
	(R_str_table_t *)calloc(1,sizeof(*R_StringTable));
    memuse += sizeof(R_str_table_t);
    if (!R_StringTable)
    	R_Suicide("couldn't allocate memory for string table");
    R_StringTable->slot = (R_str_slot_t **)calloc(STSIZE,sizeof(R_str_slot_t *));
    if (!R_StringTable->slot)
    	R_Suicide("couldn't allocate memory for string table");
    memuse += STSIZE * sizeof(R_str_slot_t *);
    R_StringTable->len = STSIZE;
    
    R_SymbolTable = NULL;

    /* Note: we don't want NA_STRING to be in the CHARSXP cache, so that
       mkChar("NA") is distinct from NA_STRING */
    /* NA_STRING */
    NA_STRING = allocVector(intCHARSXP,strlen("NA"));
    strcpy((char *)DATAPTR(NA_STRING),"NA");
    SETLENGTH(NA_STRING, 2);
    SET_CACHED(NA_STRING);
    R_print.na_string = NA_STRING;

    InformGCofMemUsage(memuse,TRUE);
}

static inline size_t st_multiple(size_t requested_size){
    size_t multiple = 8;
    return (
	( 
	    requested_size/multiple + 
	    (((requested_size % multiple)>0)? 1: 0)
	) * multiple
    );
}
static void *st_malloc(size_t size){
    void *memptr;
    if (posix_memalign(&memptr, 8, size)==0){
	InformGCofMemUsage(size,TRUE);
	return memptr;
    }
    return NULL;
}

static void *st_realloc(void *oldmem, size_t oldsize, size_t newsize){
    void *newmem;

    if (posix_memalign(&newmem, 8, newsize)!=0) return NULL;

    memcpy(newmem, oldmem, oldsize);

    InformGCofMemUsage(newsize-oldsize,TRUE);

    free(oldmem);

    return newmem;
}

static void fix_char_encoding(SEXP c,const char *name, int len){
    Rboolean embedNul = FALSE, is_ascii = TRUE;
    for (int slen = 0; slen < len; slen++) {
	if ((unsigned int) name[slen] > 127) is_ascii = FALSE;
	if (!name[slen]) embedNul = TRUE;
    }
    if (is_ascii) SET_ASCII(c);
    if (embedNul) {
	error(_("embedded nul in string: '%s'"),
	      EncodeString(c, 0, 0, Rprt_adj_none));
    }
}
static R_str_elem_t *SlotInsElem(R_len_t slotn, const char *name, R_len_t len)
{
    R_str_slot_t *slot;
    R_str_elem_t *e;
    size_t esize, size;
    SEXP c;

    /* Make non-volatile string if it can fit in two VECREC's */
    if (BYTE2VEC(len+1) <= 2){
	PROTECT(c = allocVector(intCHARSXP,len));
    } else {
	PROTECT(c = allocVector(intCHARSXP,1));
	SET_STRING_VOLATILE(c);
    }
    SETLENGTH(c, len);
    /* HASHVALUE is slot */
    SET_HASHVALUE(c,slotn);
    SET_HASHASH(c,1);

    slot = R_StringTable->slot[slotn];
    /* TODO: determine best word alignment: 8, 16, ? */
    esize = st_multiple(sizeof(R_str_elem_t) + BYTE2VEC(len + 1)*sizeof(VECREC));

    if (slot){
	R_str_elem_t *elem, *last_elem;
	size_t oldsize = slot->size;
	size = st_multiple(oldsize + esize);
	slot = (R_str_slot_t *)st_realloc(slot,oldsize, size);
	if (!slot){
	    UNPROTECT(1);
	    error("couldn't allocate memory for string table");
	}
	slot->size = size;
	slot->dsize = 0;
	last_elem = (R_str_elem_t *)((char *)slot + oldsize);

	/* Fix up charsxp pointers to str for all but last element*/
	elem = FIRST_ELEMENT(slot);
	while (elem != last_elem){
	    if (elem->charsxp){
		if (IS_STRING_VOLATILE(elem->charsxp))
		    SET_CHAR_PTR(elem->charsxp,ELEMENT_CHAR(elem));
	    }
	    elem = NEXT_ELEMENT(elem);
	}
	e = last_elem;
    } else {
	size = st_multiple(sizeof(R_str_slot_t) + esize);
	esize = size - sizeof(R_str_slot_t);
	slot = (R_str_slot_t *)st_malloc(size);
	if (!slot){
	    UNPROTECT(1);
	    error("couldn't allocate memory for string table");
	}
	slot->size = size;
	slot->dsize = 0;
	e = FIRST_ELEMENT(slot);
    }
    R_StringTable->slot[slotn] = slot;

    e->charsxp = c;
    e->symsxp = NULL;
    e->size = esize;
    e->len = len;
    e->encoding = 0;
    if (IS_STRING_VOLATILE(c)){
	SET_CHAR_PTR(c, ELEMENT_CHAR(e));
    } else {
	memcpy((char *)DATAPTR(c),name,len);
	((char *)DATAPTR(c))[len] = 0;
    }
    memcpy(ELEMENT_CHAR(e),name,len);
    ELEMENT_CHAR(e)[len] = 0;

    UNPROTECT(1);

    return e;
}

/* Insert CHARSXP using String. Called from allocCharsxp in memory.c */
SEXP R_STInsChrStr(const char *name, R_len_t len)
{
    R_str_elem_t *e = SlotInsElem(ST_SLOT(name,len), name, len);
    fix_char_encoding(e->charsxp,name, len);
    
    return e->charsxp;
}

void R_STCompactSlot(R_str_slot_t *slot, int i)
{
    InformGCofMemUsage(slot->size,FALSE);

    if (slot->dsize == slot->size){
	free(slot);
	R_StringTable->slot[i] = NULL;
	return;
    }
    if (slot->dsize <=0) return;

    /* Would be nice to defer this until really necessary */
    size_t size = st_multiple(slot->size - slot->dsize);
    R_str_slot_t *newslot = (R_str_slot_t *)st_malloc(size);
    newslot->size = size;
    newslot->dsize = 0;
    R_str_elem_t *newelem = FIRST_ELEMENT(newslot);
    R_str_elem_t *elem = FIRST_ELEMENT(slot);
    R_str_elem_t *last_elem = (R_str_elem_t *)((char *)slot + slot->size);
    while (elem != last_elem){
	if (elem->charsxp){
	    memcpy(newelem,elem,elem->size);
	    if (IS_STRING_VOLATILE(newelem->charsxp))
		SET_CHAR_PTR(newelem->charsxp,ELEMENT_CHAR(newelem));
	    newelem = NEXT_ELEMENT(newelem);
	}
	elem = NEXT_ELEMENT(elem);
    }
    free(slot);
    R_StringTable->slot[i] = newslot;
}

inline const char *R_STCHAR(SEXP charsxp){
    /*if (TYPEOF(charsxp) != CHARSXP){
	error("Not a CHARSXP but a %d!\n",TYPEOF(charsxp));
    }*/

    if (!IS_STRING_VOLATILE(charsxp)) return ((const char *)DATAPTR(charsxp));

    const char *str = *(char **)DATAPTR(charsxp);
    /*R_str_elem_t *e;
    
    e = (R_str_elem_t *)str - 1;
    if (e->charsxp != charsxp){
	const char *str2 = *(char **)DATAPTR(e->charsxp);
	if (IS_STRING_VOLATILE(e->charsxp))
	    str2 = *(char **)DATAPTR(e->charsxp);
	else
	    str2 = ((const char *)DATAPTR(e->charsxp));
	error("CHARSXP MISMATCH! '%p:%s' '%p:%s'", e->charsxp, str, charsxp, str2);
    }*/
    
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
	printf("%lu, ", (long unsigned int)sum);
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

#include <string.h>

static char *xlate_encoding(SEXP c){
    if (IS_UTF8(c)) return "UTF8";
    if (IS_ASCII(c)) return "ASCII";
    if (IS_BYTES(c)) return "BYTES";
    return "NONE";
}
void do_write_cache()
{
    FILE *f = fopen("/tmp/CACHE", "w");
    char *str;
    char *c;
    int j=1;
    
    if (f != NULL) {
	fprintf(f,"\"slot\" \"sym\" \"encoding\" \"word\"\n");
	for (R_len_t i=0; i < R_StringTable->len; i++){
	    TRAVERSE_ST_SLOT(i)
		str = malloc(__e__->len+1);
		strcpy(str,ELEMENT_CHAR(__e__));
		while((c = strchr(str,'\n')) != NULL) {
		    *(c++) = '|';
		}
		while((c = strchr(str,'"')) != NULL) {
		    *(c++) = '\'';
		}
		fprintf(f,"\"%d\" %d %s \"%s\" \"%s\"\n", 
		    j++, i, (__e__->symsxp)? "TRUE": "FALSE", xlate_encoding(__e__->charsxp),str);
		free(str);
	    END_TRAVERSE_ST_SLOT
	}
	fclose(f);
    }
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

    /*if (len > MAX_CHAR_LEN){
	SEXP c;
	PROTECT(c = mkCharLen(name,len));
	sym = mkSYMSXP(c, R_UnboundValue);
	INSERT_SYMBOL(sym);
	UNPROTECT(1);
	return sym;
    }*/

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

    e = SlotInsElem(slotn, name, len);
    fix_char_encoding(e->charsxp,name, len);
    PROTECT(e->charsxp);
    PROTECT(sym = mkSYMSXP(e->charsxp,R_UnboundValue));
    INSERT_SYMBOL(sym);
    e->symsxp = sym;

    UNPROTECT(2);

    return sym;
}

SEXP installChar(SEXP charSXP)
{
    SEXP sym, csxp;
    const char *str;

    /* Really? */
    if (charSXP == NA_STRING){
	return mkSYMSXP(NA_STRING, R_UnboundValue);
    }

    R_len_t len = LENGTH(charSXP);

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
    }

    if (IS_STRING_VOLATILE(csxp)){
	str = *(char **)DATAPTR(csxp);
	R_str_elem_t *e = (R_str_elem_t *)str - 1;

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
	error(_("volatile CHARSXP not cached! '%s'"),CHAR(charSXP));
    } else if (HASHASH(charSXP)){
	TRAVERSE_ST_SLOT(HASHVALUE(csxp)){
	    if (__e__->charsxp == csxp){
		/* Symbol already exists */
		if (__e__->symsxp) return __e__->symsxp;

		sym = mkSYMSXP(__e__->charsxp, R_UnboundValue);
		INSERT_SYMBOL(sym);
		/* Add to cache */
		__e__->symsxp = sym;
		return sym;
	    }
	} END_TRAVERSE_ST_SLOT;
	error(_("non-volatile CHARSXP not cached! '%s'"),CHAR(charSXP));
    } else {
	error(_("non-cached CHARSXP in the wild! '%s'"),CHAR(charSXP));
	/*PROTECT(csxp);
	sym = mkSYMSXP(csxp, R_UnboundValue);
	INSERT_SYMBOL(sym);
	UNPROTECT(1);
	return sym;*/
    }
}

/* mkCharCE - make a character (CHARSXP) variable and set its
   encoding bit.  If a CHARSXP with the same string already exists in
   the global CHARSXP cache, R_StringHash, it is returned.  Otherwise,
   a new CHARSXP is created, added to the cache and then returned. */

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


/*
SEXP mkCharLenCE_NOCACHE(const char *name, int len, cetype_t enc)
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
	SEXP c = allocVector(intCHARSXP,len);
	memcpy((char *)DATAPTR(c),name,len);
	((char *)DATAPTR(c))[len] = 0;
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

    cval = allocVector(intCHARSXP,len);
    memcpy((char *)DATAPTR(cval),name,len);
    ((char *)DATAPTR(cval))[len] = 0;

    switch(enc) {
    case CE_NATIVE:
	break;
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
    SET_CACHED(cval);

    return cval;
}*/

SEXP mkCharLenCE(const char *name, int len, cetype_t enc)
{
    SEXP cval;
    int need_enc;
    Rboolean embedNul = FALSE, is_ascii = TRUE;
    R_len_t slotn;
    R_str_elem_t *e;

    /*if (len > MAX_CHAR_LEN) return mkCharLenCE_NOCACHE(name, len, enc);*/

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
	/* This is tricky: we want to make a reasonable job of
	   representing this string, and EncodeString() is the most
	   comprehensive */
	e = SlotInsElem(ST_SLOT(name,len), name, len);
	SEXP c = e->charsxp;
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
    cval = traverse_st_slot(slotn, name, len, need_enc);
    if (cval != R_NilValue) return cval;
    
    /* no cached value; need to allocate one and add to the cache */
    e = SlotInsElem(slotn, name, len);
    cval = e->charsxp;
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

    return cval;
}
