/*
cache.c
caching of Cget and Dget results in dynamically allocated memory
this file is part of LoopTools
last modified 7 Jun 01 th
*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#ifdef HAVE_UNDERSCORE
#define cachelookup cachelookup_
#define setcachelast setcachelast_
#define getcachelast getcachelast_
#endif


#ifndef REALSIZE
#define REALSIZE sizeof(double)
#endif

#define CPLXSIZE (2*REALSIZE)


int cachelookup(char *para, char **base,
  void (*calc)(char *, char *), int *npara, int *nval)
{
  char *node, *last, **next = base;
  int parasize = *npara*REALSIZE, valsize = *nval*CPLXSIZE;

  if(last = *(base + 1)) {
    do {
      node = *next;
      next = (char **)(node - sizeof(char *));
      if(memcmp(para, node + valsize, parasize) == 0) goto done;
    } while(node != last);
  }
  node = *next;
  if(node == NULL) {
    node = malloc(sizeof(char *) + valsize + parasize + CPLXSIZE);
	/* MUST have extra CPLXSIZE for alignment so that node
	   can be reached with an integer offset into base */
    if(node == NULL) {
      fputs("Out of cache memory.\n", stderr);
      exit(1);
    }
    node += sizeof(char *);
    node += ((long)base - (long)node) & (CPLXSIZE - 1);
    *next = node;
    *(char **)(node - sizeof(char *)) = NULL;
  }
  *(base + 1) = node;
  memcpy(node + valsize, para, parasize);
  calc(para, node);
done:
  return ((long)node - (long)base)/(int)CPLXSIZE;
}


void setcachelast(char **base, int *offset)
{
  *(base + 1) =
    *offset ? (char *)((long)base + *offset*CPLXSIZE) : NULL;
}


int getcachelast(char **base)
{
  return *(base + 1) ?
    ((long)*(base + 1) - (long)base)/(int)CPLXSIZE : 0;
}

