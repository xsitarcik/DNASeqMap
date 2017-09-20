#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "compression.h"

struct symbol_table *build_symbol_table(char *alphabet)
{
 int alphabet_size = strlen(alphabet);
 unsigned int i = 0;
 struct symbol_table *front = NULL;
 struct symbol_table *prev = (struct symbol_table*)malloc(sizeof(struct symbol_table));
 struct symbol_table *next = NULL;
 prev->index = i;
 prev->symbol = alphabet[i];
 prev->previous = NULL;
 if (prev==NULL)
  {
   printf("Error when allocating memory for symbol table\n");
   exit(1);
  }
 front = prev; //for returning
 for (i=1;i<alphabet_size;i++){
  next = (struct symbol_table*)malloc(sizeof(struct symbol_table));
  if (next==NULL)
  {
   printf("Error when allocating memory for symbol table\n");
   exit(1);
  }
  next->index = i;
  next->symbol = alphabet[i];
  next->previous = prev;
  prev->next = next;
  prev = next;
 }
 next->next = NULL;
 return front;
}

unsigned int move_to_front_encode (struct symbol_table *front, char *s)
{
 struct symbol_table *current = front;
 while (current->next!=NULL)
 {
  printf("inddex %d = %c\n",current->index,current->symbol);
  current = current->next;
 }
 printf("inddex %d = %c\n",current->index,current->symbol);
 return 1;
}
