#include "dwtnode.h"
struct jprimevec{float jprime; bool prune;};
struct bestJvec{float J; direction cdir; char cshift;};
class orienttree
{
public:
  orienttree(bool setdir, char setshift, bool setleaf)
  { 
    dir=setdir;
    shift=setshift;
    leaf=setleaf;
    jprimecand=NULL;
    Jvec=NULL;
    for (int i=0;i<4;i++)
      children[i]=NULL;
  }
  ~orienttree()
  {
    if (jprimecand!=NULL) delete[] jprimecand;
    if (Jvec!=NULL) delete[] Jvec;
    for (int i=0;i<4;i++)
      if (children[i]!=NULL) delete children[i];
  }
  bool dir; // 0==false==vertical, 1==true==horizontal
  char shift;
  bool leaf;
  jprimevec *jprimecand; // array storing the best J' for each possible shift for this node
  // array storing the best J and associated shift for this node for each possible parent shift
  bestJvec *Jvec;
  orienttree *children[4]; // children[1] is offset in the x direction, child[2] in y
};
void codetree(char *fname, orienttree *root);