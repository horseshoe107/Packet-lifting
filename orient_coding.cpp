#include "stdafx.h"
#include "dwtnode.h"
#include "orient.h"
class bitstream_output
{
public:
  ~bitstream_output();
  void open(char *filename);
  void push_data(bool *, int size);
private:
  int position;
  unsigned char buffer;
  std::ofstream bitsout;
} bout;
int bitscount;
bitstream_output::~bitstream_output()
{
  if (position != 0)
    bitsout << buffer;
  bitsout.close();
}
void bitstream_output::open(char *fname)
{
  position=0;
  buffer=0;
  bitsout.open(fname,std::ios::binary);
}
void bitstream_output::push_data(bool *bits,int size)
{
  if ((size<=0)||(size>20))
  {
    std::cerr << size << " number of bits cannot be written" << std::endl;
    exit(1);
  }
  for (int i=0;i<size;i++)
  {
    buffer |= (bits[i] << (position%8)); // note the bit must be initialised 0 first
    position++;
    if (position%8 == 0)
    { // write out a char
      bitsout << buffer;
      position = 0;
      buffer=0; // clear all bits to 0
    }
  }
  bitscount+=size;
  return;
}
// begin by simply encoding each element of the orientation field independently
void orientationfield::orientencode(char *fname)
{
  bout.open(fname);
  bool bits[20];
  int size, shift;
  for (int i=0;i<numblks;i++)
  {
    if ((orientvec[i].hshift!=0)&&(orientvec[i].vshift!=0))
    {
      std::cerr << "Error in orientation field: non zero shifts for both transform directions" << std::endl;
      exit(2);
    }
    else
    {
      if ((orientvec[i].hshift==0)&&(orientvec[i].vshift==0))
      {
        size = 1;
        bits[0] = 0;
      }
      else
      {
        bits[0] = 1;
        bits[1] = (orientvec[i].vshift==0);
        shift = (orientvec[i].vshift==0) ? orientvec[i].hshift:orientvec[i].vshift;
        bits[2] = (shift<0); // sign bit
        shift = abs(shift);
        // encode shift with golomb coding/exponential golomb
        // comma code implementation
        int j=1;
        for (;j<shift;j++)
          bits[2+j]=0;
        bits[2+j]=1;
        size = shift+3;
      }
    }
    bout.push_data(bits,size);
  }
}
// encodes the current node, calling itself recursively for children
void recursiveencode(orienttree *node, char parentdata, bool parentdir)
{
  static bool bits[20];
  int size,shift;
  if (parentdata==0) // encode independently
  {
    if (node->shift == 0)
    {
      bits[0]=0;
      size=1;
    }
    else
    {
      bits[0]=1;
      bits[1]=node->dir; // indicate direction
      bits[2]=(node->shift)<0;
      shift=abs(node->shift);
      // comma code implementation
      int j=1;
      for (;j<shift;j++)
        bits[2+j]=0;
      bits[2+j]=1;
      size = shift+3;
    }
  }
  else
  {
    if (parentdir == node->dir) //code differentially
    {
      bits[0]=0; // signal same direction
      shift = (node->shift) - parentdata;
    }
    else // code independently
    {
      bits[0]=1; // signal orthogonal direction
      shift = node->shift;
    }
    if (shift == 0)
    {
      bits[1]=0; // signal zero and end
      size=2;
    }
    else
    {
      bits[1]=1; // signal non-zero differential
      bits[2] = shift<0; // signal sign
      shift = abs(shift);
      // comma code implementation
      int j=1;
      for (;j<shift;j++)
        bits[2+j]=0;
      bits[2+j]=1;
      size = shift+3;
    }
  }
  // output a final bit to indicate whether this is a leaf node
  if (node->leaf)
  {
    bits[size]=1;
    bout.push_data(bits,size+1);
  }
  else
  {
    bits[size]=0;
    bout.push_data(bits,size+1);
    recursiveencode(node->children[0],node->shift,node->dir);
    recursiveencode(node->children[1],node->shift,node->dir);
    recursiveencode(node->children[2],node->shift,node->dir);
    recursiveencode(node->children[3],node->shift,node->dir);
  }
  return;
}
void codetree(char *fname, orienttree *root)
{
  bitscount=0;
  bout.open(fname);
  recursiveencode(root,0,true);
  std::cout << bitscount << std::endl;
}