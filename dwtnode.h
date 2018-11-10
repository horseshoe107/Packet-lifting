enum direction {vertical,horizontal,both};
enum dwttype {w5x3, w9x7, disabled};
enum testmode {base,pyramid3x2,pyramid2x3,packlift,orient,aaorient,hpfprelift};
class orientationfield
{
  friend class dwtnode;
  friend class estorient;
  friend class est2orient;
public:
  // io functions (defined in support_io.cpp)
  orientationfield(){h=0, w=0, blksz=0, numblks=0, oprec=0,
    maxshift=0, fieldtype=blockgrid, orientvec=NULL;};
  void init_orient(int setblksz, int setprec, int setmaxshift,
    int defaulthorzshft=0, int defaultvertshft=0);
  void init_orient(char *fname);
  void copy(const orientationfield &);
  void inherit(orientationfield parent);
  ~orientationfield(){if (orientvec!=NULL) delete[] orientvec;}
  void orientwrite(char *fname="orientout.dat");
  // manipulation functions (defined in orient_tx.cpp)
  void clearfield(int seth, int setw);
  inline void setaffinefield(){fieldtype = affinegrid;}
  // returns the shift (in 1/oprec units) corresponding to
  // the oriented transform at pixel location (y,x).
  int retrieve(int y, int x, direction=vertical);
  inline int block_retrieve(int y, int x, direction dir)
  {
    if (y<0)        y=0;
    else if (y>=h)  y=h-1;
    if (x<0)        x=0;
    else if (x>=w)  x=w-1;
    int nblk = (y/blksz)*((w-1)/blksz+1) +x/blksz;
    return (dir==vertical)?orientvec[nblk].hshift:orientvec[nblk].vshift;
  }
  int affine_retrieve(int y, int x, direction dir);
  void transpose();
private:
  int h,w; // note the h,w dimensions must be the same as those
           // of the dwtnode container
  int blksz, numblks;
  int oprec, maxshift;
  enum {blockgrid, affinegrid} fieldtype;
  // This vector should have (N/blksz^2) elements, where N is the
  // number of elements in image, and blksz^2 is the size of an
  // "orientation block" - 1 pixel in the most precise case.
  // The value at orientvec[n].hshift refers to the relative shift
  // estimated between pairs of rows belonging to orientation block n.
  // For a block with blksz number of rows, the orientation applies to
  // rows pairs (0,1) through to (blksz-1,blksz). Note that (-1,0)
  // "belongs" to the previous orientation block
  struct orientation{char hshift; char vshift;} *orientvec;
};
class dwtnode
{
  friend class estorient;
public:
  // constructor functions (defined in dwtnode_io.cpp)
  dwtnode(int hset, int wset, dwttype, bool initzero=false);
  dwtnode(char *fname, dwttype);
  dwtnode(char *fname, int hset, int wset, dwttype, int expi=6);
	dwtnode(const dwtnode &copy);
  virtual ~dwtnode()
  {
    for (int n=0;n<4;n++)
      if (subbands[n] != NULL)
        delete subbands[n];
    if (pixels!=NULL)
      delete[] pixels;
  }
  // io functions (defined in dwtnode_io.cpp)
  void pgmread(char *fname);
  void rawlread(char *fname, int hset, int wset, int expi=6);
  void rawlread(char *fname, int expi=6);
  bool yuvstreamread(ifstream &yuvin);
  void pgmwrite(char *fname);
  void rawlwrite(char *fname, int expi=6, bool allbands=false);
  void csvwrite(char *fname);
  void yuvstreamwrite(ofstream &yuvout);
  // ordinary manipulation functions (defined in dwtnode_tx.cpp)
  void extract_subband(int band);
  void insert_subband(int band, bool suppress_warnings=false);
  void interleave(bool suppress_warnings=false);
  void transpose();
  double filt(double *f, int pixelloc, int offset, int fN,
    direction, bool forward=true);
  // ordinary transform functions (defined in dwtnode_tx.cpp)
  void apply_LHlift(double, direction);
  void apply_HLlift(double, direction);
  void apply_gain_factors(double k0, double k1, direction);
  void analysis(direction dir);
  void synthesis(direction dir);
  void packet_analysis(direction dir);
  void packet_synthesis(direction dir);
  // laplacian pyramid functions (defined in pyramid.cpp)
  void upsample_lift(bool analysis);
  void downsample_lift(bool analysis);
  void lp3x2_halfres();
  void lp3x2_encode(bool halfres, bool adapt);
  void lp3x2_decode(char *bitrate, bool adapt);
  void lp2x3_halfres();
  void lp2x3_encode(bool halfres, bool adapt);
  void lp2x3_decode(char *bitrate, bool adapt);
  // antialiasing transform (defined in antialias_tx.cpp)
  friend void packet_transfer(dwtnode &donor, dwtnode &receiver,
    bool analysis, direction);
  friend void packet_cancel(dwtnode &donor, dwtnode &receiver,
    bool analysis, direction);
  friend void packet_transfer_adaptive(dwtnode &donor, dwtnode &receiver,
    bool analysis, direction);
  friend void packswap(dwtnode &donor, dwtnode &receiver,
    bool analysis, bool adaptive, direction);
	friend void average3x3abs(dwtnode &in, dwtnode &out);
  friend double alphaTlookup(dwtnode &hE, dwtnode &donorE, int y, int x, direction);
  void packlift(direction dim, bool analysis, bool adaptive=false);
  // oriented transform functions (defined in orient_tx.cpp)
	void kernel_selection(int n, int sigma, direction, int &intshift,
                        int &LUT_index, bool &shiftdirection);
  void apply_oriented_LHlift(double, direction);
  void apply_oriented_HLlift(double, direction);
	void oriented_analysis(direction);
	void oriented_synthesis(direction);
	void oriented_packet_analysis(direction);
	void oriented_packet_synthesis(direction);
  // defined in hpfprelift.cpp
  void hpf_HLlift(double a, direction dir, bool adaptive);
  void hpf_oriented_analysis(direction, bool adaptive=true);
  void hpf_oriented_synthesis(direction, bool adaptive=true);
  // experiment testing functions (defined in experiment.cpp)
  friend double mse(dwtnode &a, dwtnode &b);
  void shift(int sigma);
  void halveimage(bool);
  void call_batch(testmode mode, bool halfres);
  void rawl_encode(bool halfres=false, bool adapt=false);
  void antialias_encode(bool halfres=false, bool adapt=false);
  void orient_encode(bool halfres=false, bool adapt=false);
  void aa_orient_encode( bool halfres=false, bool adapt=false);
  void rawl_decode(char *bitrate, bool adapt=false);
  void antialias_decode(char *bitrate, bool adapt=false);
  void orient_decode(char *bitrate, bool adapt=false);
  void aa_orient_decode(char *bitrate, bool adapt=false);
  void orient2_encode(bool out=false, bool est=false);
  void orient2_decode(char *bitrate,  bool est=false);
  void aa_orient2_encode(bool out=false, bool est=false);
  void aa_orient2_decode(char *bitrate, bool est=false);
	void hpfprelift_encode(bool halfres=false, bool adapt=true);
	void hpfprelift_decode(char *bitrate, bool adapt=true);
// data members
protected:
  int h,w; // height and width of the image
  // records the depth of analysis that has been carried out so far
  // on the image data. assignment/retrieval of dwtlevel should be
  // done through indexing dwtlevel[vertical] and dwtlevel[horizontal]
  // as appropriate.
  int dwtlevel[2];
  enum dwttype dwtbase; // indicates the wavelet filters used, eg w5x3 or w9x7
  double *pixels;
public:
  dwtnode *subbands[4]; // order of subbands: LL,HL,LH,HH
  orientationfield ofield;
};
class estorient : public dwtnode
{ // functions defined in orient_est.cpp
public:
  estorient(char* fname,dwttype type):dwtnode(fname,type){}
  estorient(char* fname,int hset,int wset,dwttype type):dwtnode(fname,hset,wset,type){}
  estorient(dwtnode *target);
  void init_orient(int setblksz, int setprec, int setmaxshift);
  void transpose();
  void calc_energies();
  double fetch_residual(int orientblk, direction dir, char shift);
  bool test_residual(int orientblk, direction dir, char shift, double &best);
  void choose_shift(int orientblk, direction dir, double thresh);
  void choose_orient();
  void legacy_choose_orient();  
protected:
  struct tuple{double **hJ; double **vJ;} orientenergy;
};