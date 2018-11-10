#include "mex.h"

int divc(int,int);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // define local variables
    int spacing, offset, length, rows, cols, x, y1, y2;
    double *in, *out;
    double *ivals, *ovals;
    mxArray *ivalsM, *ovalsM;
    
    in = mxGetPr(prhs[0]); // define input matrix
    rows = mxGetM(prhs[0]); // number of rows
    cols = mxGetN(prhs[0]); // number of cols
    spacing = (int) *mxGetPr(prhs[1]);
    offset = (int) *mxGetPr(prhs[2]);
    plhs[0] = mxDuplicateArray(prhs[0]);
    out = mxGetPr(plhs[0]);
    
    // initialise ivals matrix
    length = divc((rows-offset),spacing); // true length of desired vector
    ivalsM = mxCreateDoubleMatrix(length+8,cols,mxREAL);
    ivals = mxGetPr(ivalsM);
    // copy indexed component of in -> ivals
    for (x=0;x<cols;x++) // going along columns
        for (y1=4,y2=offset;y2<rows;y1++,y2+=spacing)
            ivals[x*(length+8)+y1] = in[x*rows+y2];
    // symmetric extension
    for (x=0;x<cols;x++)
    {
        for (y1=1;y1<=4;y1++)
        {
            ivals[x*(length+8)+4-y1] = ivals[x*(length+8)+4+y1];
            ivals[x*(length+8)+(3+length)+y1] = ivals[x*(length+8)+(3+length)-y1];
        }
    }
    // initialise ovals matrix
    ovalsM = mxCreateDoubleMatrix(length+16,cols,mxREAL);
    ovals = mxGetPr(ovalsM);
    // apply low pass synthesis basis vector
    for (x=0;x<cols;x++)
    {
        for (y1=0;y1<(length+8);y1+=2)
        {
            ovals[x*(length+16)+y1+4] += 2*0.557543*ivals[x*(length+8)+y1];
            ovals[x*(length+16)+y1+4-1] += 2*0.295636*ivals[x*(length+8)+y1];
            ovals[x*(length+16)+y1+4+1] += 2*0.295636*ivals[x*(length+8)+y1];
            ovals[x*(length+16)+y1+4-2] -= 2*0.028772*ivals[x*(length+8)+y1];
            ovals[x*(length+16)+y1+4+2] -= 2*0.028772*ivals[x*(length+8)+y1];
            ovals[x*(length+16)+y1+4-3] -= 2*0.045636*ivals[x*(length+8)+y1];
            ovals[x*(length+16)+y1+4+3] -= 2*0.045636*ivals[x*(length+8)+y1];
        }
    }
    // apply high pass synthesis basis vector
    for (x=0;x<cols;x++)
    {
        for (y1=1;y1<(length+8);y1+=2)
        {
            ovals[x*(length+16)+y1+4] += 2*0.602949*ivals[x*(length+8)+y1];
            ovals[x*(length+16)+y1+4-1] -= 2*0.266864*ivals[x*(length+8)+y1];
            ovals[x*(length+16)+y1+4+1] -= 2*0.266864*ivals[x*(length+8)+y1];
            ovals[x*(length+16)+y1+4-2] -= 2*0.078223*ivals[x*(length+8)+y1];
            ovals[x*(length+16)+y1+4+2] -= 2*0.078223*ivals[x*(length+8)+y1];
            ovals[x*(length+16)+y1+4-3] += 2*0.016864*ivals[x*(length+8)+y1];
            ovals[x*(length+16)+y1+4+3] += 2*0.016864*ivals[x*(length+8)+y1];
            ovals[x*(length+16)+y1+4-4] += 2*0.026749*ivals[x*(length+8)+y1];
            ovals[x*(length+16)+y1+4+4] += 2*0.026749*ivals[x*(length+8)+y1];
        }
    }
    // write pertinent section of ovals to out
    for (x=0;x<cols;x++)
        for (y1=8,y2=offset;y2<rows;y1++,y2+=spacing)
            out[x*(rows)+y2] = ovals[x*(length+16)+y1];
    // cleanup
    mxDestroyArray(ivalsM);
    mxDestroyArray(ovalsM);
    return;
}


// integer division where the quotient is rounded up
int divc(int num, int den)
{
    int quotient, sign=1;
    if (num<0)
    { // force only positive num and den
        num=-num;
        sign=-sign;
    }
    if (den<0)
    {
        den=-den;
        sign=-sign;
    }
    for (quotient=0;num>0;num-=den,quotient++) {}
    return (sign*quotient);
}
