#include "stdafx.h"
#include "dwtnode.h"
int testreconstruction(int argc, _TCHAR* argv[])
{
  shker shftbase("sideinf\\baseker.dat"); // create shift lut objects
  shker shftpoly2("sideinf\\v2poly.dat");
  shker shftpoly4("sideinf\\v4packnormalised.dat");
  packlift_filters filts("sideinf\\icip_aa_filters.dat");
  //char currfile[] = "C:\\Data\\Images\\nonstandard\\thin_horz_gen.pgm";
  char currfile[] = "C:\\Data\\Images\\barbara.pgm";
  //char currfile[] = "C:\\Data\\Images\\nonstandard\\nyquist_horz_gen.pgm";
  dwtnode in(currfile,w9x7);
  //in.ofield.init_orient("sideinf\\horz_generated.dat");
  //in.ofield.init_orient("sideinf\\horz_generated_wrong.dat");
  //in.ofield.init_orient("sideinf\\barb4.dat");

  // option 2: reconstructs LL1 using only its child subbands
  in.oriented_packet_analysis(shftbase,shftpoly4);
  in.extract_subband(0);
  in.subbands[0]->ofield.inherit(in.ofield);
  in.subbands[0]->oriented_synthesis(shftbase,shftpoly2);
  in.subbands[0]->pgmwrite("barb_rightleg.pgm");

  //in.oriented_packet_analysis(shftbase,shftpoly4);
  //in.operating_depth++;
  //in.oriented_synthesis(shftbase,shftpoly4);
  //in.pgmwrite("test_option1.pgm");

  //in.packet_synthesis(both);
  //in.rawlwrite("tmp\\out.rawl");
	return 0;
}
int antialiasing(int argc, _TCHAR* argv[])
{
  packlift_filters filts("sideinf\\icip_aa_filters.dat");
  char currfile[] = "C:\\Data\\Images\\barbara.pgm";
  //char currfile[] = "C:\\Data\\Images\\nonstandard\\nyquist_horz_gen.pgm";
  dwtnode in(currfile,w9x7);
  in.analysis(both);
  in.extract_subband(0);
  in.extract_subband(1);
  in.extract_subband(2);
  in.subbands[0]->analysis(both);
  in.subbands[1]->analysis(both);
  in.subbands[2]->analysis(both);
  in.packlift(both,true,true); // adaptive packet lifting
  in.subbands[0]->synthesis(both);
  in.subbands[0]->pgmwrite("adaptiveantialias.pgm");

  //// reverse steps to check reconstruction
  //in.packlift(both,false,true);
  //in.subbands[0]->synthesis(both);
  //in.subbands[1]->synthesis(both);
  //in.subbands[2]->synthesis(both);
  //in.interleave();
  //in.synthesis(both);
  //dwtnode comparison(currfile);
  //cout << mse(in,comparison);
  return 0;
}
int compresstest(int argc, _TCHAR* argv[])
{
  enum testmode {base,pyramid3x2,pyramid2x3,packlift,orient,aaorient};
  shker shft1("sideinf\\baseker.dat"); // create shift lut objects
  shker shft2("sideinf\\v2poly.dat");
  shker shft4("sideinf\\v4packnormalised.dat");
  char currfile[] = "C:\\Data\\Images\\city0.pgm";
  ofstream dout("results\\dumpout.txt",ios::app);
  if (false){ // orientation estimation flow
    estorient est(currfile,w9x7);
    est.init_orient(4,8,16);
    est.calc_energies(shft1);
    //est.legacy_choose_orient();
    est.choose_orient();
    est.ofield.orientwrite("sideinf\\tmp.dat");
    est.ofield.orient_csvout();
  }
  if (false){
    est2orient est2(currfile,w9x7);
    est2.init_orient(4,8,16);
    est2.choose_orient();
    est2.ofield.orientwrite("sideinf\\tmp2.dat");
    est2.ofield.orient_csvout();
  }
  dwtnode ref(currfile,w9x7);
  dwtnode in(currfile,w9x7);
  in.load_packfilts("sideinf\\icip_aa_filters.dat");
  //in.ofield.init_orient("sideinf\\barb4.dat");
  //in.ofield.setaffinefield();
  bool adapt=false; // select adaptive mode
  bool halfres=true; // compute mse for half resolution recon instead
  bool imageout=true; // dump out reconstructions and collate
  testmode mode=orient;
  std::stringstream batch; // set batch file and resolution arguments
  batch << (halfres?"halfres.bat":"out.bat")<<" "<<ref.geth()<<" "<<ref.getw();
  void (dwtnode::*encode_ptr)(shker &,shker &,shker &,bool,bool) = NULL;
  void (dwtnode::*decode_ptr)(char *,shker &, shker &, shker &, bool) = NULL;
  switch (mode)
  {
  case base:
    encode_ptr = &dwtnode::rawl_encode;
    dout << "Ordinary dwt MSEs";
    decode_ptr = &dwtnode::rawl_decode;
    if (halfres)
    {
      dout << " half resolution";
      ref.analysis(both);
    }
    dout << endl;
    break;
  case packlift:
    encode_ptr = &dwtnode::antialias_encode;
    decode_ptr = &dwtnode::antialias_decode;
    dout << "Antialiased MSEs";
    if (halfres)
    {
      dout << " half resolution";
      decode_ptr = &dwtnode::rawl_decode;
      ref.load_packfilts("sideinf\\icip_aa_filters.dat");
      ref.analysis(both);
      ref.extract_subband(0);
      ref.extract_subband(1);
      ref.extract_subband(2);
      ref.subbands[0]->analysis(both);
      ref.subbands[1]->analysis(both);
      ref.subbands[2]->analysis(both);
      ref.packlift(both,true,adapt);
      ref.subbands[0]->synthesis(both);
      ref.interleave(true);
    }
    if (adapt) dout << " adaptive";
    dout << endl;
    break;
  case pyramid3x2:
    encode_ptr = &dwtnode::lp3x2_encode;
    dout << "Flierl pyramid MSEs";
    decode_ptr = &dwtnode::lp3x2_decode;
    if (halfres)
    {
      dout << " half resolution";
      ref.lp3x2_halfres(); // replace image with laplacian half res
      decode_ptr = &dwtnode::rawl_decode;
    }
    else
    {
      batch.str("");
      batch << "pyramid.bat " << ref.geth()<<" "<<ref.getw();
    }
    dout << endl;
    break;
  case pyramid2x3:
    encode_ptr = &dwtnode::lp2x3_encode;
    dout << "Tran pyramid MSEs";
    decode_ptr = &dwtnode::lp2x3_decode;
    if (halfres)
    {
      dout << " half resolution";
      ref.lp2x3_halfres(); // replace image with laplacian half res
      decode_ptr = &dwtnode::rawl_decode;
    }
    else
    {
      batch.str("");
      batch << "pyramid.bat " << ref.geth()<<" "<<ref.getw();
    }
    dout << endl;
    break;
  case orient:
    encode_ptr = &dwtnode::orient_encode;
    dout << "Oriented wavelet MSEs";
    if (halfres)
    {
      dout << " half resolution";
      decode_ptr = &dwtnode::rawl_decode;
      ref.oriented_analysis(shft1,shft2);
      cerr << "not yet implemented!";
      exit(1);
    }
    else  decode_ptr = &dwtnode::orient_decode;
    dout << endl;
    break;
  case aaorient:
    encode_ptr = &dwtnode::aa_orient_encode;
    dout << "Oriented + antialiased MSEs" << endl;
    if (halfres)
    {
      dout << " half resolution";
      decode_ptr = NULL;
      cerr << "not yet implemented!";
      exit(1);
    }
    else  decode_ptr = &dwtnode::aa_orient_decode;
    dout << endl;
    break;
  }
  {
    (in.*encode_ptr)(shft1,shft2,shft4,halfres,adapt);
    system(batch.str().c_str());
    in.halveimage(halfres);
    //// test perfect reconstruction
    (in.*decode_ptr)("",shft1,shft2,shft4,adapt);
    cout << mse(ref,in) << endl;
    //'(' << 10*log10(255^2/mse) << ')'; // output psnr instead
    (in.*decode_ptr)("0.1",shft1,shft2,shft4,adapt);
    if (imageout) in.pgmwrite("tmp\\recon0.1.pgm");
    dout << mse(ref,in) << ' ';
    (in.*decode_ptr)("0.2",shft1,shft2,shft4,adapt);
    if (imageout) in.pgmwrite("tmp\\recon0.2.pgm");
    dout << mse(ref,in) << ' ';
    (in.*decode_ptr)("0.4",shft1,shft2,shft4,adapt);
    if (imageout) in.pgmwrite("tmp\\recon0.4.pgm");
    dout << mse(ref,in) << ' ';
    (in.*decode_ptr)("0.6",shft1,shft2,shft4,adapt);
    if (imageout) in.pgmwrite("tmp\\recon0.6.pgm");
    dout << mse(ref,in) << ' ';
    (in.*decode_ptr)("0.8",shft1,shft2,shft4,adapt);
    if (imageout) in.pgmwrite("tmp\\recon0.8.pgm");
    dout << mse(ref,in) << ' ';
    (in.*decode_ptr)("1.0",shft1,shft2,shft4,adapt);
    if (imageout)
    {
      if (halfres)
        ref.pgmwrite("tmp\\halfres.pgm");
      in.pgmwrite("tmp\\recon1.0.pgm");
      system("kdu_compress -i tmp\\recon0.1.pgm,tmp\\recon0.2.pgm,tmp\\recon0.4.pgm,"
        "tmp\\recon0.6.pgm,tmp\\recon0.8.pgm,tmp\\recon1.0.pgm"
        " -o tmp\\recon.jp2 Creversible=yes Cycc=no");
    }
    dout << mse(ref,in) << endl;
  }
  dout.close();
  return 0;
}
int hpftest(int argc, _TCHAR* argv[])
{
  shker shft1("sideinf\\baseker.dat"); // create shift lut objects
  shker shft2("sideinf\\v2poly.dat");
  char currfile[] = "C:\\Data\\Images\\nonstandard\\vert_generated.pgm";
  //char currfile[] = "C:\\Data\\Images\\barbara.pgm";
  //char currfile[] = "dctest.pgm";
  dwtnode in(currfile,w5x3);
  //in.ofield.init_orient("sideinf\\barb4.dat");
  //in.ofield.init_orient(4,8,8);
  in.ofield.init_orient(4,8,8,4,0);
  in.ofield.setaffinefield();

  //in.transpose();
  //in.analysis(vertical);
  //in.hpf_oriented_analysis(shft2,horizontal);
  //in.transpose();
  //in.pgmwrite("");

  //in.analysis(both);
  //in.oriented_analysis(shft1,shft2);
  in.hpf_oriented_analysis(shft1,shft2);

  in.rawlwrite("C:\\Program Files\\Matlab7\\work\\temp\\test_hpforient.rawl",6,true);
  return 0;
}
int _tmain(int argc, _TCHAR* argv[])
{
  //system("del sideinf\\alpha_transfer.dat");
  //system("del sideinf\\alpha_cancel.dat");
  //system("del sideinf\\cancel_energy.dat");
  //compresstest(argc,argv);
  //antialiasing(argc,argv);
  hpftest(argc,argv);
  return 0;
}