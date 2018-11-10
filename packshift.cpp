#include "stdafx.h"
#include "dwtnode.h"
int antialiasing(int argc, _TCHAR* argv[])
{
  packlift_filters filts("sideinf\\icip_aa_filters.dat");
  char currfile[] = "D:\\Work\\Images\\barbara.pgm";
  //char currfile[] = "D:\\Work\\Images\\nonstandard\\nyquist_horz_gen.pgm";
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
  char currfile[] = "D:\\Work\\Images\\city0.pgm";
  ofstream dout("results\\dumpout.txt",ios::app);
  if (false){ // orientation estimation flow
    estorient est(currfile,w9x7);
    est.init_orient(4,8,16);
    est.calc_energies();
    //est.legacy_choose_orient();
    est.choose_orient();
    est.ofield.orientwrite("sideinf\\tmp.dat");
    est.ofield.orient_csvout();
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
  void (dwtnode::*encode_ptr)(bool,bool) = NULL;
  void (dwtnode::*decode_ptr)(char *, bool) = NULL;
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
      ref.oriented_analysis(both);
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
    (in.*encode_ptr)(halfres,adapt);
    system(batch.str().c_str());
    in.halveimage(halfres);
    //// test perfect reconstruction
    (in.*decode_ptr)("",adapt);
    cout << mse(ref,in) << endl;
    //'(' << 10*log10(255^2/mse) << ')'; // output psnr instead
    (in.*decode_ptr)("0.1",adapt);
    if (imageout) in.pgmwrite("tmp\\recon0.1.pgm");
    dout << mse(ref,in) << ' ';
    (in.*decode_ptr)("0.2",adapt);
    if (imageout) in.pgmwrite("tmp\\recon0.2.pgm");
    dout << mse(ref,in) << ' ';
    (in.*decode_ptr)("0.4",adapt);
    if (imageout) in.pgmwrite("tmp\\recon0.4.pgm");
    dout << mse(ref,in) << ' ';
    (in.*decode_ptr)("0.6",adapt);
    if (imageout) in.pgmwrite("tmp\\recon0.6.pgm");
    dout << mse(ref,in) << ' ';
    (in.*decode_ptr)("0.8",adapt);
    if (imageout) in.pgmwrite("tmp\\recon0.8.pgm");
    dout << mse(ref,in) << ' ';
    (in.*decode_ptr)("1.0",adapt);
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
  //char currfile[] = "D:\\Work\\Images\\barbara.pgm";
  char currfile[] = "D:\\Work\\Images\\nonstandard\\vert_generated.pgm";
  dwtnode in(currfile,w5x3);
  //in.ofield.init_orient("sideinf\\barb4.dat");
  in.ofield.init_orient(4,8,8,4,0);
  in.ofield.setaffinefield();

  //in.analysis(both);
  in.oriented_analysis(both);
	//in.synthesis(both);

	in.pgmwrite("base.pgm");

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