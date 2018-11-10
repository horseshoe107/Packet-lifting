#include "stdafx.h"
#include "dwtnode.h"
int estimate(int argc, _TCHAR* argv[])
{
	//char currfile[] = "D:\\Work\\Images\\barbara.pgm";
	char currfile[] = "D:\\Work\\Images\\nonstandard\\horz_generated.pgm";
  estorient est(currfile,w5x3);
  est.init_orient(4,8,16);
  est.calc_energies();
  //est.legacy_choose_orient();
  est.choose_orient();
  est.ofield.orientwrite("sideinf\\test.dat");
	return 0;
}
int compresstest(int argc, _TCHAR* argv[])
{
  char currfile[] = "D:\\Work\\Images\\barbara.pgm";
  ofstream dout("results\\dumpout.txt",ios::app);
  dwtnode in(currfile,w5x3);
  in.ofield.init_orient("sideinf\\barb4.dat");
	in.ofield.setaffinefield();
	dwtnode ref=in;
  bool adapt=true; // select adaptive mode
  bool halfres=false; // compute mses for half resolution instead
  bool imageout=true; // dump out compressed, decoded images and collate
  testmode mode=orient;
  void (dwtnode::*encode_ptr)(bool, bool) = NULL;
  void (dwtnode::*decode_ptr)(char *, bool) = NULL;
  switch (mode)
  {
  case base:
    encode_ptr = &dwtnode::rawl_encode;
    dout << "Ordinary dwt MSEs";
    decode_ptr = &dwtnode::rawl_decode;
    if (halfres)
    {
      ref.analysis(both);
    }
    break;
  case packlift:
    encode_ptr = &dwtnode::antialias_encode;
    decode_ptr = &dwtnode::antialias_decode;
    dout << "Packet lifting MSEs";
    if (halfres)
    {
      decode_ptr = &dwtnode::rawl_decode;
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
    break;
  case pyramid3x2:
    encode_ptr = &dwtnode::lp3x2_encode;
    dout << "Flierl pyramid MSEs";
    decode_ptr = &dwtnode::lp3x2_decode;
    if (halfres)
    {
      ref.lp3x2_halfres(); // replace image with laplacian half res
      decode_ptr = &dwtnode::rawl_decode;
    }
    break;
  case pyramid2x3:
    encode_ptr = &dwtnode::lp2x3_encode;
    dout << "Tran pyramid MSEs";
    decode_ptr = &dwtnode::lp2x3_decode;
    if (halfres)
    {
      ref.lp2x3_halfres(); // replace image with laplacian half res
      decode_ptr = &dwtnode::rawl_decode;
    }
    break;
  case orient:
    encode_ptr = &dwtnode::orient_encode;
    dout << "Oriented wavelet MSEs";
		decode_ptr = &dwtnode::orient_decode;
    if (halfres)
    {
      decode_ptr = &dwtnode::rawl_decode;
      ref.oriented_analysis(both);
      cerr << " half res not yet implemented!";
      exit(1);
    }
    break;
  case aaorient:
    encode_ptr = &dwtnode::aa_orient_encode;
    dout << "Oriented + antialiased MSEs";
    decode_ptr = &dwtnode::aa_orient_decode;
    if (halfres)
    {
      decode_ptr = NULL;
      cerr << " half res not yet implemented!";
      exit(1);
    }
    break;
  default:
    exit(1);
  }
  {
    (in.*encode_ptr)(halfres,adapt);
    in.call_batch(mode,halfres,dout); // run kdu_compress
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
int hpfcompresstest(int argc, _TCHAR* argv[])
{
  char currfile[] = "D:\\Work\\Images\\barbara.pgm";
  ofstream dout("results\\dumpout.txt",ios::app);
  dwtnode in(currfile,w5x3);
  in.ofield.init_orient("sideinf\\barb4.dat");
	in.ofield.setaffinefield();
	dwtnode ref=in;
  bool adapt=true; // select adaptive mode
  bool halfres=false; // compute mses for half resolution instead
  bool imageout=true; // dump out compressed, decoded images and collate
  testmode mode=packlift;
  void (dwtnode::*encode_ptr)(bool, bool) = NULL;
  void (dwtnode::*decode_ptr)(char *, bool) = NULL;
  switch (mode)
  {
  case base:
    encode_ptr = &dwtnode::rawl_encode;
    dout << "Ordinary dwt MSEs";
    decode_ptr = &dwtnode::rawl_decode;
    if (halfres)
    {
      ref.analysis(both);
    }
    break;
  case packlift:
    encode_ptr = &dwtnode::antialias_encode;
    decode_ptr = &dwtnode::antialias_decode;
    dout << "Packet lifting MSEs";
    if (halfres)
    {
      decode_ptr = &dwtnode::rawl_decode;
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
    break;
  case orient:
    encode_ptr = &dwtnode::orient_encode;
    dout << "Oriented wavelet MSEs";
		decode_ptr = &dwtnode::orient_decode;
    if (halfres)
    {
      decode_ptr = &dwtnode::rawl_decode;
      ref.oriented_analysis(both);
      cerr << " half res not yet implemented!";
      exit(1);
    }
    break;
  case hpfprelift:
    encode_ptr = &dwtnode::hpfprelift_encode;
    dout << "HPF prelift MSEs";
    decode_ptr = &dwtnode::hpfprelift_decode;
    if (halfres)
    {
      decode_ptr = &dwtnode::rawl_decode;
      ref.hpf_oriented_analysis(both,adapt);
    }
    break;
  default: exit(1);
  }
  {
    (in.*encode_ptr)(halfres,adapt);
    in.call_batch(mode,halfres,dout); // run kdu_compress
    in.halveimage(halfres);
    //// test perfect reconstruction
    (in.*decode_ptr)("",adapt);
    cout << "testing perfect reconstruction: " << mse(ref,in) << endl;
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
  char currfile[] = "D:\\Work\\Images\\barbara.pgm";
  //char currfile[] = "D:\\Work\\Images\\nonstandard\\vert_generated.pgm";
	//char currfile[] = "D:\\Work\\Images\\city0.pgm";
  dwtnode in(currfile,w5x3);
  in.ofield.init_orient("sideinf\\barb4.dat");
  //in.ofield.init_orient(4,8,8,4,0);
  in.ofield.setaffinefield();

	in.hpf_oriented_analysis(vertical,true);
  //in.oriented_analysis(vertical);
  in.oriented_analysis(horizontal);
	in.pgmwrite("hpf_everywhere_updatewithoutaliasing.pgm");

  return 0;
}
int orienttest(int argc, _TCHAR* argv[])
{
	enum testmode {orient, recon1, recon2, recon3} mode;
  char currfile[] = "D:\\Work\\Images\\barbara.pgm";
	//char currfile[] = "D:\\Work\\Images\\nonstandard\\vert_generated.pgm";
  for (int i=0;i<3;i++)
  {
    mode = (testmode)i;
    dwtnode in(currfile,w5x3);
    in.ofield.init_orient("sideinf\\barb4.dat");
    //in.ofield.init_orient(4,8,8,4,0);
    in.ofield.setaffinefield();
    switch (mode)
    {
    case orient:
		  in.oriented_analysis(both);
		  in.pgmwrite("LL1.pgm");
      break;
	  case recon1: // reconstruct w/ all subbands
		  in.oriented_packet_analysis(both);
		  in.oriented_synthesis(horizontal);
		  in.oriented_synthesis(vertical);
		  in.pgmwrite("LL1recon1.pgm");
		  break;
	  case recon2: // reconstruct only w/ LL1
		  in.oriented_packet_analysis(both);
		  in.extract_subband(0);
		  in.subbands[0]->oriented_synthesis(both);
		  in.subbands[0]->pgmwrite("LL1recon2.pgm");
		  break;
	  case recon3: // reconstruct w/ LL1 and lower res ofield
		  break;
    }
  }
	return 0;
}
int _tmain(int argc, _TCHAR* argv[])
{
  //system("del sideinf\\alpha_transfer.dat");
  //system("del sideinf\\alpha_cancel.dat");
  //system("del sideinf\\cancel_energy.dat");
  //compresstest(argc,argv);
  orienttest(argc,argv);
	//estimate(argc,argv);
  //hpftest(argc,argv);
  //hpfcompresstest(argc,argv);
  return 0;
}