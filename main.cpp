#include "stdafx.h"
#include "dwtnode.h"
#include <assert.h>
int estimate(int argc, char* argv[])
{
	char currfile[] = "D:\\Work\\Images\\lena.pgm";
  estorient est(currfile,w5x3);
  est.init_orient(4,8,16);
  est.calc_energies();
  //est.legacy_choose_orient();
  est.choose_orient();
  est.ofield.orientwrite("sideinf\\lena.dat");
	return 0;
}
int quadtree_estimate(int argc, char* argv[])
{
  char currfile[] = "D:\\Work\\Images\\barbclean.pgm";
  estorient2 est(currfile,w5x3);
  est.ofield.init_orient(1,8,8,0,0);
  est.calc_energies(1,400);
  est.quadtree_estimate();
  est.quadtree_flatten();
  est.ofield.orientwrite("sideinf\\quadtree.dat");
  return 0;
}
int yuvstreamprocess(int argc, char* argv[])
{
  char fname[] = "D:\\Work\\Videos\\MOBILE_352x288_30_orig_01.yuv";
  std::ifstream fin(fname,std::ios::binary);
  std::ofstream half("mobile_176x144.yuv",std::ios::binary);
  dwtnode node(288,352,w5x3);
  node.yuvstreamread(fin);
  node.pgmwrite("mobile0.pgm");
  //while (node.yuvstreamread(fin))
  //{
  //  //cout << '.';
  //  node.analysis(both);
  //  node.yuvstreamwrite(half);
  //  //node.pgmwrite("test.pgm");
  //}
  return 0;
}
int shifttest(int argc, char* argv[])
{
  char currfile[] = "D:\\Work\\Images\\city0.pgm";
  std::ofstream yuvout("shifttest.yuv",std::ios::binary);
  if (!yuvout.good())
  {
    std::cerr << "Access of file unsuccessful." << std::endl;
    exit(1);
  }
  for (int sigma=0;sigma<80;sigma++)
  {
    dwtnode in(currfile,w5x3);
    in.ofield.init_orient(4,8,32,0,0);
    in.analysis(horizontal);
    in.analysis(horizontal);
    in.shift(sigma,horizontal);
    in.synthesis(horizontal);
    in.synthesis(horizontal);
    in.yuvstreamwrite(yuvout);
  }
  yuvout.close();
  return 0;
}
int orienttest(int argc, char* argv[])
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
		  in.oriented_synthesis(both);
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
int compresstest(int argc, char* argv[])
{
  char Cdecomp[] = "B(BH-H-:BVV--:-),B(H:V:B)";
  //char Cdecomp[] = "B(-:-:-),B(-:-:-)";
  char currfile[] = "D:\\Work\\Images\\bikecrop.pgm";
  std::ofstream dout("results\\dumpout.txt",std::ios::app);
  txtype txtest = w5x3;
  bool adapt=true; // select adaptive mode
  bool halfres=false; // compute mses for half resolution instead
  int depth=2; // number of layers of scalability to create
  int layer=0; // compute mses for encoding layer (0 is full resolution)
  bool imageout=(layer>0); // dump out compressed, decoded images and collate
  testmode mode=orient;
  dwtnode in(currfile,txtest);
  in.ofield.init_orient("sideinf\\bikecrop_11241.dat");
	in.ofield.setaffinefield();
	dwtnode ref=in;
  void (dwtnode::*encode_ptr)(int, int, bool) = NULL;
  void (dwtnode::*decode_ptr)(char *, int, int, bool) = NULL;
  assert(layer<=depth);
  switch (mode)
  {
  case base:{
    encode_ptr = &dwtnode::rawl_encode;
    dout << "Ordinary dwt MSEs";
    decode_ptr = &dwtnode::rawl_decode;
    for (int i=0;i<layer;i++)
      ref.analysis(both);
    break;}
  case pyramid_test:{
    encode_ptr = &dwtnode::pyramid_encode;
    decode_ptr = &dwtnode::pyramid_decode;
    if (txtest == pyramid)
      dout << "Laplacian pyramid MSEs";
    else if (txtest == pyramid3x2)
      dout << "Flierl pyramid MSEs";
    else if (txtest == pyramid2x3)
      dout << "Tran pyramid MSEs";
    else
    {
      std::cerr << "Non-pyramid transform selected with pyramid encode test" << std::endl;
      exit(1);
    }
    for (int i=0;i<layer;i++)
      ref.pyramid_analysis();
    break;}
  case packlift:{
    encode_ptr = &dwtnode::packlift_encode;
    dout << "Packet lifting MSEs";
    decode_ptr = &dwtnode::packlift_decode;
    if (halfres)
    {
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
    break;}
  //case packlift_2layer:{
  //  encode_ptr = &dwtnode::packlift_2layer_encode;
  //  dout << "2 level packet lifting MSEs";
  //  decode_ptr = &dwtnode::packlift_2layer_decode;
  //  if (halfres){
  //    cerr << "quarter res currently unsupported" << endl;
  //    exit(1);
  //  }
  //  if (adapt) dout << " adaptive";
  //  break;}
  //case packliftorient2:{
  //  encode_ptr = &dwtnode::packlift_orient2_encode;
  //  dout << "2 level oriented + packlift MSEs";
  //  decode_ptr = &dwtnode::packlift_orient2_decode;
  //  if (halfres)
  //  {
  //    ref.oriented_analysis(both);
  //    ref.extract_subband(0);
  //    ref.extract_subband(1);
  //    ref.extract_subband(2);
  //    ref.subbands[0]->analysis(both);
  //    ref.subbands[1]->analysis(both);
  //    ref.subbands[2]->analysis(both);
  //    ref.packlift(both,true,adapt);
  //    ref.subbands[0]->synthesis(both);
  //    ref.interleave(true);
  //  }
  //  break;}
  case orient:{
    encode_ptr = &dwtnode::orient_encode;
    dout << "Oriented wavelet MSEs";
		decode_ptr = &dwtnode::orient_decode;
    dwtnode *curr = &ref;
    for (int i=0;i<layer;i++,curr=curr->subbands[0])
    {
      curr->oriented_analysis(both);
      curr->extract_subband(0);
    }
    break;}
  //case orient2packet:{
  //  encode_ptr = &dwtnode::orient2_packet_encode;
  //  dout << "Oriented packet wavelet MSEs";
		//decode_ptr = &dwtnode::orient2_packet_decode;
  //  if (halfres)
  //  {
  //    ref.oriented_packet_analysis(both);
  //    ref.extract_subband(0);
  //    ref.subbands[0]->oriented_synthesis(both);
  //    dwtnode tmp(*ref.subbands[0]);
  //    ref = tmp;
  //  }
  //  break;}
  case hpfprelift:{
    encode_ptr = &dwtnode::hpfprelift_encode;
    dout << "HPF prelift MSEs";
    decode_ptr = &dwtnode::hpfprelift_decode;
    dwtnode *curr = &ref;
    for (int i=0;i<layer;i++,curr=curr->subbands[0])
    {
      curr->hpf_oriented_analysis(both,adapt);
      curr->extract_subband(0);
    }
    break;}
  default:
    std::cerr << "This test mode is unsupported" << std::endl;
    exit(1);
  }
  {
    (in.*encode_ptr)(depth,layer,adapt);
    in.call_batch(mode,Cdecomp,depth,layer,dout); // run kdu_compress
    in.shrink(layer); // reduce dimensions if layer>0
    //// test perfect reconstruction
    (in.*decode_ptr)("",depth,layer,adapt);
    std::cout << "testing perfect reconstruction: " << mse(ref,in) << std::endl;
    //'(' << 10*log10(255^2/mse) << ')'; // output psnr instead
    (in.*decode_ptr)("0.1",depth,layer,adapt);
    if (imageout) in.pgmwrite("tmp\\recon0.1.pgm");
    dout << mse(ref,in) << ' ';
    (in.*decode_ptr)("0.2",depth,layer,adapt);
    if (imageout) in.pgmwrite("tmp\\recon0.2.pgm");
    dout << mse(ref,in) << ' ';
    (in.*decode_ptr)("0.4",depth,layer,adapt);
    if (imageout) in.pgmwrite("tmp\\recon0.4.pgm");
    dout << mse(ref,in) << ' ';
    (in.*decode_ptr)("0.6",depth,layer,adapt);
    if (imageout) in.pgmwrite("tmp\\recon0.6.pgm");
    dout << mse(ref,in) << ' ';
    (in.*decode_ptr)("0.8",depth,layer,adapt);
    if (imageout) in.pgmwrite("tmp\\recon0.8.pgm");
    dout << mse(ref,in) << ' ';
    (in.*decode_ptr)("1.0",depth,layer,adapt);
    dout << mse(ref,in) << std::endl;
    if (imageout)
    {
      if (layer==1)
        ref.pgmwrite("results\\halfres.pgm");
      if (layer==2)
        ref.pgmwrite("results\\quartres.pgm");
      in.pgmwrite("tmp\\recon1.0.pgm");
      system("kdu_compress -i tmp\\recon0.1.pgm,tmp\\recon0.2.pgm,tmp\\recon0.4.pgm,"
        "tmp\\recon0.6.pgm,tmp\\recon0.8.pgm,tmp\\recon1.0.pgm"
        " -o tmp\\recon.jp2 Creversible=yes Cycc=no");
    }
  }
  dout.close();
  return 0;
}
int hpftest(int argc, char* argv[])
{
  char currfile[] = "D:\\Work\\Images\\barbara.pgm";
  //char currfile[] = "D:\\Work\\Images\\nonstandard\\vert_generated.pgm";
	//char currfile[] = "D:\\Work\\Images\\city0.pgm";
  dwtnode in(currfile,w5x3);
  in.ofield.init_orient("sideinf\\barb4.dat");
  //in.ofield.init_orient(4,8,8,4,0);
  in.ofield.setaffinefield();

	in.hpf_oriented_analysis(both,true);
	in.pgmwrite("test.pgm");
  return 0;
}
int main(int argc, char* argv[])
{
  //system("del sideinf\\alpha_transfer.dat");
  compresstest(argc,argv);
  //orienttest(argc,argv);
  //quadtree_estimate(argc,argv);
  //hpftest(argc,argv);
  //yuvstreamprocess(argc,argv);
  return 0;
}