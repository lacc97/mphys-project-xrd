//
// Copyright © 2001-2007 Rigaku Americas Corp.
//                       9009 New Trails Drive
//                       The Woodlands, TX, USA  77381
//
// The contents are unpublished proprietary source
// code of Rigaku Americas
//
// All rights reserved
//
// ReadRAXISImage.cpp    Initial author: T.L.Hendrixson           Aug 2001
//
//
// Description:
//
// This program reads in an R-AXIS image written using the Rigaku-style
// format and writes the values of a few of the pixels to standard output.
// It additionally writes out some of the information contained in the
// header, such as axis positions, exposure time, pixel information.
//
// This program is primarily intended for illustrative purposes, and so
// there have been no attempts made to optimize either the reads or the
// byte swapping used to get the data in the correct endianness for the
// computer being used.
//
// Although the routines used for reading the image suggest that an RAXIS
// class should be developed, that is left as an exercise for the reader.
//
// This program assumes that the image was written using either SGI Control
// software or the instrument server provided with CrystalClear.  In these
// programs, the data is written in the file so that it is big-endian.
//
//
// Calling syntax:
//
// ReadRAXISImage.exe filename
//
// where filename is the name of the image file to be read.
//
// LEGAL ISSUES
// ============
//
// In plain English:
//
// 1. We don't promise that this software works.  (But if you find any bugs,
//    please let us know!)
//
// 2. You can use this software for whatever you want.  You don't have to
//    pay us.
//
// 3. You may not pretend that you wrote this software.  If you use it in a
//    program, you must acknowledge somewhere in your documentation that
//    you've used Rigaku code.
//
// In legalese:
//
// Rigaku makes NO WARRANTY or representation, either express or implied,
// with respect to this software, its quality, accuracy, merchantability, or
// fitness for a particular purpose.  This software is provided "AS IS", and
// you, its user, assume the entire risk as to its quality and accuracy.
//
// Permission is hereby granted to use, copy, modify, and distribute this
// software (or portions thereof) for any purpose, without fee, subject to
// these conditions:
//
// 1. If any part of the source code for this software is distributed, then
//    this notice must be included, with this copyright and no-warranty
//    notice unaltered; and any additions, deletions, or changes to the
//    original files must be clearly indicated in accompanying documentation.
//
// 2. If only executable code is distributed, then the accompanying
//    documentation must state that "this software is based in part on
//    software developed at Rigaku".
//
// 3. Permission for use of this software is granted only if the user accepts
//    full responsibility for any undesirable consequences; Rigaku accepts
//    NO LIABILITY for damages of any kind.
//

/****************************************************************************
 *                              Include Files                               *
 ****************************************************************************/

#ifdef WIN32
#pragma warning(disable:4786) // so warning about debug term > 255 chars ignored
#endif

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

/****************************************************************************
 *                               Definitions                                *
 ****************************************************************************/

/****************************************************************************
 *                                 Typedefs                                 *
 ****************************************************************************/

/****************************************************************************
 *                          Structure definitions                           *
 ****************************************************************************/

struct RigakuHeader {
  char  name[10];     /* type of instrument */
  char  vers[10];     /* version */
  char  crnm[20];     /* crystal name */
  char  crsy[12];     /* crystal system */
  float alng;         /* a */
  float blng;         /* b */
  float clng;         /* c */
  float alfa;         /* alpha */
  float beta;         /* beta */
  float gamm;         /* gamma */
  char  spgr[12];     /* space group symbol */
  float mosc;         /* mosaic spread */
  char  memo[80];     /* memo, comments */
  char  res1[84];     /* reserved space for future use */

  char  date[12];     /* date of measurement */
  char  opnm[20];     /* username of account collecting image */
  char  trgt[4];      /* type of X-ray target (Cu, Mo, etc.) */
  float wlng;         /* X-ray wavelength */
  char  mcro[20];     /* type of monochromator */
  float m2ta;         /* monochromator 2theta (deg) */
  char  colm[20];     /* collimator size and type */
  char  filt[4];      /* filter type (Ni, etc.) */
  float camr;         /* crystal-to-detector distance */
  float vltg;         /* generator voltage (kV) */
  float crnt;         /* generator current (mA) */
  char  focs[12];     /* focus info */
  char  optc[80];     /* xray memo */
  long  cyld;         /* IP shape, 0=flat,1=cylinder */
  float weis;         /* Weissenberg oscillation 1 */
  char  res2[56];     /* reserved space for future use */

  char  mnax[4];      /* crystal mount axis closest to spindle axis */
  char  bmax[4];      /* crystal mount axis closest to beam axis */
  float phi0;         /* datum phi angle (deg) */
  float phis;         /* phi oscillation start angle (deg) */
  float phie;         /* phi oscillation end angle (deg) */
  long  oscn;         /* frame number */
  float fext;         /* exposure time (min) */
  float drtx;         /* direct beam X position */
  float drtz;         /* direct beam Z position */
  float omga;         /* goniostat angle omega */
  float fkai;         /* goniostat angle chi */
  float thta;         /* goniostat angle 2theta */
  float mu;           /* spindle inclination angle */
  char  scan[204];    /* reserved space for future use */
  /* This space is now used for storing the scan
     template information - tlh, 01 Feb 1999 */

  long  xpxl;         /* number of pixels in X direction */
  long  zpxl;         /* number of pixels in Z direction */
  float xsiz;         /* size of pixel in X direction (mm) */
  float zsiz;         /* size of pixel in Z direction (mm) */
  long  rlng;         /* record length (bytes) */
  long  rnum;         /* number of records (lines) in image */
  long  ipst;         /* starting line number */
  long  ipnm;         /* IP number */
  float rato;         /* photomultiplier output hi/lo ratio */
  float ft_1;         /* fading time, end of exposure to start of read */
  float ft_2;         /* fading time, end of exposure to end of read */
  char  host[10];     /* type of computer (IRIS, VAX) => endian */
  char  ip[10];       /* type of IP */
  long  dr_x;         /* horizontal scanning code: 0=left->right, 1=>right->left */
  long  dr_z;         /* vertical scanning code: 0=down->up, 1=up->down */
  long  drxz;          /* front/back scanning code: 0=front, 1=back */
  float shft;         /* pixel shift, R-AXIS V */
  float ineo;         /* intensity ratio E/O R-AXIS V */
  long  majc;         /* magic number to indicate next values are legit */
  long  naxs;         /* Number of goniometer axes */
  float gvec[5][3];   /* Goniometer axis vectors */
  float gst[5];       /* Start angles for each of 5 axes */
  float gend[5];      /* End angles for each of 5 axes */
  float goff[5];      /* Offset values for each of 5 axes */
  long  saxs;         /* Which axis is the scan axis? */
  char  gnom[40];     /* Names of the axes (space or comma separated?) */

/*
 * Most of below is program dependent.  Different programs use
 * this part of the header for different things.  So it is essentially
 * a big "common block" area for dumping transient information.
 */
  char  file[16];     /* */
  char  cmnt[20];     /* */
  char  smpl[20];     /* */
  long  iext;         /* */
  long  reso;         /* */
  long  save;         /* */
  long  dint;         /* */
  long  byte;         /* */
  long  init;         /* */
  long  ipus;         /* */
  long  dexp;         /* */
  long  expn;         /* */
  long  posx[20];     /* */
  long  posy[20];     /* */
  int   xray;         /* */
  char  res5[768];    /* reserved space for future use */
};

struct axis {

  axis() : start(0.0),end(0.0),offset(0.0),name("") {}

  double start,end,offset;
  std::string name;
};

/****************************************************************************
 *                               Enumerations                               *
 ****************************************************************************/

/****************************************************************************
 *                                Constants                                 *
 ****************************************************************************/

const int MinimumRigakuHeaderLength = 1400;

std::string ERROR_OPENING_FILE   = "ERROR Can not open file: ";
std::string ERROR_READING_HEADER = "ERROR Can not read Rigaku-style header.";
std::string ERROR_READING_IMAGE  = "ERROR Can not read image data.";

/****************************************************************************
 *                           Function prototypes                            *
 ****************************************************************************/

// This routine returns a string containing a description of the program
// usage.
std::string ProgramUsage(void);

// This routine reads a Rigaku-style header from the input stream
// <VAR>InputFile</VAR>, and returns a boolean indicating if the header
// could be read successfully.  If the value returned is <I>true</I>,
// then the header was read successfully and the header information will
// be placed in the header structure <VAR>Header</VAR>.  If the value
// returned is <I>false</I>, then an error occurred while attempting to
// read the header and the contents of <VAR>Header</VAR> are undefined.
bool ReadHeader(std::ifstream &InputFile, RigakuHeader &Header);

// This routine reads a Rigaku-style image from the input stream
// <VAR>InputFile</VAR>, and returns a boolean indicating if the image
// could be read successfully.  The routine assumes that the header
// structure <VAR>Header</VAR> contains valid header information for
// the image being read.  If the value returned is <I>true</I>,
// then the image was read successfully and the image will be
// placed as (possibly packed) unsigned shorts at the memory location
// indicated by pImageBuffer.  If the value returned is <I>false</I>,
// then an error occurred while attempting to read the image and the
// contents of <VAR>pImageBuffer</VAR> are undefined.  It is up to the
// calling routine to ensure that <VAR>pImageBuffer</VAR> contains
// enough space to hold the image.
bool ReadImageBuffer(std::ifstream &InputFile,
                     const RigakuHeader &Header,
                     unsigned short * const pImageBuffer);

// This routine swaps the unsigned short <VAR>Value</VAR> so that the
// endianness of <VAR>Value</VAR> is reversed and then returns the
// swapped value.  This routine assumes that the size of a unsigned short
// is 2 bytes.
unsigned short Swap(unsigned short Value);

// This routine swaps the long <VAR>Value</VAR> so that the endianness
// of <VAR>Value</VAR> is reversed and then returns the swapped value.
// This routine assumes that the size of a long is 4 bytes.
long Swap(long Value);

// This routine swaps the float <VAR>Value</VAR> so that the endianness
// of <VAR>Value</VAR> is reversed and then returns the swapped value.
// This routine assumes that the size of a float is 4 bytes and that the
// IEEE standard format for floats is used.
float Swap(float Value);

// This routine swaps the values in a Rigaku-style header,
// <VAR>Header</VAR> so that the endianness of the <CODE>long</CODE>
// and </CODE>float</CODE> is reversed.
void Swap(RigakuHeader &Header);

// This routine returns a long containing the actual value of a packed
// R-AXIS image pixel, <VAR>PixelValue</VAR>.  The Rigaku-style header
// contains information concerning the packing scheme used.
long UnPackPixelValue(const RigakuHeader &Header,
                      const unsigned short PixelValue);

// This routine writes out pixel values for some prefined pixel positions
// in the image.  The values are written to standard output.  The header
// of the image is contained in <VAR>Header</VAR> and
// <VAR>pImageBuffer</VAR> points to the first pixel in the image.
void OutputSamplePixelValues(const RigakuHeader &Header,
                             unsigned short * const pImageBuffer);

// This routine writes out some basic information contained in the
// header <VAR>Header</VAR>.
void OutputHeaderInfo(const RigakuHeader &Header);

// This routine writes out the axis information, <VAR>Axis</VAR>, on
// the output stream <VAR>Stream</VAR>
std::ostream &operator<<(std::ostream &Stream, axis &Axis);

/****************************************************************************
 ****************************************************************************/
int main(int argc,
         char *argv[])
{
  // Make sure that there is a filename on the command line.

  if(2 != argc){
    std::cerr << ProgramUsage() << std::endl;
    return 1;
  }

  // Open the file

  std::ifstream inp(argv[1],std::ios::in|std::ios::binary);
  if(inp.fail() || inp.bad()){
    std::cerr << ERROR_OPENING_FILE << argv[1] << std::endl;
    if(inp.is_open())
      inp.close();
    return 2;
  }

  // Read the header of the image

  RigakuHeader header;
  if(!ReadHeader(inp,header)){
    std::cerr << ERROR_READING_HEADER << std::endl;
    if(inp.is_open())
      inp.close();
    return 3;
  }

  OutputHeaderInfo(header);

  // Allocate space to store the entire image

  long size = header.rlng*header.rnum;
  unsigned short *pBuffer = new unsigned short [size/2];

  // Read in the image

  if(!ReadImageBuffer(inp,header,pBuffer)){
    std::cerr << ERROR_READING_IMAGE << std::endl;
    if(inp.is_open())
      inp.close();
    delete [] pBuffer;
    return 4;
  }

  if(inp.is_open())
    inp.close();

  OutputSamplePixelValues(header,pBuffer);

  delete [] pBuffer;

  return 0;
}

/****************************************************************************
 ****************************************************************************/
std::string ProgramUsage(void)
{
  std::string s = "ReadRAXISImage.exe image_filename\n\n";
  s += "image_filename : Filename of image to be read.\n";

  return s;
}

/****************************************************************************
 ****************************************************************************/
bool ReadHeader(std::ifstream &InputFile,
                RigakuHeader &Header)
{
  // Read in the minimal amount of header information.

  InputFile.read(reinterpret_cast<char *>(&Header),MinimumRigakuHeaderLength);

  if(InputFile.bad() || InputFile.fail())
    return false;

  // Check to see if endianness of the header matches endianness of this
  // architecture.  If it doesn't, then we have to swap the non-character
  // elements of the header.

  unsigned short one = 1;
  char *p = reinterpret_cast<char *>(&one);
  if(1 == static_cast<int>(*p))
    Swap(Header);

  // Read the additional padding in the header

  long size = Header.rlng-MinimumRigakuHeaderLength;
  p = new char [size];
  InputFile.read(p,size);
  delete [] p;

  return !(InputFile.bad() || InputFile.fail());
}

/****************************************************************************
 ****************************************************************************/
bool ReadImageBuffer(std::ifstream &InputFile,
                     const RigakuHeader &Header,
                     unsigned short * const pImageBuffer)
{
  // Figure out how much of the image is contained in the file.  Note that
  // if the image file contains only a partial image, we will NOT construct
  // an entire image.  The header element ipst should be used to determine
  // the line in the entire image corresponding to a line in the image file.
  //
  // Note that rlng is the size of a record in bytes, and we are assuming
  // 2 bytes/pixel.

  long size = Header.rlng*Header.rnum/2;

  // Read in the data.

  InputFile.read(reinterpret_cast<char *>(pImageBuffer),size);

  if(InputFile.bad() || InputFile.fail())
    return false;

  // Do we need to swap the data?

  unsigned short one = 1;
  char *p = reinterpret_cast<char *>(&one);

  if(1 == static_cast<int>(*p)){

    // Swap each pixel value in the image buffer

    unsigned short *p = const_cast<unsigned short *>(pImageBuffer);
    unsigned short *pEnd = p+size;
    for(NULL; p < pEnd; p++)
      *p = Swap(*p);

  }

  return true;

}

/****************************************************************************
 ****************************************************************************/
unsigned short Swap(const unsigned short Value)
{
  unsigned short swap = Value;
  char *p = reinterpret_cast<char *>(&swap);

  char tmp = p[0];
  p[0] = p[1];
  p[1] = tmp;

  return swap;
}

/****************************************************************************
 ****************************************************************************/
long Swap(const long Value)
{
  long swap = Value;
  char *p = reinterpret_cast<char *>(&swap);

  char tmp = p[0];
  p[0] = p[3];
  p[3] = tmp;

  tmp = p[1];
  p[1] = p[2];
  p[2] = tmp;

  return swap;
}

/****************************************************************************
 ****************************************************************************/
float Swap(const float Value)
{
  float swap = Value;
  char *p = reinterpret_cast<char *>(&swap);

  char tmp = p[0];
  p[0] = p[3];
  p[3] = tmp;

  tmp = p[1];
  p[1] = p[2];
  p[2] = tmp;

  return swap;
}

/****************************************************************************
 ****************************************************************************/
void Swap(RigakuHeader &Header)
{

// Swap the bytes in the long elements of the header.

  Header.cyld = Swap(Header.cyld);
  Header.oscn = Swap(Header.oscn);
  Header.xpxl = Swap(Header.xpxl);
  Header.zpxl = Swap(Header.zpxl);
  Header.rlng = Swap(Header.rlng);
  Header.rnum = Swap(Header.rnum);
  Header.ipst = Swap(Header.ipst);
  Header.ipnm = Swap(Header.ipnm);
  Header.dr_x = Swap(Header.dr_x);
  Header.dr_z = Swap(Header.dr_z);
  Header.drxz = Swap(Header.drxz);

// Swap the bytes in the float elements of the header.

  Header.alng = Swap(Header.alng);
  Header.blng = Swap(Header.blng);
  Header.clng = Swap(Header.clng);
  Header.alfa = Swap(Header.alfa);
  Header.beta = Swap(Header.beta);
  Header.gamm = Swap(Header.gamm);
  Header.mosc = Swap(Header.mosc);
  Header.wlng = Swap(Header.wlng);
  Header.m2ta = Swap(Header.m2ta);
  Header.camr = Swap(Header.camr);
  Header.vltg = Swap(Header.vltg);
  Header.crnt = Swap(Header.crnt);
  Header.weis = Swap(Header.weis);
  Header.phi0 = Swap(Header.phi0);
  Header.phis = Swap(Header.phis);
  Header.phie = Swap(Header.phie);
  Header.fext = Swap(Header.fext);
  Header.drtx = Swap(Header.drtx);
  Header.drtz = Swap(Header.drtz);
  Header.omga = Swap(Header.omga);
  Header.fkai = Swap(Header.fkai);
  Header.thta = Swap(Header.thta);
  Header.mu   = Swap(Header.mu);
  Header.xsiz = Swap(Header.xsiz);
  Header.zsiz = Swap(Header.zsiz);
  Header.rato = Swap(Header.rato);
  Header.ft_1 = Swap(Header.ft_1);
  Header.ft_2 = Swap(Header.ft_2);

// New variables for RAPID image headers

  Header.shft = Swap(Header.shft);
  Header.ineo = Swap(Header.ineo);
  Header.majc = Swap(Header.majc);
  Header.naxs = Swap(Header.naxs);
  for(int i = 0; i < 5; i++){
    for(int j = 0; j < 3; j++)
      Header.gvec[i][j] = Swap(Header.gvec[i][j]);
    Header.gst[i] = Swap(Header.gst[i]);
    Header.gend[i] = Swap(Header.gend[i]);
    Header.goff[i] = Swap(Header.goff[i]);
  }
  Header.saxs = Swap(Header.saxs);

  return;
}

/****************************************************************************
 ****************************************************************************/
long UnPackPixelValue(const RigakuHeader &Header,
                      const unsigned short PixelValue)
{
  // Pixels are packed so that if the most significant bit is set, then
  // the value has been bit-shifted by the amount specified in the rato
  // header element.
  //
  // Note that this form of packing is actually quite inspired as the
  // amount of data lost due to packing is much less than the error in
  // the measurement due to normal counting statistics (assuming that
  // they are Poissonian).

  if(0x7fff >= PixelValue)
    return PixelValue;

  // We only need to use the lower 15 bits of the value to uncompress it.

  long Value = 0x7fff&PixelValue;

  // Shift the value by the ratio in the header.

  Value *= Header.rato;

  return Value;
}

/****************************************************************************
 ****************************************************************************/
void OutputSamplePixelValues(const RigakuHeader &Header,
                             unsigned short * const pImageBuffer)
{
  // Output a few pixel values.  Note that the first line in the image
  // buffer may not correspond to the first line in a "complete" image.

  std::cout << std::endl << "Sample pixel information: " << std::endl;

  std::cout << "   First pixel (0," << Header.ipst << ") in image = " <<
            UnPackPixelValue(Header,*pImageBuffer) << std::endl;

  long index = (Header.rnum/2-1)*Header.xpxl+(Header.xpxl/2-1);
  std::cout << "   Middle pixel (" << Header.xpxl/2 << ',' <<
            Header.rnum/2+Header.ipst << ") in image = " <<
            UnPackPixelValue(Header,pImageBuffer[index]) << std::endl;

  index = Header.rnum*Header.xpxl-1;
  std::cout << "   Last pixel (" << Header.xpxl << ',' <<
            Header.rnum+Header.ipst << ") in image = " <<
            UnPackPixelValue(Header,pImageBuffer[index]) << std::endl;

  return;
}

/****************************************************************************
 ****************************************************************************/
std::string fromFORTRAN(char *string, int size)
{
  std::ostringstream s;
  for(int i = 0; i < size; i++)
    s << string[i];
  return s.str();
}

/****************************************************************************
 ****************************************************************************/
void OutputHeaderInfo(const RigakuHeader &Header)
{
  std::cout << "Header information: " << std::endl;

  std::cout << "   InstrumentType: " << fromFORTRAN((char *)Header.name,10) << std::endl;
  std::cout << "   Version: " << fromFORTRAN((char *)Header.vers,10) << std::endl;
  std::cout << "   Crystal Name: " << fromFORTRAN((char *)Header.crnm,20) << std::endl;
  std::cout << "   Crystal System: " << fromFORTRAN((char *)Header.crsy,12) << std::endl;
  std::cout << "   a: " << Header.alng << std::endl;
  std::cout << "   b: " << Header.blng << std::endl;
  std::cout << "   c: " << Header.clng << std::endl;
  std::cout << "   alpha: " << Header.alfa << std::endl;
  std::cout << "   beta: " << Header.beta << std::endl;
  std::cout << "   gamma: " << Header.gamm << std::endl;
  std::cout << "   Space Group: " << fromFORTRAN((char *)Header.spgr,12) << std::endl;
  std::cout << "   Mosaicity: " << Header.mosc << std::endl;
  std::cout << "   Memo: " << fromFORTRAN((char *)Header.memo,80) << std::endl;

  std::cout << "   Date: " << fromFORTRAN((char *)Header.date,12) << std::endl;
  std::cout << "   User: " << fromFORTRAN((char *)Header.opnm,20) << std::endl;
  std::cout << "   Xray Target: " << fromFORTRAN((char *)Header.trgt,4) << std::endl;
  std::cout << "   Wavelength: " << Header.wlng << std::endl;
  std::cout << "   Monochromator: " << fromFORTRAN((char *)Header.mcro,20) << std::endl;
  std::cout << "   Monochromator 2theta: " << Header.m2ta << std::endl;
  std::cout << "   Collimator: " << fromFORTRAN((char *)Header.colm,20) << std::endl;
  std::cout << "   Filter: " << fromFORTRAN((char *)Header.filt,4) << std::endl;
  std::cout << "   Generator Voltage: " << Header.vltg << std::endl;
  std::cout << "   Generator Current: " << Header.crnt << std::endl;
  std::cout << "   Focus: " << fromFORTRAN((char *)Header.focs,12) << std::endl;
  std::cout << "   Xray Memo: " << fromFORTRAN((char *)Header.optc,80) << std::endl;
  std::cout << "   IP shape: " << ((1 == Header.cyld) ? "cylindrical" :
                                   "flat") << std::endl;
  std::cout << "   OscillationType: " << ((1 == Header.weis) ? "weissenberg" : "regular") << std::endl;

  std::cout << "   Crystal Mount (spindle axis): " << fromFORTRAN((char *)Header.mnax,4) << std::endl;
  std::cout << "   Crystal Mount (beam axis): " << fromFORTRAN((char *)Header.bmax,4) << std::endl;
  std::cout << "   Frame Number: " << Header.oscn << std::endl;
  std::cout << "   Exposure time: " << Header.fext << std::endl;
  std::cout << "   Direct beam position: (" << Header.drtx << ',' <<
            Header.drtz << ')' << std::endl;
  std::cout << "   Spindle Inclication Angle: " << Header.mu << std::endl;
  std::cout << "   Image Template: " << fromFORTRAN((char *)Header.scan,204) << std::endl;

  std::cout << "   Pixel information" << std::endl <<
            "      size: (" << Header.xsiz << ',' << Header.zsiz << ')' <<
            std::endl <<
            "      number: (" << Header.xpxl << ',' << Header.zpxl << ')' <<
            std::endl;
  std::cout << "   Line information" << std::endl <<
            "      record length: " << Header.rlng << std::endl <<
            "      total: " << Header.rnum << std::endl <<
            "      starting: " << Header.ipst << std::endl;
  std::cout << "   IP Number: " << Header.ipnm << std::endl;
  std::cout << "   Photomultiplier Ratio: " << Header.rato << std::endl;
  std::cout << "   Fade Time (to start of read): " << Header.ft_1 << std::endl;
  std::cout << "   Fade Time (to end of read): " << Header.ft_2 << std::endl;
  std::cout << "   Host Type/Endian: " << fromFORTRAN((char *)Header.host,10) << std::endl;
  std::cout << "   IP Type: " << fromFORTRAN((char *)Header.ip,10) << std::endl;
  std::cout << "   Horizontal Scan: " << ((0 == Header.dr_x) ? "left to right" : "right to left") << std::endl;
  std::cout << "   Vertical Scan: " << ((0 == Header.dr_z) ? "down to up" : "up to down") << std::endl;
  std::cout << "   Front/Back Scan: " << ((0 == Header.drxz) ? "front" : "back") << std::endl;
  std::cout << "   Pixel Shift (RAXIS V): " << Header.shft << std::endl;
  std::cout << "   Even/Odd Intensity Ratio (RAXIS V): " << Header.ineo << std::endl;


  axis Axes[5];

  std::string s = Header.name;
  bool IsRapid = std::string::npos != s.find("R-AXIS-CS");
  std::string ScanAxis = "Phi";

  if(0 != Header.majc){
    std::istringstream is(Header.gnom);
    for(int i = 0; i < Header.naxs; i++){
      is >> Axes[i].name;
      Axes[i].start = Header.gst[i];
      Axes[i].end = Header.gend[i];
      Axes[i].offset = Header.goff[i];
    }
    ScanAxis = Axes[Header.saxs].name;
  }
  else{
    Axes[0].name = "Omega";
    Axes[2].name = "Phi";
    if(!IsRapid){
      Axes[0].start = Header.omga;
      Axes[0].end = Header.omga;
      Axes[0].offset = 0.0;
      Axes[2].start = Header.phis;
      Axes[2].end = Header.phie;
      Axes[2].offset = Header.phi0;
      ScanAxis = Axes[2].name;
    }
    else{
      Axes[2].start = Header.omga;
      Axes[2].end = Header.omga;
      Axes[2].offset = 0.0;
      Axes[0].start = Header.phis;
      Axes[0].end = Header.phie;
      Axes[0].offset = Header.phi0;
      ScanAxis = Axes[0].name;
    }
    Axes[1].name = "Chi";
    Axes[1].start = Header.fkai;
    Axes[1].end = Header.fkai;
    Axes[1].offset = 0.0;
  }

  std::cout << "   Instrument axes (using " <<
            ((0 != Header.majc) ? "extended" : "conventional") <<
            " section of header):" << std::endl;
  std::cout << "      Scan axis: " << ScanAxis << std::endl;
  for(int i = 0; i < 5 && "" != Axes[i].name; i++)
    std::cout << "      " << Axes[i] << std::endl;

  axis theta,dist;
  theta.name = "2Theta";
  theta.start = Header.thta;
  theta.end = Header.thta;
  theta.offset = 0.0;
  dist.name = "Distance";
  dist.start = Header.camr;
  dist.end = Header.camr;
  dist.offset = 0.0;
  std::cout << "      " << theta << std::endl;
  std::cout << "      " << dist << std::endl;

  return;
}

/****************************************************************************
 ****************************************************************************/
std::ostream &operator<<(std::ostream &a,
                         axis &b)
{
  a << b.name << ": start = " << b.start << ", end = " << b.end <<
    ", offset = " << b.offset;
  return a;
}
