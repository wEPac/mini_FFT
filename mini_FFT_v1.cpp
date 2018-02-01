#include "mini_FFT_v1.h"


/* based on the following... */
/* fix_fft.c - Fixed-point in-place Fast Fourier Transform  */
/*
  All data are fixed-point short integers, in which -32768
  to +32768 represent -1.0 to +1.0 respectively. Integer
  arithmetic is used for speed, instead of the more natural
  floating-point.

  For the forward FFT (time -> freq), fixed scaling is
  performed to prevent arithmetic overflow, and to map a 0dB
  sine/cosine wave (i.e. amplitude = 32767) to two -6dB freq
  coefficients. The return value is always 0.

  For the inverse FFT (freq -> time), fixed scaling cannot be
  done, as two 0dB coefficients would sum to a peak amplitude
  of 64K, overflowing the 32k range of the fixed-point integers.
  Thus, the fix_fft() routine performs variable scaling, and
  returns a value which is the number of bits LEFT by which
  the output must be shifted to get the actual amplitude
  (i.e. if fix_fft() returns 3, each value of fr[] and fi[]
  must be multiplied by 8 (2**3) for proper scaling.
  Clearly, this cannot be done within fixed-point short
  integers. In practice, if the result is to be used as a
  filter, the scale_shift can usually be ignored, as the
  result will be approximately correctly normalized as is.

  Written by:  Tom Roberts  11/8/89
  Made portable:  Malcolm Slaney 12/15/94 malcolm@interval.com
  Enhanced:  Dimitrios P. Bouras  14 Jun 2006 dbouras@ieee.org
  Modified for 8bit values David Keller  10.10.2010
  Mofified: 20 Jan 2018 Eric Paquot:
  - reduced size of sine wave to 1/4,
  - optimized values in the sine wave array,
  - wiped the fft reverse,
  - imaginaries values are created by the fft routine and start
  with zero value,
  - changed process to record values, using temporary array for
  the imaginaries and returning imaginaries result into the last
  half array part of real
*/


#define N_WAVE      256    // full length of Sinewave[]
#define LOG2_N_WAVE 8      // log2(N_WAVE), to know how many bits from N_WAVE



/*
  Since we only use 3/4 of N_WAVE, we define only
  this many samples, in order to conserve data space.
  
  From the symetric curve of the wave, we only need
  to keep 1/4 of N_WAVE
*/
//const int8_t Sinewave[N_WAVE - N_WAVE / 4] PROGMEM = {
const int8_t Sinewave[N_WAVE / 4] PROGMEM = {
    0,   2,   3,   5,   6,   8,   9,  11,
   12,  14,  15,  17,  18,  20,  21,  23,
   24,  25,  27,  28,  30,  31,  32,  34,
   35,  36,  38,  39,  40,  41,  42,  44,
   45,  46,  47,  48,  49,  50,  51,  52,
   53,  54,  54,  55,  56,  57,  57,  58,
   59,  59,  60,  60,  61,  61,  62,  62,
   62,  63,  63,  63,  63,  63,  63,  63/*,
  // ^ already divided by 2 to reduce calculs
    0,   3,   6,   9,  12,  15,  18,  21,
   24,  28,  31,  34,  37,  40,  43,  46,
   48,  51,  54,  57,  60,  63,  65,  68,
   71,  73,  76,  78,  81,  83,  85,  88,
   90,  92,  94,  96,  98, 100, 102, 104,
  106, 108, 109, 111, 112, 114, 115, 117,
  118, 119, 120, 121, 122, 123, 124, 124,
  125, 126, 126, 127, 127, 127, 127, 127/*,

  127, 127, 127, 127, 127, 127, 126, 126,
  125, 124, 124, 123, 122, 121, 120, 119,
  118, 117, 115, 114, 112, 111, 109, 108,
  106, 104, 102, 100,  98,  96,  94,  92,
   90,  88,  85,  83,  81,  78,  76,  73,
   71,  68,  65,  63,  60,  57,  54,  51,
   48,  46,  43,  40,  37,  34,  31,  28,
   24,  21,  18,  15,  12,   9,   6,   3,

     0,   -3,   -6,   -9,  -12,  -15,  -18,  -21,
   -24,  -28,  -31,  -34,  -37,  -40,  -43,  -46,
   -48,  -51,  -54,  -57,  -60,  -63,  -65,  -68,
   -71,  -73,  -76,  -78,  -81,  -83,  -85,  -88,
   -90,  -92,  -94,  -96,  -98, -100, -102, -104,
  -106, -108, -109, -111, -112, -114, -115, -117,
  -118, -119, -120, -121, -122, -123, -124, -124,
  -125, -126, -126, -127, -127, -127, -127, -127,

  -127, -127, -127, -127, -127, -127, -126, -126,
  -125, -124, -124, -123, -122, -121, -120, -119,
  -118, -117, -115, -114, -112, -111, -109, -108,
  -106, -104, -102, -100,  -98,  -96,  -94,  -92,
   -90,  -88,  -85,  -83,  -81,  -78,  -76,  -73,
   -71,  -68,  -65,  -63,  -60,  -57,  -54,  -51,
   -48,  -46,  -43,  -40,  -37,  -34,  -31,  -28,
   -24,  -21,  -18,  -15,  -12,   -9,   -6,   -3*/
};



/*
  FIX_MPY() - fixed-point multiplication & scaling.
  Substitute inline assembly for hardware-specific
  optimization suited to a particluar DSP processor.
  Scaling ensures that result remains 16-bit.
*/
inline char FIX_MPY(char a, char b)
{
  int c = ((int)a * (int)b) >> 6; // shift right one less bit (i.e. 15-1)
  b     = c & 0x01;               // last bit shifted out = rounding-bit
  a     = (c >> 1) + b;           // last shift + rounding bit

  return a;
}

/*
  fast_FFT() - perform forward/inverse fast Fourier transform.
  fr[n],fi[n] are real and imaginary arrays, both INPUT AND
  RESULT (in-place FFT), with 0 <= n < 2**m; set inverse to
  0 for forward transform (FFT), or 1 for iFFT.
*/
void mini_FFT(char fr[], int m)
{
  int  mr, nn, i, j, k, l, istep, n_fft;
  char qr, qi, tr, ti, wr, wi;

  n_fft = 1 << m;
  if (n_fft > N_WAVE) return; // max FFT size = N_WAVE

  // create imaginary parts then fill it up the with zero
  char fi[n_fft];
  memset(fi, 0, n_fft);
  
  mr    = 0;
  nn    = n_fft - 1;

  /* re-order data
     ie for FFT_N = 4, permutations are:
         1 <=>  8 :  2 <=>  4 :  3 <=> 12 :
         5 <=> 10 :           :  7 <=> 14 :
                  :           : 11 <=> 13 : 
     ie for FFT_N = 6, permutations are:
         1 <=> 32 :  2 <=> 16 :  3 <=> 48 :  4 <=>  8 :  5 <=> 40 :  6 <=> 24 :
         7 <=> 56 :           :  9 <=> 36 : 10 <=> 20 : 11 <=> 52 : 
        13 <=> 44 : 14 <=> 28 : 15 <=> 60 :           : 17 <=> 34 :            
        19 <=> 50 :           : 21 <=> 42 : 22 <=> 26 : 23 <=> 58 :
        25 <=> 38 :           : 27 <=> 54 :           : 29 <=> 46 :             
        31 <=> 62 :           :           :           : 35 <=> 49 : 
        37 <=> 41 :           : 39 <=> 57 :           :           :
                  : 43 <=> 53 :           :           :           : 47 <=> 61 :
                  :           :
                  : 55 <=> 59 : 
  */
  for (m = 1; m <= nn; ++m) {
    l = n_fft;
    do { l >>= 1; } while (mr + l > nn);
    
    mr = (mr & (l - 1)) + l; 
    if (mr <= m) continue;
    
    tr     = fr[m];
    fr[m]  = fr[mr];
    fr[mr] = tr;
    /* removed, cause we only need zero imaginaries at begining
    ti     = fi[m];
    fi[m]  = fi[mr];
    fi[mr] = ti;//*/
  }

  l = 1;
  k = LOG2_N_WAVE - 1;
  while (l < n_fft) {
    /*
      fixed scaling, for proper normalization --
      there will be log2(n) passes, so this results
      in an overall factor of 1/n, distributed to
      maximize arithmetic accuracy.
    */

    istep = l << 1;
    for (m = 0; m < l; ++m) {
      j = m << k;       // j = [0, N_WAVE / 2[
      /* change here to use only 1/4 from sine wave
      wr = ( pgm_read_byte_near(Sinewave + j + N_WAVE / 4)) >> 1;
      wi = (-pgm_read_byte_near(Sinewave + j)) >> 1;//*/
      if (j < N_WAVE / 4) {
        wr    =  pgm_read_byte_near(Sinewave + N_WAVE / 4 - j);
        wi    = -pgm_read_byte_near(Sinewave + j);
      }
      else {
        wr    = -pgm_read_byte_near(Sinewave + j - N_WAVE / 4);
        wi    = -pgm_read_byte_near(Sinewave + N_WAVE / 2 - j);
      }
      
      for (i = m; i < n_fft; i += istep) {
        j     = i + l;
        tr    = FIX_MPY(wr, fr[j]) - FIX_MPY(wi, fi[j]);
        ti    = FIX_MPY(wr, fi[j]) + FIX_MPY(wi, fr[j]);
        qr    = fr[i];
        qi    = fi[i];
        
        qr  >>= 1;
        qi  >>= 1;
        
        fr[j] = qr - tr;
        fi[j] = qi - ti;
        fr[i] = qr + tr;
        fi[i] = qi + ti;
      }
    }
    --k;
    l = istep;
  }

  // copy 1st half imaginaries parts (fi) into last half real parts (fr)
  mr = n_fft / 2;
  memcpy(fr + mr, fi, mr);
}


