  /*
    mini_FFT() - perform forward Fast Fourier Transform.
      fr[n] are real array, filled with n integers from [-127, 127]
      with m [1, 8], this will allows   n = [0, 2**m[ (ie with m = 2, n from [0, 3])
    result in fr[n]
      reals       from fr[0]   to fr[n/2 - 1]
      imaginaries from fr[n/2] to fr[n - 1]
      each value with an amplitude = 127, so:
          cos(i) = fr|i] / 127
          sin(i) = fi[i] / 127
      
    void mini_FFT(char fr[], int m);
  */



#include  <avr/pgmspace.h>
#include  <mini_FFT_v1.h>



PROGMEM const unsigned int  FFTe = 7;          // FFT factor to perform amount of frequencies bins (7 => FFTe = 128)
PROGMEM const unsigned int  FFTn = bit(FFTe);  // 2**FFTe, how many bins we want, final is FFTn / 2



char          FFT_Sample[FFTn];                // FFT bins array



void setup() {
  sampleSound();
  
  /****************************************************
    Performes FFT

    Since we filled FFT_Sample[] with sound samples,
    mini_FFT will compute and provide the bins.

    Return result into:
    - fReal[0, n/2 -1] for real parts
    - fReal[n/2, n -1] for imaginary parts
  *****************************************************/
  mini_FFT(FFT_Sample, FFTe);
}

void loop() {
  ;
}

void sampleSound() {
  // create a rountine here to sample and to fill up FFT_Sample[]
}
