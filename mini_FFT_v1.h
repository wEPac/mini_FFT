#ifndef   MINIFFT_H
  #define   MINIFFT_H
  #include  <avr/pgmspace.h>

  #if ARDUINO >= 100
    #include  "Arduino.h"
  #else
    #include  "WProgram.h"
  #endif

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
  */
  void mini_FFT(char fr[], int m);

#endif
