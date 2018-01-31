mini_FFT
Fork from 'fix_fft' and 'fast_fft' to reduce code size and needed RAM 

Alike fix_fft, this works with int8 to compute faster, here the only change are:
  - reduced size of sine wave to 1/4,
  - wiped the fft reverse,
  - imaginaries values are created by the fft routine and start with zero value,
  - changed process to record values, using temporary array for the imaginaries and returning imaginaries result into the last half array part of real
