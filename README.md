# mini_FFT

### Fork from 'fix_fft' to reduce code size and needed RAM

---

Alike fix_fft, this works with int8 to compute faster, here the only change are:
  - reduced size of sine wave to 1/4,
  - optimized values in the sine wave array,
  - wiped the fft reverse,
  - imaginaries values are created by the mini_FFT routine and start with zero value,
  - changed process for the values:
  
      - INPUT, only send sampled signal (real parts array),
      - temporary using an array for the imaginaries (same size than the real part array)
      - OUTPUT, both reals and imaginaries are recorded into the real array
      
            r = [0 to n/2[    and i = [n/2 to n[
            
