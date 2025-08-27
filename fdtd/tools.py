import numpy as np
import matplotlib.pyplot as plt


def fft_plot(signal: list, time_diff: float, output_file="fft.png"):
    
  N = len(signal)
  fft_signal = np.fft.fft(signal)
  fft_signal = np.abs(fft_signal) / N
  
  freqs = np.fft.fftfreq(N, time_diff)
  freqs_ghz = freqs / 1e9
  
  mask = freqs_ghz >= 0
  freqs_ghz_pos = freqs_ghz[mask]
  fft_signal_pos = fft_signal[mask]
  
  # Greatest magnitude
  idx_max = np.argmax(fft_signal_pos)
  freq_max = freqs_ghz_pos[idx_max]
  mag_max = fft_signal_pos[idx_max]
  
  plt.figure(figsize=(8,4))
  plt.plot(freqs_ghz_pos, fft_signal_pos)
  plt.xlabel('Frequência (GHz)')
  plt.ylabel('Magnitude FFT')
  plt.title('Espectro de Frequência')
  
  # Marks the greatest magnitude
  plt.plot(freq_max, mag_max, 'ro')
  plt.annotate(f'Max: {freq_max:.3f} GHz',
               xy=(freq_max, mag_max),
               xytext=(freq_max, mag_max*1.1),
               arrowprops=dict(facecolor='red', shrink=0.05),
               ha='center')
  
  plt.grid(True)
  plt.savefig(output_file, dpi=300)
