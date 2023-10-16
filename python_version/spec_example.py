'''
Example of how to use tools.spectogram_cal function to create spectrograms using the same Matlab implementation
from a time waveform, x, i.e., 'spectrogram_PSD_matlab'
'''
import os
import numpy as np
from scipy.signal import chirp
from scipy.io import savemat
import matplotlib.pyplot as plt
import tools.spectogram_gen as spectrogram
import scipy.io.wavfile as wav

def main():
    read_audio=False # True to read in a wav file or False to generate a sample signal
    if read_audio:
        fs, x = wav.read("chirp_example.wav")
        t = np.arange(0,len(x))/fs
    else:
        fs = 8000 # in Hz, sampling frequency 
        duration=10 # Duration in seconds 
        N = fs * duration # Number of samples
        t = np.arange(0,N+1)/fs
        
        ## Example 
        # Generate a linear frequency modulated signal from 100 to 2000 Hz with a duration of 50 seconds
        x = chirp(t, f0=100, f1=2000, t1=10, method='linear') 
        noise = np.random.normal(0, 0.2, x.shape) 
        x += noise # Add some random noise
        
        ### Add a tone
        tone_freq=3500.0 # in Hz
        x2 = np.sin(2 * np.pi * tone_freq*t)
        x +=x2
        x = x/(np.max(np.abs(x))*1.001)
    # Calculate spectrogram
    nfft = 2**12 # number of samples for each FFT (2^12=4096). df=fs/nfft 
    tWindow = 0.1 #  Integration time in seconds, ideally it is nfft/fs
    NWindow = fs*tWindow # Window size, i.e., number of points to be considered in the fft
    window = np.hamming(NWindow) # Apply the window you would like to use.
    overlap=0.5 # Overlap (0.0 - 0.99).Typical value: 0.5
    # The final time step will be: tWindow*(1-overlap) in [sec]
    # The frequency step is: fs/nfft
    NOverlap=NWindow*overlap # Number of overlap points
    [P, F, T] = spectrogram.spectrogram_PSD_matlab(x,window,NOverlap,nfft,fs)
    spec = 10*np.log10(np.abs(P)) # In dB
    file_root = "chirp" 
    folder_name='Out_spectrograms/'
    # Create output directory
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    # Plot time series
    plt.figure()
    plt.plot(t,x)
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title('Time Series')
    plt.savefig(folder_name+file_root+'_time_series.png')
    
    # Plot spectrogram
    plt.figure()
    plt.pcolormesh( F, T,  spec, shading="auto")
    plt.clim(-100,-20) ######## Make sure you modify these values of dynamic range for the plot based on your intensity levels. 
    cbar = plt.colorbar()
    cbar.set_label('Power Spectral Density (dB)')
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Time (sec)")
    plt.title('Spectrogram')
    plt.savefig(folder_name+file_root+'_spec.png')            

    # Example for saving a dictionary as a .mat file using scipy
    mdic={'spec_py':spec,'f_py':F,'t_py':T,'overlap':NOverlap,'x_py':x,'NFFT_py':nfft,'fs':fs}
    savemat('Test_spectrogram.mat',mdic)    
    # Save audiofile
    filename = "chirp_example.wav"
    wav.write(filename, fs, x)
    
    plt.show()
if __name__ == "__main__":
    main()
    
    