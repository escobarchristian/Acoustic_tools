""" spectogram calculations for measured data
"""
import numpy as np

def spectrogram_PSD_matlab(x,window,noverlap,nfft,fs):
    """
    This function uses the same Matlab implementation to return the power spectral density from the time series x
    Inputs:
        x (np.ndarray): Time series
        fs (float): sampling frequency, Hz
        nfft (int): number of samples for each FFT
            frequency bin size, df = nfft/fs
        window (np.ndarray): Time window. Use window to divide the signal into segments.
        noverlap (int) : number of samples to overlap in time
    Outputs:
    
        P : ndarray of float
            Power spectral density in linear units
        t : ndarray of float
            Time array corresponding to the time steps
        f : ndarray of float
            Frequency array            
    """
    x = np.squeeze(x.astype(np.float32)) # Format the input time series
    W = np.mean(np.multiply(window, window)) # normalizing factor
    nx=len(window) # Window size
    time_step=nx/fs # Time spacing
    total_time=len(x)/fs
    # determine how many unique times will be in the final spectrogram
    perc_noverlap=noverlap/nx # Overlap in percentage 0-99
    time_step=time_step*(1-perc_noverlap) # Time spacing
    nt = np.floor(total_time/time_step)  # Number of time points
    # make time array of final spectrogram
    t = np.arange(0,nt) * time_step
    # make nt an integer
    nt = int(nt.item())
    f = np.fft.fftfreq(nfft, d=1/fs)[:int(np.ceil(nfft/2)+1)]
    # frequency bin size
    df = f[2] - f[1]
    # number of frequencies returned from the FFT
    nf = len(f)

    P = np.zeros([nt,nf])
    nstart=0   

    for time_idx in np.arange(nt):
        # beginning at index nstart sample, do the FFT on a signal ns samples in length
        nstop = nstart + nx
        x1 = x[int(nstart):int(nstop)]

        if len(x1)==nx:
            # apply the window
            x1 = np.multiply(x1.T, window)
            if nx > nfft: # If the number of time steps is greater than nfft then do a cyclic average
                n_aux=int(np.ceil(nx/nfft)) 
                new_size=n_aux*nfft
                x1=np.resize(x1,new_size)
                if new_size>nx:
                    x1[nx:]=0
                x1=x1.reshape(n_aux,nfft)
                x1=np.sum(x1,0)                
            # calculate the fft and obtain the single sided part
            Y = np.fft.fft(x1, nfft)[:int(nfft/2)+1]
            # Scale the signal based on the normalizing factor 
            Scale = (2/nfft/fs/W)/(nx/nfft)
            Y = Scale * np.multiply(np.conj(Y), Y)    
            P[time_idx] = np.real(Y)
            # update nstart to the next time step
            nstart = nstart + nx - noverlap
        else:
            P=P[0:time_idx,:]
            t=t[0:time_idx]
            break
    P[:,0]=P[:,0]/2 # Unscale the factor of two.  % Don't double unique Nyquist point
    P[:,nf-1]=P[:,nf-1]/2 # Unscale the factor of two.  % Don't double unique Nyquist point
    f=np.abs(f)
    return P,f,t
