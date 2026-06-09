function [omega,DFT] = FourierTransform(Signal,Time)
%%Perform FFT Given a signal and a time vector (these should have same
%%length)

    dft = fft(Signal);
    Fs  = 1/mean(diff(Time));
    L   = length(Time);
    
    P2  = abs(dft/L);
    DFT = fftshift(P2);
    
    if mod(L,2)==1
        omega = 2*pi*Fs.*linspace(-L/2,L/2,L)./L;
    end
    if mod(L,2)==0
        omega = 2*pi*Fs*(-L/2:L/2-1)./L;
    end

end
