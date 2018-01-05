function y = filt_common(z, Fs, Fpass, Fstop, dBpass, dBstop)
Bw = Fs/2;
[n, Wc] = buttord(Fpass/Bw,Fstop/Bw,dBpass,dBstop);

if Fpass<Fstop
    highlow = 'low';
else
    highlow = 'high';
end

fprintf(1, 'Applying a %s-pass filter, order = %d.\n', highlow, n);


[zf, p, k] = butter(n, Wc, highlow);
[sos, gg] = zp2sos(zf, p, k);
Hd = dfilt.df2tsos(sos, gg);
[B, A] = sos2tf(Hd.sosMatrix, Hd.Scalevalues);
y = filtfilt(B, A, z')';

end