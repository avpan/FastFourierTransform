clear all;clc;clf;

l = [8;16;32];
N1 = [32;64;128];
m = 1;

for i = 1:3
    for j = 1:3
        L = l(i);
        N = N1(j);
        x = linspace(0,L,N+1);
        x = x(1:N);
        n = 1-N/2:N/2;
        kn = 2*pi/L*n;

        x0 = L/2;
        fn = exp(-(x-x0).^2);

        %fast fourier transform
        ffn = fft(fn)/N;
        ffn_hold = ffn(N/2+1:N);
        ffn(N/2+1:N) = ffn(1:N/2);
        ffn(1:N/2) = ffn_hold;

        power = ffn.*ffn;
        
        hold on;subplot(3,3,m);semilogy(kn,power,'-o');  
        xlabel('wave number (k_n)'); ylabel('|f_{n}|^{2}');
        m = m+1;
    end
end
