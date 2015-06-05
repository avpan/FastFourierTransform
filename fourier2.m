clear all;clc;%clf;

L = 16;
N1 = [32;64;128;256];
m = 1;

for j = 4:4
        N = N1(j);
        x = linspace(0,L,N+1);
        x = x(1:N);
        n = 1-N/2:N/2;
        nhold = n(N/2+1:N);
        n(N/2+1:N) = n(1:N/2);
        n(1:N/2) = nhold;
        kn = 2*pi/L*n;

        x0 = L/2;
        fn = exp(-(x-x0).^2);
        fnd = 2*(x0-x).*fn;
        
        %fast fourier transform
        ffn = fft(fn)./N;
        orig = ffn;
        ffn_hold = ffn(N/2+1:N);
        ffn(N/2+1:N) = ffn(1:N/2);
        ffn(1:N/2) = ffn_hold;

        power = ffn.*ffn;
        
        %deriffiate fft
        ffnd = i*kn.*orig;
        df = real(ifft(ffnd).*N);
        
        Error = abs(df-fnd)./fnd.*100;
        hold on; 
        subplot(2,1,m); plot(x,df,'x',x,fnd,'o')
        legend('FFT','Actual')
        xlabel('x'); ylabel('df(x)');
        subplot(2,1,m+1); plot(x,Error,'x');
         xlabel('x'); ylabel('P.Error (%)');
        m = m+2;
end

%plot(x,df,'*',x,fnd,'o')