clc;
clear all;
close all;

% Susigeneruojam dviejų sluoksnių neuronų kolekcijas su rand reikšmėmis
wx_1 = [randn(1),randn(1),randn(1),randn(1),randn(1),randn(1)];
wx_2 = [randn(1),randn(1),randn(1),randn(1),randn(1),randn(1)];

%pirmo sluoksnio rezultatai
bx_1 = [randn(1),randn(1),randn(1),randn(1),randn(1),randn(1)];

%antro sluoksnio rezultatas
b1_2 = randn(1);

for k = 0:1:10000
    for x = 0.1:1/50:1
        % Sugeneruojame d d
        d = (1 + 0.6*sin(2*pi*x/0.7)) + 0.3*sin(2*pi*x)/2;
    
        % Skaiciuojam pirmo sluoksnio isejimus
        vx_1 = [];
        for i = 1:1:length(wx_1)
            vx_1(i)=x*wx_1(i)+bx_1(i);
        end        
        
        % Pritaikom aktyviaja funkcija
        yx_1 = [];
        for i = 1:1:length(wx_1)
            yx_1(i)=1/(1+exp(-vx_1(i)));
        end
    
        % Skaiciuojam antro sluoksnio isejimus
        v = 0;
        for i = 1:1:length(wx_2)
            v = v + yx_1(i)*wx_2(i);
        end

        y=v;
        
        % Skaiciuojam klaida, delta_out = e
        e=d-y;
    
        % Atnaujinam koeficientus: w = w + n * phi * y, n bet koks
        n=0.15;
        for i = 1:1:length(wx_2)
            wx_2(i) = wx_2(i)+n*e*yx_1(i);
        end

        b1_2 = b1_2 + n * e;
        
        % Atnaujinam pirmo sluoksnio koeficientus        
        deltax_1 = [];
        for i = 1:1:length(wx_2)
            deltax_1(i)=(yx_1(i)*(1-yx_1(i)))*e*wx_2(i);
        end
        
        for i = 1:1:length(wx_1)           
            wx_1(i)=wx_1(i)+n*deltax_1(i)*x;
        end
    
        for i = 1:1:length(deltax_1)           
        bx_1(i)=bx_1(i)+n*deltax_1(i);
        end
    end
end

yf = [];
i = 0;

for x = 0.1:1/50:1
    % Apsirašome ir apskaičiuojame d
    d = (1 + 0.6 * sin(2*pi*x/0.7)) + 0.3 * sin(2*pi*x) / 2;

    % Skaiciuojam pirmo sluoksnio isejimus
    vx_1 = [];
    for i = 1:1:length(wx_1)
        vx_1(i)=x*wx_1(i)+bx_1(i);
    end

    % Pritaikom aktyviaja funkcija
    yx_1 = [];
    for i = 1:1:length(vx_1)
        yx_1(i)=1/(1+exp(-vx_1(i)));
    end
    
    % Skaiciuojam antro sluoksnio isejimus
    y = 0;
    for i = 1:1:length(wx_2)
        y = y + yx_1(i) * wx_2(i);
    end
    yf = [yf, y];
end

x = 0.1:1/50:1;
d = (1 + 0.6 * sin(2*pi*x/0.7)) + 0.3 * sin(2*pi*x) / 2;
plot ( x, d, 'r*', x, yf, 'bx')