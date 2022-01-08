% Script for OpAmp

function [] = script_OpAmp()
    j=33;

    % Values related to j
    Cl = (2 + 0.01*j) * 10^(-12);
    SR = (18 + 0.01*j) * 10^6;
    Vdd = 1.8 + 0.003*j;
    Vss = -(1.8 + 0.003*j);
    GB = (7 + 0.01*j) * 10^6;
    A = (20 + 0.01*j);
    P = (50 + 0.01*j) * 10^(-3);

    % Other values that we need
    uop = 180.2 * 10^(-4);
    kp = 2.9352 * 10^(-5);
    Cox = kp / uop;
    Ip = 0.05;
    Vtop = -0.9056;
    uon = 591.7 * 10^(-4);
    kn = 9.6379 * 10^(-5);
    In = 0.04;
    Vton = 0.786;
    Vin_max = 0.1;
    Vin_min = -0.1;
    L = 10^(-6);

    % Calculations
    Cc = 0.22 * Cl;
    I5 = Cc * SR;
    S3 = I5 / (kp * (Vdd - Vin_max - abs(Vtop) -0.15 + Vton -0.15)^2);
    S4 = S3;

    p3 = sqrt(2 * kp * S3 * I5/2)/(2 * 0.667 * S3 * (L^2) * Cox);
    gm1 = 2 * pi * GB * Cc;
    gm2 = gm1;
    S1 = gm1^2 / (kn * I5);
    S2 = S1;

    Vds5 = Vin_min - Vss - sqrt(I5 / (S1*kn)) - Vton -0.15;
    S5 = 2 * I5 / (kn * Vds5^2);
    S8 = S5;

    gm6 = 2.2 * gm2 * Cl/Cc;
    gm4 = sqrt(2 * kp * S4 * I5/2);
    S6 = S4 * gm6 / gm4;
    I6 = gm6^2 / (2 * kp * S6);
    S7 = I6 * S5 / I5;

    Av = 2 * gm2 * gm6 / (I5 * (In + Ip) * I6 * (In + Ip));
    Pdiss = (I5 + I6) * (Vdd + abs(Vss));

    W1 = S1 * L;
    W2 = S2 * L;
    W3 = S3 * L;
    W4 = S4 * L;
    W5 = S5 * L;
    W6 = S6 * L;
    W7 = S7 * L;
    W8 = S8 * L;

    % Show the results we need 
    
    Cc
    I5
    S3
    p3
    gm1
    S1
    Vds5
    S5
    gm6
    gm4
    S6
    I6
    S7
    Av
    Pdiss
    W1
    W2
    W3
    W4
    W5
    W6
    W7
    W8

end
