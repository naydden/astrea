function rho = rho_formules(a)
RE=6.378e6; %Earth Radius [m]
    H = (a - RE)*1e-3;
    % Compute exospheric temperature [K]
    T = 900 + 2.5*(100-70);
    % Compute effective atmospheric molecular mass [km/K], valid 180<H<500
    M = 27 - 0.012 * (H - 200);
    % Compute atmospheric scale height [km]
    SH = T / M;
    % Compute atmospheric density [kg/m3]
    rho = 6E-10 * exp(-(H - 175) / SH);
end