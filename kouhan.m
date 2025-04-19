%各種パラメータ
l = 2.000e-3;  % 探針長さ (m)
lc = 0.15e-3;     % 探針エッチング部分の長さ (m)
Eq = 8.00e10;     % クオーツのヤング率 (Pa)
L = 2.353e-3; %QTFprongの長さ (m)
Le = 1.16e-3; %QTFの電極の長さ (m)
h = 2.23e-4; %QTFの幅 (m)
b_h = 1.11e-4; %QTFの厚さ (m)
lb = l - h - lc ;
rho_q = 2.650e3;  % クオーツの密度 (kg/m^3)
Aq = h * b_h;  % クオーツ断面積 (m^2)
Iq = b_h * h^3 / 12;  % クオーツ断面二次モーメント (m^4) 
Ew = 3.45e11;    % タングステンのヤング率 (Pa)
d = 1.00e-4;    % 探針の直径 (m)  
rho_w = 1.93e4; % タングステンの密度 (kg/m^3)
Aw = pi * d^2 / 4; % タングステン断面積 (m^2)
Ma = Aw * h * rho_w; %探針接着部の質量
Mb = Aw * lb * rho_w; %探針変形部の質量
Mc = 1/3 * rho_w * Aw * lc;     % 探針のエッチングされた部分の質量 (kg)
Mtip = Ma + Mb + Mc;   % 探針質量 (kg)  
Iw = pi * d^4 / 64; % タングステン断面二次モーメント (m^4)
Ia = Ma * (d^2/16 + h^2 /12);     % 探針接着部の慣性モーメント
Ic = 3/80 * Mc * (d^2 + lc^2);    % 探針先端部の慣性モーメント
Itip = Ia + 1/3 * Mb * (lb^2 + 3/2 * lb * h + 3/4 * h^2) + ...
    Mc * (h/2 + lb + 1/4 * lc)^2 + Ic + (d/2)^2 * Mtip; X1 = [-1.2568, 0.8859, -0.3625, -0.0562, 0.4328];
f_fre_1 = 9712.98616869;
omega = 2 * pi * f_fre_1;
a = ((omega^2 * rho_q * Aq) / (Eq * Iq))^(1/4);
b = ((omega^2 * rho_w * Aw) / (Ew * Iw))^(1/4);


C1 = 1;
C2 = X1(1);
D1 = X1(2);
D2 = X1(3);
D3 = X1(4);
D4 = X1(5);
syms x1 z1;
wL = C1 * (sin(a*L))-sinh(a*L) + C2 * (cos(a*L)-cosh(a*L));
w = (C1*(sin(a*x1)-sinh(a*x1)))+ (C2*(cos(a*x1)-cosh(a*x1)))/wL;
u = (D1*sin(b*z1)+D2*cos(b*z1)+D3*sinh(b*z1)+D4*cosh(b*z1))/wL;
w_num = subs(w, {x1}, {x1});
u_num = subs(u, {z1}, {z1});
% プロット
%figure;
%fplot(matlabFunction(w_num), [0 L]);
%figure;
%fplot(matlabFunction(u_num), [0 l]);

A1 = 1 + (1/2)*a*(C1*(cos(a*L) - cosh(a*L)) + C2*(-sin(a*L) - sinh(a*L)))*d/wL;
dw = @(Le) a*(C1*(cos(a*Le) - cosh(a*Le)) + C2*(-sin(a*Le) - sinh(a*Le)))/wL;
% A2 の計算
A2 = @(lb) (D1*sin(b*lb)+D2*cos(b*lb)+D3*sinh(b*lb)+D4*cosh(b*lb))/wL + (b*(D1*cos(b*lb) - D2*sin(b*lb) + D3*cosh(b*lb) + D4*sinh(b*lb)))/wL*lc;
A2_val = A2(lb);
% dw1 と dw2 の計算
dw1 = dw(Le) / A1;
dw2 = dw(Le) / A2_val;

% d2w と d2u の計算
d2w = @(x) a^2*(C1*(-sin(a*x) - sinh(a*x)) + C2*(-cos(a*x) - cosh(a*x))) / wL;
d2u = @(z) b^2*(-D1*sin(b*z) - D2*cos(b*z) + D3*sinh(b*z) + D4*cosh(b*z)) / wL;

% V の計算
V = 1/2*Eq*Iq*integral(@(x) d2w(x).^2, 0, L) + 1/2*Ew*Iw*integral(@(z) d2u(z).^2, 0, lb);

% k1 と k2 の計算
k1 = 2*V / A1^2;
k2 = 2*V / A2_val^2;

% angle の計算
angle = 180 * atan(A1 / A2_val) / pi;

% 結果の表示
fprintf('k1 = %.3f\n', k1);
fprintf('k2 = %.3f\n', k2);
fprintf('Angle = %.3f\n', angle);



