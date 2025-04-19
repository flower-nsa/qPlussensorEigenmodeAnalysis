%-------- compute_detMの定義 ------------
function val = compute_detM(fre)
% --- 基本パラメータの設定 ---
% 探針長さを固定 (3 mm)
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
    Mc * (h/2 + lb + 1/4 * lc)^2 + Ic + (d/2)^2 * Mtip;    
% 探針全体の慣性モーメント(m^4, 近似(0 < x2 and x2 < lb)の慣性モーメント無視, 中心;
%x1 = L, z1 = 0 (x2 = -h/2);
omega = 2 * pi * fre;
a = ((omega^2 * rho_q * Aq) / (Eq * Iq))^(1/4);
b = ((omega^2 * rho_w * Aw) / (Ew * Iw))^(1/4);
sina = sin(a*L);
cosa = cos(a*L);
sinha = sinh(a*L);
cosha = cosh(a*L);
sinb = sin(b*lb);
cosb = cos(b*lb);
sinhb = sinh(b*lb);
coshb = cosh(b*lb);

% 各成分を数値として計算 (付録A.2に基づく)
M11 = -a * (h/2)*(cosa - cosha);
M12 = -a * (h/2)*(-sina - sinha);
M13 = 0; M14 = 1; M15 = 0; M16 = 1;

M21 = -a * (cosa - cosha);
M22 = -a * (-sina - sinha);
M23 = b; M24 = 0; M25 = b; M26 = 0;

M31 = Eq * Iq * a^3 * (-cosa - cosha) + ...
           Mtip * omega^2 * (sina - sinha) + ...
           Mtip * omega^2 * (d / 2) * a * (cosa - cosha);
M32 = Eq * Iq * a^3 * (sina - sinha) + Mtip * omega^2 * (cosa - cosha) + Mtip * omega^2 * (d /2) * a * (-sina - sinha);
M33 = 0; M34 = 0; M35 = 0; M36 = 0;

M41 = -Eq * Iq * a^2 * (-sina - sinha) + d/2 * Mtip * omega^2 * (sina - sinha)...
          +((d/2)^2 * Mtip + Ia) * omega^2 * a * (cosa - cosha);
M42 = -Eq * Iq * a^2 * (-cosa - cosha) + d/2 * Mtip * omega^2 * (cosa - cosha)...
          +((d/2)^2 * Mtip + Ia) * omega^2 * a * (-sina - sinha);
M43 = (Ic + 0.25 * lc * Mc * (lb + h / 2 + 0.25 * lc)) * omega^2 * b * cosb + ...
          (lb + h / 2 + 0.25 * lc)* Mc * omega^2 * sinb +...
          rho_w * Aw * omega^2 * (-1 / b * (lb + h / 2) * cosb+ ...
          1 / b^2 * sin(b * lb) + h / (2 * b));
M44 = -(Ic + 0.25 * lc * Mc * (lb + h / 2 + 0.25 * lc)) * omega^2 * b * sinb + (lb + h / 2 + 0.25 * lc)* Mc * omega^2 * cosb +...
          rho_w * Aw * omega^2 * (1 / b * (lb + h / 2) * sinb + 1 / b^2 * cos(b * lb) - 1 / b^2 );
M45 = (Ic + 0.25 * lc * Mc * (lb + h / 2 + 0.25 * lc)) * omega^2 * b * coshb + ...
          (lb + h / 2 + 0.25 * lc)* Mc * omega^2 * sinhb +...
          rho_w * Aw * omega^2 * (1 / b * (lb + h / 2) * coshb - 1 / b^2 * sinhb - h / (2 * b) );
M46 = (Ic + 0.25 * lc * Mc * (lb + h / 2 + 0.25 * lc)) * omega^2 * b * sinhb + ...
          (lb + h / 2 + 0.25 * lc)* Mc * omega^2 * coshb +...
          rho_w * Aw * omega^2 * (1 / b * (lb + h / 2) * sinhb - 1 / b^2 * coshb + 1 /  b^2 );
M51 = 0;M52 = 0;
M53 = -Ew * Iw * b^3 * cosb+ 0.25 * lc * Mc * omega ^2 * b * cosb + Mc * omega^2 * sinb; 
M54 = Ew * Iw * b^3 * sinb - 0.25 * lc * Mc * omega ^2 * b * sinb + Mc * omega^2 * cosb; 
M55 = (Ew * Iw * b^3 + 0.25 * lc * Mc * omega ^2 * b) * coshb + Mc * omega^2 * sinhb; 
M56 = (Ew * Iw * b^3 + 0.25 * lc * Mc * omega ^2 * b) * sinhb + Mc * omega^2 * coshb;

M61 = 0;
M62 = 0;
M63 = (Ew * Iw * b^2 + 0.25 * lc * Mc * omega^2) * sinb + b * ((0.25 * lc)^2 * Mc + Ic) * omega^2 * cosb; 
M64 = (Ew * Iw * b^2 + 0.25 * lc * Mc * omega^2) * cosb - b * ((0.25 * lc)^2 * Mc + Ic) * omega^2 * sinb; 
M65 = (-Ew * Iw * b^2 + 0.25 * lc * Mc * omega^2) * sinhb + b * ((0.25 * lc)^2 * Mc + Ic) * omega^2 * coshb; 
M66 = (-Ew * Iw * b^2 + 0.25 * lc * Mc * omega^2) * coshb + b * ((0.25 * lc)^2 * Mc + Ic) * omega^2 * sinhb;

    % 行列 M の構築
M = [M11 M12 M13 M14 M15 M16;
     M21 M22 M23 M24 M25 M26;
     M31 M32 M33 M34 M35 M36;
     M41 M42 M43 M44 M45 M46;
     M51 M52 M53 M54 M55 M56;
     M61 M62 M63 M64 M65 M66];
   
    val = double(det(M));
end

% ----------- root 探索本体 ------------
% 探索範囲
f_min = 1e3;
f_max = 300e3;
N_points = 100000000;
f_list = linspace(f_min, f_max, N_points);

tol = 1e-9;
roots_found = [];
prev_val = compute_detM(f_list(1));
prev_sign = sign(prev_val);

for i = 2:N_points
    f_curr = f_list(i);
    D_curr = compute_detM(f_curr);
    curr_sign = sign(D_curr);

    if prev_sign ~= 0 && curr_sign ~= 0 && curr_sign ~= prev_sign
        f_left = f_list(i-1);
        f_right = f_list(i);

        while (f_right - f_left) > tol
            f_mid = (f_left + f_right) / 2;
            D_mid = compute_detM(f_mid);

            if sign(D_mid) == sign(compute_detM(f_left))
                f_left = f_mid;
            else
                f_right = f_mid;
            end
        end

        root = (f_left + f_right) / 2;
        roots_found(end+1) = root;
        fprintf('%.8f Hz の根を発見しました。\n', root);

        if length(roots_found) == 3
            break;
        end
    end
    prev_sign = curr_sign;
end

% 見つかった根を小さい順にソート
roots_found = sort(roots_found);

% f_fre_1, f_fre_2, f_fre_3 に代入
f_fre_1 = roots_found(1);
f_fre_2 = roots_found(2);
f_fre_3 = roots_found(3);

fprintf('\n--- 見つかった根（周波数） ---\n');
fprintf('f_fre_1: %.8f Hz\n', f_fre_1);
fprintf('f_fre_2: %.8f Hz\n', f_fre_2);
fprintf('f_fre_3: %.8f Hz\n', f_fre_3);
%ここまで共振周波数の計算





function [M3x3, B] = compute_M3x3_and_B(fre)
    % 上記と同じパラメータをここでも定義
    l = 2.000e-3;
    lc = 0.15e-3;
    Eq = 8.00e10;
    L = 2.353e-3;
    h = 2.23e-4;
    b_h = 1.11e-4;
    lb = l - h - lc ;
    rho_q = 2.650e3;
    Aq = h * b_h;
    Iq = b_h * h^3 / 12;
    Ew = 3.45e11;
    d = 1.00e-4;
    rho_w = 1.93e4;
    Aw = pi * d^2 / 4;
    Ma = Aw * h * rho_w;
    Mb = Aw * lb * rho_w;
    Mc = 1/3 * rho_w * Aw * lc;
    Mtip = Ma + Mb + Mc;
    Iw = pi * d^4 / 64;
    Ia = Ma * (d^2/16 + h^2 /12);
    Ic = 3/80 * Mc * (d^2 + lc^2);
    omega = 2 * pi * fre;

    a = ((omega^2 * rho_q * Aq) / (Eq * Iq))^(1/4);
    b_val = ((omega^2 * rho_w * Aw) / (Ew * Iw))^(1/4);

    % 基本三角関数
    sina = sin(a*L); cosa = cos(a*L); sinha = sinh(a*L); cosha = cosh(a*L);
    sinb = sin(b_val*lb); cosb = cos(b_val*lb); sinhb = sinh(b_val*lb); coshb = cosh(b_val*lb);

    % 使いたい成分のみを再定義
    M11 = -a * (h/2)*(cosa - cosha);
    M21 = -a * (cosa - cosha);
    M31 = Eq * Iq * a^3 * (-cosa - cosha) + ...
           Mtip * omega^2 * (sina - sinha) + ...
           Mtip * omega^2 * (d / 2) * a * (cosa - cosha);

    M12 = -a * (h/2)*(-sina - sinha);
    M22 = -a * (-sina - sinha);
    M32 = Eq * Iq * a^3 * (sina - sinha) + Mtip * omega^2 * (cosa - cosha) + Mtip * omega^2 * (d /2) * a * (-sina - sinha);

    M53 = -Ew * Iw * b_val^3 * cosb + 0.25 * lc * Mc * omega^2 * b_val * cosb + Mc * omega^2 * sinb; 
    M54 = Ew * Iw * b_val^3 * sinb - 0.25 * lc * Mc * omega^2 * b_val * sinb + Mc * omega^2 * cosb; 
    M55 = (Ew * Iw * b_val^3 + 0.25 * lc * Mc * omega^2 * b_val) * coshb + Mc * omega^2 * sinhb; 
    M56 = (Ew * Iw * b_val^3 + 0.25 * lc * Mc * omega^2 * b_val) * sinhb + Mc * omega^2 * coshb;

    M63 = (Ew * Iw * b_val^2 + 0.25 * lc * Mc * omega^2) * sinb + b_val * ((0.25 * lc)^2 * Mc + Ic) * omega^2 * cosb; 
    M64 = (Ew * Iw * b_val^2 + 0.25 * lc * Mc * omega^2) * cosb - b_val * ((0.25 * lc)^2 * Mc + Ic) * omega^2 * sinb; 
    M65 = (-Ew * Iw * b_val^2 + 0.25 * lc * Mc * omega^2) * sinhb + b_val * ((0.25 * lc)^2 * Mc + Ic) * omega^2 * coshb; 
    M66 = (-Ew * Iw * b_val^2 + 0.25 * lc * Mc * omega^2) * coshb + b_val * ((0.25 * lc)^2 * Mc + Ic) * omega^2 * sinhb;

    % Bベクトルと M3x3行列の構築（5x5）
    B = [-M11; -M21; -M31; 0; 0];
    M3x3 = [M12, 0, 1, 0, 1;
            M22, b_val, 0, b_val, 0;
            M32, 0, 0, 0, 0;
            0, M53, M54, M55, M56;
            0, M63, M64, M65, M66];
end

% f_fre_1 に対して X を求める
[M3x3_1, B1] = compute_M3x3_and_B(f_fre_1);
X1 = M3x3_1 \ B1;
disp('X1 ='); disp(X1);

% f_fre_2 に対して
[M3x3_2, B2] = compute_M3x3_and_B(f_fre_2);
X2 = M3x3_2 \ B2;
disp('X2 ='); disp(X2);

% f_fre_3 に対して
[M3x3_3, B3] = compute_M3x3_and_B(f_fre_3);
X3 = M3x3_3 \ B3;
disp('X3 ='); disp(X3);

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

%共振1
omega1 = 2 * pi * f_fre_1;
a = ((omega1^2 * rho_q * Aq) / (Eq * Iq))^(1/4);
b = ((omega1^2 * rho_w * Aw) / (Ew * Iw))^(1/4);
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
disp('eigenmode1')
fprintf('k1 = %.3f\n', k1);
fprintf('k2 = %.3f\n', k2);
fprintf('Angle = %.3f\n', angle);

%共振2
omega2 = 2 * pi * f_fre_2;
a = ((omega2^2 * rho_q * Aq) / (Eq * Iq))^(1/4);
b = ((omega2^2 * rho_w * Aw) / (Ew * Iw))^(1/4);
C1 = 1;
C2 = X2(1);
D1 = X2(2);
D2 = X2(3);
D3 = X2(4);
D4 = X2(5);
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
disp('eigenmode2')
fprintf('k1 = %.3f\n', k1);
fprintf('k2 = %.3f\n', k2);
fprintf('Angle = %.3f\n', angle);

%共振3
omega3 = 2 * pi * f_fre_3;
a = ((omega3^2 * rho_q * Aq) / (Eq * Iq))^(1/4);
b = ((omega3^2 * rho_w * Aw) / (Ew * Iw))^(1/4);
C1 = 1;
C2 = X3(1);
D1 = X3(2);
D2 = X3(3);
D3 = X3(4);
D4 = X3(5);
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
disp('eigenmode3')
fprintf('k1 = %.3f\n', k1);
fprintf('k2 = %.3f\n', k2);
fprintf('Angle = %.3f\n', angle);
