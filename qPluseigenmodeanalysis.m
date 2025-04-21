% --- 基本パラメータの設定 ---
% 探針長さを固定 (3 mm)
startl_point_end_l = 3.0e-3:0.1e-3:3.5e-3;%探針の長さ(m) 初期値：間隔：終了値
param.lc = 0.15e-3;     % 探針エッチング部分の長さ (m)
param.Eq = 8.00e10;     % クオーツのヤング率 (Pa)
param.L = 2.353e-3; %QTFprongの長さ (m)
param.Le = 1.16e-3; %QTFの電極の長さ (m)
param.h = 2.23e-4; %QTFの幅 (m)
param.b_h = 1.11e-4; %QTFの厚さ (m)
param.rho_q = 2.650e3;  % クオーツの密度 (kg/m^3)
param.Ew = 3.45e11;    % タングステンのヤング率 (Pa)
param.d = 1.00e-4;    % 探針の直径 (m)  
param.rho_w = 1.93e4; % タングステンの密度 (kg/m^3)    
% 結果保存用
results = [];

for l_i = startl_point_end_l
    param.l = l_i;
    param.Aq = param.h * param.b_h;  % クオーツ断面積 (m^2)
    param.lb = param.l - param.h - param.lc ;
    param.Iq = param.b_h * param.h^3 / 12;  % クオーツ断面二次モーメント (m^4) 
    param.Aw = pi * param.d^2 / 4; % タングステン断面積 (m^2)
    param.Ma = param.Aw * param.h * param.rho_w; %探針接着部の質量
    param.Mb = param.Aw * param.lb * param.rho_w; %探針変形部の質量
    param.Mc = 1/3 * param.rho_w * param.Aw * param.lc;     % 探針のエッチングされた部分の質量 (kg)
    param.Mtip = param.Ma + param.Mb + param.Mc;   % 探針質量 (kg) 
    param.Iw = pi * param.d^4 / 64; % タングステン断面二次モーメント (m^4)
    param.Ia = param.Ma * (param.d^2/16 + param.h^2 /12);     % 探針接着部の慣性モーメント
    param.Ic = 3/80 * param.Mc * (param.d^2 + param.lc^2);    % 探針先端部の慣性モーメント
    param.Itip = param.Ia + 1/3 * param.Mb * (param.lb^2 + 3/2 * param.lb * param.h + 3/4 * param.h^2) + ...
         param.Mc * (param.h/2 + param.lb + 1/4 * param.lc)^2 + param.Ic + (param.d/2)^2 * param.Mtip;
    % 探針全体の慣性モーメント(m^4, 近似(0 < x2 and x2 < lb)の慣性モーメント無視, 中心;
   %x1 = L, z1 = 0 (x2 = -h/2);
% ----------- root 探索本体 ------------
% 探索範囲
f_min = 1e3;
f_max = 1000e3;
N_points = 1000000;
f_list = linspace(f_min, f_max, N_points);

tol = 1e-8;
roots_found = [];

prev_val = compute_detM66(f_list(1),param);
prev_sign = sign(prev_val);

for i = 2:N_points
    f_curr = f_list(i);
    D_curr = compute_detM66(f_curr,param);
    curr_sign = sign(D_curr);

    if prev_sign ~= 0 && curr_sign ~= 0 && curr_sign ~= prev_sign
        f_left = f_list(i-1);
        f_right = f_list(i);

        while (f_right - f_left) > tol
            f_mid = (f_left + f_right) / 2;
            D_mid = compute_detM66(f_mid,param);

            if sign(D_mid) == sign(compute_detM66(f_left,param))
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
% f_fre_1 に対して X を求める
[M3x3_1, B1] = compute_M3x3_and_B(f_fre_1,param);
X1 = M3x3_1 \ B1;
disp('X1 ='); disp(X1);

% f_fre_2 に対して
[M3x3_2, B2] = compute_M3x3_and_B(f_fre_2,param);
X2 = M3x3_2 \ B2;
disp('X2 ='); disp(X2);

% f_fre_3 に対して
[M3x3_3, B3] = compute_M3x3_and_B(f_fre_3,param);
X3 = M3x3_3 \ B3;
disp('X3 ='); disp(X3);

% f_fre_1 に対してk_eff_v,k_eff_l,Angleを求める
[k1_eff_v, k1_eff_l,angle1] = compute_k_theta(X1(1),X1(2),X1(3),X1(4),X1(5),f_fre_1,param);

% 結果の表示
disp('eigenmode1')
fprintf('k1_eff_v = %.3f\n', k1_eff_v);
fprintf('k1_eff_v = %.3f\n', k1_eff_v);
fprintf('Angle1 = %.3f\n', angle1);

% f_fre_2 に対してk_eff_v,k_eff_l,Angleを求める
[k2_eff_v, k2_eff_l,angle2] = compute_k_theta(X2(1),X2(2),X2(3),X2(4),X2(5),f_fre_2,param);

% 結果の表示
disp('eigenmode2')
fprintf('k2_eff_v = %.3f\n', k2_eff_v);
fprintf('k2_eff_v = %.3f\n', k2_eff_v);
fprintf('Angle2 = %.3f\n', angle2);

% f_fre_3 に対してk_eff_v,k_eff_l,Angleを求める
[k3_eff_v, k3_eff_l,angle3] = compute_k_theta(X3(1),X3(2),X3(3),X3(4),X3(5),f_fre_3,param);

% 結果の表示
disp('eigenmode3')
fprintf('k3_eff_v = %.3f\n', k2_eff_v);
fprintf('k3_eff_v = %.3f\n', k2_eff_v);
fprintf('Angle3 = %.3f\n', angle3);

% 結果を保存
results = [results; double([param.l*1e3, f_fre_1, k1_eff_v, k1_eff_l,angle1, X1(1),X1(2), ...
    X1(3),X1(4),X1(5),f_fre_2, k2_eff_v, k2_eff_l,angle2,X2(1),X2(2), X2(3),X2(4),X2(5)...
    f_fre_3, k3_eff_v, k3_eff_l,angle3,X3(1),X3(2), X3(3),X3(4),X3(5)])];  % 単位:mm, Hz

end

% テーブル表示
T = array2table(results, 'VariableNames', {'l_mm', 'f1_Hz', 'k1_eff_v', 'k1_eff_l','angle1','X1(1)','X1(2)', ...
    'X1(3)','X1(4)','X1(5)','f2_Hz', 'k2_eff_v', 'k2_eff_l','angle2','X2(1)','X2(2)', ...
    'X2(3)','X2(4)','X2(5)','f3_Hz','k3_eff_v', 'k3_eff_l','angle3','X3(1)','X3(2)', ...
    'X3(3)','X3(4)','X3(5)'});
disp(T)

%-------- compute_detMの定義 ------------
function val_M66 = compute_detM66(freq_M,param)
omega_M = 2 * pi * freq_M;
a_M = ((omega_M^2 * param.rho_q * param.Aq) / (param.Eq * param.Iq))^(1/4);
b_M = ((omega_M^2 * param.rho_w * param.Aw) / (param.Ew * param.Iw))^(1/4);
sina_M = sin(a_M*param.L);
cosa_M = cos(a_M*param.L);
sinha_M = sinh(a_M*param.L);
cosha_M = cosh(a_M*param.L);
sinb_M = sin(b_M*param.lb);
cosb_M = cos(b_M*param.lb);
sinhb_M = sinh(b_M*param.lb);
coshb_M = cosh(b_M*param.lb);

% 各成分を数値として計算 (付録A.2に基づく)
M11_M = -a_M * (param.h/2)*(cosa_M - cosha_M);
M12_M = -a_M * (param.h/2)*(-sina_M - sinha_M);
M13_M = 0; M14_M = 1; M15_M = 0; M16_M = 1;

M21_M = -a_M * (cosa_M - cosha_M);
M22_M = -a_M * (-sina_M - sinha_M);
M23_M = b_M; M24_M = 0; M25_M = b_M; M26_M = 0;

M31_M = param.Eq * param.Iq * a_M^3 * (-cosa_M - cosha_M) + ...
           param.Mtip * omega_M^2 * (sina_M - sinha_M) + ...
           param.Mtip * omega_M^2 * (param.d/2)*a_M*(cosa_M-cosha_M);
M32_M = param.Eq * param.Iq * a_M^3 * (sina_M - sinha_M) + param.Mtip * omega_M^2 * (cosa_M - cosha_M) ...
    + param.Mtip * omega_M^2 * (param.d/2)*a_M*(-sina_M-sinha_M);
M33_M = 0; M34_M = 0; M35_M = 0; M36_M = 0;

M41_M = -param.Eq * param.Iq * a_M^2 * (-sina_M - sinha_M) + param.d/2 * param.Mtip * omega_M^2 * (sina_M - sinha_M)...
          +((param.d/2)^2 * param.Mtip + param.Ia) * omega_M^2 * a_M * (cosa_M - cosha_M);
M42_M = -param.Eq * param.Iq * a_M^2 * (-cosa_M-cosha_M)+param.d/2 * param.Mtip * omega_M^2 * (cosa_M - cosha_M)...
           +((param.d/2)^2 * param.Mtip + param.Ia) * omega_M^2 * a_M * (-sina_M - sinha_M);
M43_M = (param.Ic + 0.25 * param.lc * param.Mc * (param.lb + param.h / 2 + 0.25 * param.lc)) * omega_M^2 * b_M * cosb_M + ...
          (param.lb + param.h / 2 + 0.25 * param.lc)* param.Mc * omega_M^2 * sinb_M +...
          param.rho_w * param.Aw * omega_M^2 * (-1 / b_M * (param.lb + param.h / 2) * cosb_M+ ...
          1/b_M^2*sin(b_M*param.lb)+param.h/(2*b_M));
M44_M = -(param.Ic + 0.25 * param.lc * param.Mc * (param.lb + param.h / 2 + 0.25 * param.lc)) * omega_M^2 * b_M * sinb_M...
        + (param.lb + param.h / 2 + 0.25 * param.lc)* param.Mc * omega_M^2 * cosb_M +...
          param.rho_w * param.Aw * omega_M^2 * (1 / b_M * (param.lb + param.h / 2) * sinb_M + 1 / b_M^2 * cosb_M - 1 / b_M^2 );
M45_M = (param.Ic + 0.25 * param.lc * param.Mc * (param.lb + param.h / 2 + 0.25 * param.lc)) * omega_M^2 * b_M * coshb_M + ...
          (param.lb + param.h / 2 + 0.25 * param.lc)* param.Mc * omega_M^2 * sinhb_M +...
          param.rho_w * param.Aw * omega_M^2 * (1 / b_M * (param.lb + param.h / 2) * coshb_M - 1 / b_M^2 * sinhb_M - param.h / (2 * b_M) );
M46_M = (param.Ic + 0.25 * param.lc * param.Mc * (param.lb + param.h / 2 + 0.25 * param.lc)) * omega_M^2 * b_M * sinhb_M + ...
          (param.lb + param.h / 2 + 0.25 * param.lc)* param.Mc * omega_M^2 * coshb_M +...
          param.rho_w * param.Aw * omega_M^2 * (1 / b_M * (param.lb + param.h / 2) * sinhb_M - 1 / b_M^2 * coshb_M + 1 /  b_M^2 );
M51_M = 0;M52_M = 0;
M53_M = -param.Ew * param.Iw * b_M^3 * cosb_M+ 0.25 * param.lc * param.Mc * omega_M^2 * b_M * cosb_M + param.Mc * omega_M^2 * sinb_M; 
M54_M = param.Ew * param.Iw * b_M^3 * sinb_M - 0.25 * param.lc * param.Mc * omega_M^2 * b_M * sinb_M + param.Mc * omega_M^2 * cosb_M; 
M55_M = (param.Ew * param.Iw * b_M^3 + 0.25 * param.lc * param.Mc * omega_M^2 * b_M) * coshb_M + param.Mc * omega_M^2 * sinhb_M; 
M56_M = (param.Ew * param.Iw * b_M^3 + 0.25 * param.lc * param.Mc * omega_M^2 * b_M) * sinhb_M + param.Mc * omega_M^2 * coshb_M;

M61_M = 0;M62_M = 0;
M63_M = (param.Ew * param.Iw * b_M^2 + 0.25 * param.lc * param.Mc * omega_M^2) * sinb_M + b_M * ((0.25 * param.lc)^2 * param.Mc + param.Ic) * omega_M^2 * cosb_M; 
M64_M = (param.Ew * param.Iw * b_M^2 + 0.25 * param.lc * param.Mc * omega_M^2) * cosb_M - b_M * ((0.25 * param.lc)^2 * param.Mc + param.Ic) * omega_M^2 * sinb_M; 
M65_M = (-param.Ew * param.Iw * b_M^2 + 0.25 * param.lc * param.Mc * omega_M^2) * sinhb_M + b_M * ((0.25 * param.lc)^2 * param.Mc + param.Ic) * omega_M^2 * coshb_M; 
M66_M = (-param.Ew * param.Iw * b_M^2 + 0.25 * param.lc * param.Mc * omega_M^2) * coshb_M + b_M * ((0.25 * param.lc)^2 * param.Mc + param.Ic) * omega_M^2 * sinhb_M;

    % 行列 M の構築
M66 = [M11_M M12_M M13_M M14_M M15_M M16_M;
     M21_M M22_M M23_M M24_M M25_M M26_M;
     M31_M M32_M M33_M M34_M M35_M M36_M;
     M41_M M42_M M43_M M44_M M45_M M46_M;
     M51_M M52_M M53_M M54_M M55_M M56_M;
     M61_M M62_M M63_M M64_M M65_M M66_M];
   
    val_M66 = double(det(M66));
end

function [M3x3, B] = compute_M3x3_and_B(freq_M,param)
    omega3 = 2 * pi * freq_M;

    a3= ((omega3^2 * param.rho_q * param.Aq) / (param.Eq * param.Iq))^(1/4);
    b3= ((omega3^2 * param.rho_w * param.Aw) / (param.Ew * param.Iw))^(1/4);

    % 基本三角関数
    sina3 = sin(a3*param.L);
    cosa3 = cos(a3*param.L); 
    sinha3 = sinh(a3*param.L); 
    cosha3 = cosh(a3*param.L);
    sinb3 = sin(b3*param.lb); 
    cosb3 = cos(b3*param.lb); 
    sinhb3 = sinh(b3*param.lb); 
    coshb3 = cosh(b3*param.lb);

    % 使いたい成分のみを再定義
    M11_3 = -a3 * (param.h/2)*(cosa3-cosha3);
    M21_3 = -a3 * (cosa3 - cosha3);
    M31_3 = param.Eq * param.Iq * a3^3 * (-cosa3 - cosha3) + ...
           param.Mtip * omega3^2 * (sina3 - sinha3) + ...
           param.Mtip * omega3^2 * (param.d / 2) * a3 * (cosa3 - cosha3);

    M12_3 = -a3 * (param.h/2)*(-sina3 - sinha3);
    M22_3 = -a3 * (-sina3 - sinha3);
    M32_3 = param.Eq * param.Iq * a3^3 * (sina3 - sinha3) + param.Mtip * omega3^2 * (cosa3 - cosha3)...
        + param.Mtip * omega3^2 * (param.d /2) * a3 * (-sina3 - sinha3);

    M53_3 = -param.Ew * param.Iw * b3^3 * cosb3 + 0.25 * param.lc * param.Mc * omega3^2 * b3 * cosb3 + param.Mc * omega3^2 * sinb3; 
    M54_3 = param.Ew * param.Iw * b3^3 * sinb3 - 0.25 * param.lc * param.Mc * omega3^2 * b3 * sinb3 + param.Mc * omega3^2 * cosb3; 
    M55_3 = (param.Ew * param.Iw * b3^3 + 0.25 * param.lc * param.Mc * omega3^2 * b3) * coshb3 + param.Mc * omega3^2 * sinhb3; 
    M56_3 = (param.Ew * param.Iw * b3^3 + 0.25 * param.lc * param.Mc * omega3^2 * b3) * sinhb3 + param.Mc * omega3^2 * coshb3;

    M63_3 = (param.Ew * param.Iw * b3^2 + 0.25 * param.lc * param.Mc * omega3^2) * sinb3 + b3 * ((0.25 * param.lc)^2 * param.Mc + param.Ic) * omega3^2 * cosb3; 
    M64_3 = (param.Ew * param.Iw * b3^2 + 0.25 * param.lc * param.Mc * omega3^2) * cosb3 - b3 * ((0.25 * param.lc)^2 * param.Mc + param.Ic) * omega3^2 * sinb3; 
    M65_3 = (-param.Ew * param.Iw * b3^2 + 0.25 * param.lc * param.Mc * omega3^2) * sinhb3 + b3 * ((0.25 * param.lc)^2 * param.Mc + param.Ic) * omega3^2 * coshb3; 
    M66_3 = (-param.Ew * param.Iw * b3^2 + 0.25 * param.lc * param.Mc * omega3^2) * coshb3 + b3 * ((0.25 * param.lc)^2 * param.Mc + param.Ic) * omega3^2 * sinhb3;

    % Bベクトルと M3x3行列の構築（5x5）
    B = [-M11_3; -M21_3; -M31_3; 0; 0];
    M3x3 = [M12_3, 0, 1, 0, 1;
            M22_3, b3, 0, b3, 0;
            M32_3, 0, 0, 0, 0;
            0, M53_3, M54_3, M55_3, M56_3;
            0, M63_3, M64_3, M65_3, M66_3];
end
%共振1
function [k_eff_v, k_eff_l,angle_abs] = compute_k_theta(X1,X2,X3,X4,X5,f,param)
omegak = 2 * pi * f;
ak = ((omegak^2 * param.rho_q * param.Aq) / (param.Eq * param.Iq))^(1/4);
bk = ((omegak^2 * param.rho_w * param.Aw) / (param.Ew * param.Iw))^(1/4);
C1 = 1;
C2 = X1;
D1 = X2;
D2 = X3;
D3 = X4;
D4 = X5;
%syms x1 z1;
wL = C1 * (sin(ak*param.L))-sinh(ak*param.L) + C2 * (cos(ak*param.L)-cosh(ak*param.L));
w = @(x1)(C1*(sin(ak*x1)-sinh(ak*x1)))+ (C2*(cos(ak*x1)-cosh(ak*x1)))/wL;
u = @(z1)(D1*sin(bk*z1)+D2*cos(bk*z1)+D3*sinh(bk*z1)+D4*cosh(bk*z1))/wL;
%プロット
%figure;
%fplot(w, [param.L]);
%xlabel('x / mm');
%ylabel('w / mm');
%grid on;
%figure;
%fplot(u, [0 param.l]);
%xlabel('x / mm');
%ylabel('w / mm');grid on;

A1_val = 1 + (1/2)*ak*(C1*(cos(ak*param.L) - cosh(ak*param.L)) + C2*(-sin(ak*param.L) - sinh(ak*param.L)))*param.d/wL;
dw = @(Le) ak*(C1*(cos(ak*param.Le) - cosh(ak*param.Le)) + C2*(-sin(ak*param.Le) - sinh(ak*param.Le)))/wL;
% A2 の計算
A2 = @(lb) (D1*sin(bk*lb)+D2*cos(bk*lb)+D3*sinh(bk*lb)+D4*cosh(bk*lb))/wL + ...
    (bk*(D1*cos(bk*lb) - D2*sin(bk*lb) + D3*cosh(bk*lb) + D4*sinh(bk*lb)))/wL*param.lc;
A2_val = A2(param.lb);
% dw1 と dw2 の計算
dw1 = dw(param.Le) / A1_val;
dw2 = dw(param.Le) / A2_val;

% d2w と d2u の計算
d2w = @(x) ak^2*(C1*(-sin(ak*x) - sinh(ak*x)) + C2*(-cos(ak*x) - cosh(ak*x))) / wL;
d2u = @(z) bk^2*(-D1*sin(bk*z) - D2*cos(bk*z) + D3*sinh(bk*z) + D4*cosh(bk*z)) / wL;

% V の計算
V = 1/2*param.Eq*param.Iq*integral(@(x) d2w(x).^2, 0, param.L) + 1/2*param.Ew*param.Iw*integral(@(z) d2u(z).^2, 0, param.lb);

% k1 と k2 の計算
k_eff_v = 2*V / A1_val^2;
k_eff_l = 2*V / A2_val^2;

% angle の計算
angle = 180 * atan(A1_val / A2_val) / pi;
angle_abs = abs(angle);
end