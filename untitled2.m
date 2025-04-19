%-------- 基本パラメータの設定 ------------
function [param] = define_parameters()
    param.l = 2.000e-3;      % 探針長さ (m)
    param.lc = 0.15e-3;      % 探針エッチング部分の長さ (m)
    param.Eq = 8.00e10;      % クオーツのヤング率 (Pa)
    param.L = 2.353e-3;      % QTF prongの長さ (m)
    param.Le = 1.16e-3;      % QTFの電極の長さ (m)
    param.h = 2.23e-4;       % QTFの幅 (m)
    param.b_h = 1.11e-4;     % QTFの厚さ (m)
    param.lb = param.l - param.h - param.lc;
    param.rho_q = 2.650e3;   % クオーツの密度 (kg/m^3)
    param.Aq = param.h * param.b_h;  % クオーツ断面積 (m^2)
    param.Iq = param.b_h * param.h^3 / 12;  % クオーツ断面二次モーメント (m^4)
    param.Ew = 3.45e11;      % タングステンのヤング率 (Pa)
    param.d = 1.00e-4;       % 探針の直径 (m)
    param.rho_w = 1.93e4;    % タングステンの密度 (kg/m^3)
    param.Aw = pi * param.d^2 / 4;  % タングステン断面積 (m^2)
    param.Ma = param.Aw * param.h * param.rho_w;  % 探針接着部の質量
    param.Mb = param.Aw * param.lb * param.rho_w; % 探針変形部の質量
    param.Mc = 1/3 * param.rho_w * param.Aw * param.lc; % 探針のエッチング部分の質量
    param.Mtip = param.Ma + param.Mb + param.Mc;   % 探針質量 (kg)
    param.Iw = pi * param.d^4 / 64;  % タングステン断面二次モーメント (m^4)
    param.Ia = param.Ma * (param.d^2 / 16 + param.h^2 / 12);  % 探針接着部の慣性モーメント
    param.Ic = 3/80 * param.Mc * (param.d^2 + param.lc^2);  % 探針先端部の慣性モーメント
end

%-------- 共振周波数計算 ------------
function val = compute_detM(fre, param)
    omega = 2 * pi * fre;
    a = ((omega^2 * param.rho_q * param.Aq) / (param.Eq * param.Iq))^(1/4);
    b = ((omega^2 * param.rho_w * param.Aw) / (param.Ew * param.Iw))^(1/4);
    
    % 角度計算
    sina = sin(a * param.L); cosa = cos(a * param.L); sinha = sinh(a * param.L); cosha = cosh(a * param.L);
    sinb = sin(b * param.lb); cosb = cos(b * param.lb); sinhb = sinh(b * param.lb); coshb = cosh(b * param.lb);

    % M行列の成分
    M11 = -a * (param.h / 2) * (cosa - cosha);
    M12 = -a * (param.h / 2) * (-sina - sinha);
    M21 = -a * (cosa - cosha);
    M22 = -a * (-sina - sinha);
    M31 = param.Eq * param.Iq * a^3 * (-cosa - cosha) + param.Mtip * omega^2 * (sina - sinha);
    M32 = param.Eq * param.Iq * a^3 * (sina - sinha) + param.Mtip * omega^2 * (cosa - cosha);
    % 他の成分を省略 (同様の処理)

    % 行列の構築
    M = [M11 M12 0; M21 M22 0; M31 M32 0];  % 簡略化
    val = double(det(M));
end

%-------- root 探索本体 ------------
function [roots_found] = find_roots(param)
    f_min = 100;           % 最小周波数を100Hzに設定
    f_max = 1e6;           % 最大周波数を1MHzに設定
    N_points = 100000;
    f_list = linspace(f_min, f_max, N_points); % 新しい周波数範囲
    tol = 1e-5;            % 許容誤差を少し緩和
    roots_found = [];
    prev_val = compute_detM(f_list(1), param);
    prev_sign = sign(prev_val);

    for i = 2:N_points
        f_curr = f_list(i);
        D_curr = compute_detM(f_curr, param);
        curr_sign = sign(D_curr);

        if prev_sign ~= 0 && curr_sign ~= 0 && curr_sign ~= prev_sign
            f_left = f_list(i-1);
            f_right = f_list(i);
            while (f_right - f_left) > tol
                f_mid = (f_left + f_right) / 2;
                D_mid = compute_detM(f_mid, param);
                if sign(D_mid) == sign(compute_detM(f_left, param))
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
end

% 主な実行部分の更新
param = define_parameters(); % パラメータを定義
roots_found = find_roots(param); % 根を探索

% 根をチェックして不足している場合に警告を表示
if length(roots_found) < 3
    warning('期待した3つの根が見つかりませんでした。見つかった根の数: %d', length(roots_found));
end

% 根をソート
roots_found = sort(roots_found);

% 根が見つかった場合に表示
if length(roots_found) >= 1
    f_fre_1 = roots_found(1);
    fprintf('f_fre_1: %.8f Hz\n', f_fre_1);
end
if length(roots_found) >= 2
    f_fre_2 = roots_found(2);
    fprintf('f_fre_2: %.8f Hz\n', f_fre_2);
end
if length(roots_found) >= 3
    f_fre_3 = roots_found(3);
    fprintf('f_fre_3: %.8f Hz\n', f_fre_3);
end
