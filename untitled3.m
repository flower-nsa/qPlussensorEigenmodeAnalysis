% -------- compute_detMの定義 ------------
function val = compute_detM(fre)
    % --- 基本パラメータの設定 ---
    l = 2.000e-3;  % 探針長さ (m)
    lc = 0.15e-3;  % 探針エッチング部分の長さ (m)
    Eq = 8.00e10;  % クオーツのヤング率 (Pa)
    L = 2.353e-3;  % QTFprongの長さ (m)
    h = 2.23e-4;   % QTFの幅 (m)
    b_h = 1.11e-4; % QTFの厚さ (m)
    lb = l - h - lc;
    rho_q = 2.650e3;  % クオーツの密度 (kg/m^3)
    Aq = h * b_h;     % クオーツ断面積 (m^2)
    Iq = b_h * h^3 / 12;  % クオーツ断面二次モーメント (m^4)
    Ew = 3.45e11;    % タングステンのヤング率 (Pa)
    d = 1.00e-4;     % 探針の直径 (m)
    rho_w = 1.93e4;  % タングステンの密度 (kg/m^3)
    Aw = pi * d^2 / 4;  % タングステン断面積 (m^2)
    Ma = Aw * h * rho_w;   % 探針接着部の質量
    Mb = Aw * lb * rho_w;  % 探針変形部の質量
    Mc = 1 / 3 * rho_w * Aw * lc;  % 探針エッチング部分の質量
    Mtip = Ma + Mb + Mc;  % 探針質量
    Iw = pi * d^4 / 64;  % タングステン断面二次モーメント (m^4)
    Ia = Ma * (d^2 / 16 + h^2 / 12);  % 探針接着部の慣性モーメント
    Ic = 3 / 80 * Mc * (d^2 + lc^2);  % 探針先端部の慣性モーメント
    omega = 2 * pi * fre;  % 角周波数

    a = ((omega^2 * rho_q * Aq) / (Eq * Iq))^(1 / 4);
    b = ((omega^2 * rho_w * Aw) / (Ew * Iw))^(1 / 4);

    % 基本三角関数
    sina = sin(a * L); cosa = cos(a * L); sinha = sinh(a * L); cosha = cosh(a * L);
    sinb = sin(b * lb); cosb = cos(b * lb); sinhb = sinh(b * lb); coshb = cosh(b * lb);

    % 各成分の計算 (行列 M の構築)
    M11 = -a * (h / 2) * (cosa - cosha);
    M12 = -a * (h / 2) * (-sina - sinha);
    M21 = -a * (cosa - cosha);
    M22 = -a * (-sina - sinha);
    M31 = Eq * Iq * a^3 * (-cosa - cosha) + ...
           Mtip * omega^2 * (sina - sinha) + ...
           Mtip * omega^2 * (d / 2) * a * (cosa - cosha);
    M32 = Eq * Iq * a^3 * (sina - sinha) + Mtip * omega^2 * (cosa - cosha) + Mtip * omega^2 * (d / 2) * a * (-sina - sinha);

    % 他の成分も同様に計算
    % ...（省略）...

    % 行列 M の構築
    M = [M11 M12 M13 M14 M15 M16;
     M21 M22 M23 M24 M25 M26;
     M31 M32 M33 M34 M35 M36;
     M41 M42 M43 M44 M45 M46;
     M51 M52 M53 M54 M55 M56;
     M61 M62 M63 M64 M65 M66];

    % 行列式を計算
    val = det(M);  % detは元々double精度なのでここではキャストしない
end

% ----------- root 探索本体 ------------
f_min = 1e3;  % 最小周波数
f_max = 300e3;  % 最大周波数
N_points = 1000000;  % サンプル数
f_list = linspace(f_min, f_max, N_points);  % 周波数リスト

tol = 1e-8;  % 許容誤差
roots_found = [];  % 見つかった根を格納
prev_val = compute_detM(f_list(1));
prev_sign = sign(prev_val);

for i = 2:N_points
    f_curr = f_list(i);
    D_curr = compute_detM(f_curr);
    curr_sign = sign(D_curr);

    if prev_sign ~= 0 && curr_sign ~= 0 && curr_sign ~= prev_sign
        % 二分法で根を探索
        f_left = f_list(i - 1);
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

        % 中央値が根
        root = (f_left + f_right) / 2;
        roots_found(end + 1) = root;
        fprintf('%.8f Hz の根を発見しました。\n', root);

        if length(roots_found) == 3
            break;  % 3つの根が見つかったら終了
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

% 各種パラメータを使用して計算
% 以下、X1, X2, X3を求める
[M3x3_1, B1] = compute_M3x3_and_B(f_fre_1);
X1 = M3x3_1 \ B1;
disp('X1 ='); disp(X1);

[M3x3_2, B2] = compute_M3x3_and_B(f_fre_2);
X2 = M3x3_2 \ B2;
disp('X2 ='); disp(X2);

[M3x3_3, B3] = compute_M3x3_and_B(f_fre_3);
X3 = M3x3_3 \ B3;
disp('X3 ='); disp(X3);
