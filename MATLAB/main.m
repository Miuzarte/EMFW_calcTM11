% 长宽 间隔
a = 20 * 1e-3;
b = 10 * 1e-3;
h = 1 * 1e-3;

% 误差阈值 最大迭代次数
tolerance = 1e-5;
maxIter = 1e5;

% 创建矩阵
nx = fix(a / h) + 1;
ny = fix(b / h) + 1;
Ez = zeros(nx, ny);

% 内部初始条件
Ez(2:end-1, 2:end-1) = 1;

% 更新Ez
function newEz = updateEz(Ez, h, kc)
    [rows, cols] = size(Ez);
    newEz = zeros(size(Ez));
    for i = 2:rows - 1
        for j = 2:cols - 1
            newEz(i, j) = (...
                Ez(i + 1, j) +...
                Ez(i - 1, j) +...
                Ez(i, j + 1) +...
                Ez(i, j - 1)...
                ) / (4 - kc * kc * h * h);
        end
    end
end

% 更新kc
function newKc = updateKc(Ez, h)
    kc_squared = zeros(size(Ez(2:end-1, 2:end-1)));
    for i = 2:size(Ez, 1) - 1
        for j = 2:size(Ez, 2) - 1
            kc_squared(i - 1, j - 1) =...
                -(Ez(i + 1, j) +...
                Ez(i - 1, j) +...
                Ez(i, j + 1) +...
                Ez(i, j - 1) -...
                4 * Ez(i, j)...
                ) / (Ez(i, j) * h * h);
        end
    end
    newKc = sqrt(mean(kc_squared(:)));
end

% 初始kc, 要大于0
kc = 0.1;
kc_old = kc + 2*tolerance;

n = 0;
while true
    n = n + 1;
    if n >= maxIter
        break;
    end

    Ez_old = Ez;
    Ez = updateEz(Ez, h, kc);

    if mod(n, 10) == 0
        kc_old = kc;
        kc = updateKc(Ez, h);
    end

    EzDiff = norm(Ez - Ez_old);
    KcDiff = norm(kc - kc_old);
    if EzDiff < tolerance || KcDiff < tolerance
        break;
    end
end

% 理论值
kc_theory = sqrt((pi / a) ^ 2 + (pi / b) ^ 2);
lambda_c_theory = 2 * pi / kc_theory;
f_c_theory = 3e8 / lambda_c_theory;

% 计算值
kc = updateKc(Ez, h);
lambda_c = 2 * pi / kc;
f_c = 3e8 / lambda_c;

fprintf('Ez:\n');
[r, c] = size(Ez);
for i = 1:r
    for j = 1:c
        fprintf('%.6f ', Ez(i, j));
    end
    fprintf('\n');
end
fprintf('kc: %.6f\n', kc);
fprintf('截止波长: %.6f mm\n', lambda_c);
fprintf('截止频率: %.6f Hz\n', f_c);
fprintf('理论kc: %.6f\n', kc_theory);
fprintf('理论截止波长: %.6f mm\n', lambda_c_theory);
fprintf('理论截止频率: %.6f Hz\n', f_c_theory);

surf(Ez');
xlim([1, nx]);
ylim([1, ny]);
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Ez');
title('TM_{11} E_z');
