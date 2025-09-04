clc
% 这是输入磁场角α，输出磁尖偏转角θ的静态转向模型。

% 1.输入参数
B0=(1:1:70)'/1000;                                %外加磁场在真空中测得的感应强度,T为单位
% B0=70/1000;                                         %mT
H=B0*1000/0.4/pi;                                   %外加磁场的真实磁场强度，kA/m为单位.
alpha=30;                                           %磁场倾斜角度

% 2.导管参数，计算Q
D1=4*10^(-3);                                       %外管外径
D2=3.8*10^(-3);                                     %外管内径
D3=1.5*10^(-3);                                     %内管外径
D4=1*10^(-3);                                       %内管内径
L=50*10^(-3);                                       %导管长度
A=pi/4*(D2^2-D3^2);                                 %空腔底面积
I1=pi/64*(D1^4-D2^4);                               %外管截面惯性矩Iz=∫x^2*dA.
I2=pi/64*(D3^4-D4^4);                               %内管截面惯性矩Iz=∫x^2*dA.
G=0.046*10^6;                                       %基质剪切模量
niu=0.5;                                            %基质泊松比
V=A*L;                                              %导管的容积
Q=3*pi*(1+niu)*G*(D1^4-D2^4+D3^4-D4^4)/(32*L);      %自定义弹性系数
% Q=1;                                              %自定义弹性系数

% 3.磁粉参数，计算P
r=125*10^(-6);                                      %铁粉半径,最终表达式与颗粒半径无关，因此这个r多大不重要
D=2*r;                                              %两颗粒初始距离
fai=0.523;                                            %铁粉的体积分数
kai0=7.1;                                           %铁粉初始磁化率
miu0=4*pi*10^(-7);                                  %真空磁导率
Ms=2.156/miu0;                                      %铁粉的饱和磁化强度
M=kai0*Ms.*H*1000./(kai0.*H*1000+Ms);                 %Frohlich–Kennelly law
m=M*4/3*pi*r.^3;                                    %单个颗粒的磁偶极矩
miu=0.05;                                           %此处的μ，指颗粒间的摩擦系数
Np=fai*V./(4/3*pi*r.^3);                            %导管内颗粒数目
NL=L./(2*r);                                        %导管内颗粒层数
Nc=Np./NL;                                          %导管内单链根数
P=miu0*fai*pi*(D2^2-D3^2)/96*(kai0*Ms.*B0./(miu0*Ms+kai0.*B0)).^2;              %自定义磁性系数
% P=100;                                              %自定义磁性系数
Q_P_ratio=Q./P;                                      %QP比


% 4.GPT5数值求解θ，输出结果；

% % ------- direct solve: fsolve (angles in degrees) -------
% alph = deg2rad(alpha); % radians
% delta1 = (Q * sin(alph)) / max(2 * P + Q * cos(alph), eps); % initial offset
% theta0 = alph - delta1; % initial guess (rad)
% theta0 = min(max(theta0, 1e-9), alph - 1e-9); % clamp to (0, alph)
% 
% F = @(th) P * sin(2 * (alph - th)) - Q * sin(th); % residual (rad)
% 
% opts = optimoptions('fsolve', ...
% 'Display','off', ...
% 'FunctionTolerance',1e-12, ...
% 'OptimalityTolerance',1e-12, ...
% 'StepTolerance',1e-12, ...
% 'MaxIterations',400);
% 
% [th, ~, exitflag] = fsolve(F, theta0, opts);
% 
% if ~(exitflag > 0 && th > 0 && th < alph)
% error('Failed to solve for theta in (0, alpha).');
% end
% 
% theta = rad2deg(th); % result in degrees
% fprintf('theta = %.10f degrees\n', theta);
% % -------------------------------------------------------

%向量版
% ====== Solve P(i)sin(2(alpha(i)-theta(i))) = Q(i)*sin(theta(i)) for vectors ======
% alpha, P, Q: scalar or vectors; angles in degrees; theta returned in degrees (n×1)

% 1) Normalize shapes to column
makecol = @(x) reshape(x, [], 1); % forces column
alpha = makecol(alpha);
P_in = P; Q_in = Q; % keep originals for debug
if ~isscalar(P), P = makecol(P); end
if ~isscalar(Q), Q = makecol(Q); end

% 2) Determine target length n from alpha (or from the longest non-scalar)
lens = [numel(alpha), numel(P), numel(Q)];
n = max(lens);

% 3) Broadcast scalars; validate vector lengths if non-scalar
if isscalar(alpha) && n>1, alpha = repmat(alpha, n, 1); end
if isscalar(P) && n>1, P = repmat(P, n, 1); end
if isscalar(Q) && n>1, Q = repmat(Q, n, 1); end

% Recompute lens after possible broadcasting
lens = [numel(alpha), numel(P), numel(Q)];

% 4) Validate compatibility: any non-scalar must have length n
if any(lens ~= n)
fprintf('Size(alpha)=%s, Size(P)=%s, Size(Q)=%s\n', mat2str(size(alpha)), mat2str(size(P_in)), mat2str(size(Q_in)));
error('alpha, P, Q must be scalars or vectors of the same length (after row->column normalization).');
end

% 5) Now solve per element (use fsolve with fallback to fzero)
alph = deg2rad(alpha); % radians
theta = nan(n,1);

den = 2 .* P + Q .* cos(alph);
% protect tiny denominators elementwise
mask = abs(den) < eps;
den(mask) = den(mask) + eps .* sign(den(mask) + (den(mask)==0));

delta1 = (Q .* sin(alph)) ./ den;
theta0 = alph - delta1;
theta0 = min(max(theta0, 1e-9), alph - 1e-9);

opts = optimoptions('fsolve', ...
'Display','off', ...
'FunctionTolerance',1e-12, ...
'OptimalityTolerance',1e-12, ...
'StepTolerance',1e-12, ...
'MaxIterations',400);

for i = 1:n
a_i = alph(i); P_i = P(i); Q_i = Q(i);
F_i = @(th) P_i * sin(2 * (a_i - th)) - Q_i * sin(th);

[th_i, ~, exitflag] = fsolve(F_i, theta0(i), opts);

if ~(exitflag > 0 && th_i > 0 && th_i < a_i)
    % bracket + fzero fallback
    Nscan = 200;
    lo = a_i * 1e-9; hi = a_i * (1 - 1e-9);
    tt = linspace(lo, hi, Nscan);
    Ft = arrayfun(F_i, tt);
    k = find(Ft(1:end-1) .* Ft(2:end) <= 0, 1, 'first');
    if ~isempty(k)
        th_i = fzero(F_i, [tt(k), tt(k+1)]);
    else
        th_i = fzero(F_i, theta0(i));
    end
end

if ~(th_i > 0 && th_i < a_i)
    error('Failed to solve at index %d (alpha=%.6g deg, P=%.6g, Q=%.6g).', i, rad2deg(a_i), P_i, Q_i);
end

theta(i) = rad2deg(th_i);
end

% ================== 新增：在 B0 或 alpha 方向的对比绘图（容忍常数向量，兼容 MATLAB 语法） ==================
% 依赖：alpha, B0, P, Q, theta 已存在；theta 为数值解（度）
toCol = @(x) reshape(x, [], 1);
isConstVec = @(x) numel(x)>1 && all(x(:)==x(1)); % 所有元素完全相等则视作常量

% 一个安全的长度对齐函数（替代三目运算）
function y = toLen_safe(x, n)
if isscalar(x)
y = repmat(x, n, 1);
else
y = reshape(x, [], 1);
if numel(y) ~= n
% 如果长度不匹配且 x 是常数向量，允许压缩为标量再扩展
if numel(y) > 1 && all(y == y(1))
y = repmat(y(1), n, 1);
else
error('toLen_safe: length mismatch. Expected %d, got %d.', n, numel(y));
end
end
end
end

% 规范形状
a = toCol(alpha);
b = toCol(B0);
th_num = toCol(theta);

% 将"常数向量"降级为标量视角
alpha_is_const = isConstVec(a);
B0_is_const = isConstVec(b);

if alpha_is_const, a = a(1); end
if B0_is_const, b = b(1); end

isA_vec = numel(a) > 1; % 真正的扫描维度
isB_vec = numel(b) > 1;

if ~(xor(isA_vec, isB_vec))
if ~isA_vec && ~isB_vec
fprintf('alpha 与 B0 均为常量，跳过曲线绘制。theta = %.6g deg\n', th_num);
else
warning('alpha 与 B0 不能同时为非常量向量，请确保只扫描其中一个。');
end
return
end

if isA_vec
x = toCol(a); xlab = '\alpha (deg)';
n = numel(x);
% 广播 P、Q 与 theta
Pv = toLen_safe(P, n);
Qv = toLen_safe(Q, n);
thv = toLen_safe(th_num, n);
% 一阶近似
alph_rad = deg2rad(x);
den = 2 .* Pv + Qv .* cos(alph_rad);
m = abs(den) < eps; den(m) = den(m) + eps .* sign(den(m) + (den(m)==0));
th1 = rad2deg(alph_rad - (Qv .* sin(alph_rad)) ./ den);
else
x = toCol(b); xlab = 'B_0';
n = numel(x);
Pv = toLen_safe(P, n);
Qv = toLen_safe(Q, n);
thv = toLen_safe(th_num, n);
% alpha 作为常量或与 B0 对齐
if isscalar(alpha) || isConstVec(alpha)
alph_rad = deg2rad(alpha(1)) * ones(n,1);
else
alph_rad = deg2rad(toLen_safe(alpha, n));
end
den = 2 .* Pv + Qv .* cos(alph_rad);
m = abs(den) < eps; den(m) = den(m) + eps .* sign(den(m) + (den(m)==0));
th1 = rad2deg(alph_rad - (Qv .* sin(alph_rad)) ./ den);
end

% 排序与绘图
[xs, idx] = sort(x);
th_num_s = thv(idx);
th1_s = th1(idx);

abs_err = abs(th_num_s - th1_s);
rel_err = abs_err ./ max(abs(th_num_s), 1e-12);

figure('Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(xs, th_num_s, 'b-', 'LineWidth', 1.8); hold on;
plot(xs, th1_s, 'r--', 'LineWidth', 1.8);
grid on; xlabel(xlab); ylabel('\theta (deg)');
title('\theta 数值解 vs 一阶近似'); legend('数值解','一阶近似','Location','best');

nexttile;
yyaxis left; plot(xs, abs_err, 'k-', 'LineWidth', 1.6); ylabel('绝对误差 | \Delta\theta | (deg)');
yyaxis right; plot(xs, rel_err*100, 'm--', 'LineWidth', 1.6); ylabel('相对误差 (%)');
grid on; xlabel(xlab); title('误差对比');

fprintf('误差统计：max |Δθ| = %.4g deg, median |Δθ| = %.4g deg, max 相对误差 = %.3g %%\n', ...
max(abs_err), median(abs_err), max(rel_err)*100);
% ======================================================================