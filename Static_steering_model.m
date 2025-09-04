clc
% ��������ų��Ǧ�������ż�ƫת�Ǧȵľ�̬ת��ģ�͡�

% 1.�������
B0=(1:1:70)'/1000;                                %��Ӵų�������в�õĸ�Ӧǿ��,TΪ��λ
% B0=70/1000;                                         %mT
H=B0*1000/0.4/pi;                                   %��Ӵų�����ʵ�ų�ǿ�ȣ�kA/mΪ��λ.
alpha=30;                                           %�ų���б�Ƕ�

% 2.���ܲ���������Q
D1=4*10^(-3);                                       %����⾶
D2=3.8*10^(-3);                                     %����ھ�
D3=1.5*10^(-3);                                     %�ڹ��⾶
D4=1*10^(-3);                                       %�ڹ��ھ�
L=50*10^(-3);                                       %���ܳ���
A=pi/4*(D2^2-D3^2);                                 %��ǻ�����
I1=pi/64*(D1^4-D2^4);                               %��ܽ�����Ծ�Iz=��x^2*dA.
I2=pi/64*(D3^4-D4^4);                               %�ڹܽ�����Ծ�Iz=��x^2*dA.
G=0.046*10^6;                                       %���ʼ���ģ��
niu=0.5;                                            %���ʲ��ɱ�
V=A*L;                                              %���ܵ��ݻ�
Q=3*pi*(1+niu)*G*(D1^4-D2^4+D3^4-D4^4)/(32*L);      %�Զ��嵯��ϵ��
% Q=1;                                              %�Զ��嵯��ϵ��

% 3.�ŷ۲���������P
r=125*10^(-6);                                      %���۰뾶,���ձ��ʽ������뾶�޹أ�������r�����Ҫ
D=2*r;                                              %��������ʼ����
fai=0.523;                                            %���۵��������
kai0=7.1;                                           %���۳�ʼ�Ż���
miu0=4*pi*10^(-7);                                  %��մŵ���
Ms=2.156/miu0;                                      %���۵ı��ʹŻ�ǿ��
M=kai0*Ms.*H*1000./(kai0.*H*1000+Ms);                 %Frohlich�CKennelly law
m=M*4/3*pi*r.^3;                                    %���������Ĵ�ż����
miu=0.05;                                           %�˴��Ħ̣�ָ�������Ħ��ϵ��
Np=fai*V./(4/3*pi*r.^3);                            %�����ڿ�����Ŀ
NL=L./(2*r);                                        %�����ڿ�������
Nc=Np./NL;                                          %�����ڵ�������
P=miu0*fai*pi*(D2^2-D3^2)/96*(kai0*Ms.*B0./(miu0*Ms+kai0.*B0)).^2;              %�Զ������ϵ��
% P=100;                                              %�Զ������ϵ��
Q_P_ratio=Q./P;                                      %QP��


% 4.GPT5��ֵ���ȣ���������

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

%������
% ====== Solve P(i)sin(2(alpha(i)-theta(i))) = Q(i)*sin(theta(i)) for vectors ======
% alpha, P, Q: scalar or vectors; angles in degrees; theta returned in degrees (n��1)

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

% ================== �������� B0 �� alpha ����ĶԱȻ�ͼ�����̳������������� MATLAB �﷨�� ==================
% ������alpha, B0, P, Q, theta �Ѵ��ڣ�theta Ϊ��ֵ�⣨�ȣ�
toCol = @(x) reshape(x, [], 1);
isConstVec = @(x) numel(x)>1 && all(x(:)==x(1)); % ����Ԫ����ȫ�������������

% һ����ȫ�ĳ��ȶ��뺯���������Ŀ���㣩
function y = toLen_safe(x, n)
if isscalar(x)
y = repmat(x, n, 1);
else
y = reshape(x, [], 1);
if numel(y) ~= n
% ������Ȳ�ƥ���� x �ǳ�������������ѹ��Ϊ��������չ
if numel(y) > 1 && all(y == y(1))
y = repmat(y(1), n, 1);
else
error('toLen_safe: length mismatch. Expected %d, got %d.', n, numel(y));
end
end
end
end

% �淶��״
a = toCol(alpha);
b = toCol(B0);
th_num = toCol(theta);

% ��"��������"����Ϊ�����ӽ�
alpha_is_const = isConstVec(a);
B0_is_const = isConstVec(b);

if alpha_is_const, a = a(1); end
if B0_is_const, b = b(1); end

isA_vec = numel(a) > 1; % ������ɨ��ά��
isB_vec = numel(b) > 1;

if ~(xor(isA_vec, isB_vec))
if ~isA_vec && ~isB_vec
fprintf('alpha �� B0 ��Ϊ�������������߻��ơ�theta = %.6g deg\n', th_num);
else
warning('alpha �� B0 ����ͬʱΪ�ǳ�����������ȷ��ֻɨ������һ����');
end
return
end

if isA_vec
x = toCol(a); xlab = '\alpha (deg)';
n = numel(x);
% �㲥 P��Q �� theta
Pv = toLen_safe(P, n);
Qv = toLen_safe(Q, n);
thv = toLen_safe(th_num, n);
% һ�׽���
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
% alpha ��Ϊ�������� B0 ����
if isscalar(alpha) || isConstVec(alpha)
alph_rad = deg2rad(alpha(1)) * ones(n,1);
else
alph_rad = deg2rad(toLen_safe(alpha, n));
end
den = 2 .* Pv + Qv .* cos(alph_rad);
m = abs(den) < eps; den(m) = den(m) + eps .* sign(den(m) + (den(m)==0));
th1 = rad2deg(alph_rad - (Qv .* sin(alph_rad)) ./ den);
end

% �������ͼ
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
title('\theta ��ֵ�� vs һ�׽���'); legend('��ֵ��','һ�׽���','Location','best');

nexttile;
yyaxis left; plot(xs, abs_err, 'k-', 'LineWidth', 1.6); ylabel('������� | \Delta\theta | (deg)');
yyaxis right; plot(xs, rel_err*100, 'm--', 'LineWidth', 1.6); ylabel('������ (%)');
grid on; xlabel(xlab); title('���Ա�');

fprintf('���ͳ�ƣ�max |����| = %.4g deg, median |����| = %.4g deg, max ������ = %.3g %%\n', ...
max(abs_err), median(abs_err), max(rel_err)*100);
% ======================================================================