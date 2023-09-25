function [gx,dgdx,dgdp] = g_shape_quick(x,P,u,in)
%G_SHAPE
% Launch with parameter (x,P,u,in) to use as a VBA observation function
% Launch without parameter to start a recovery analysis
% Launch with two parameter x,y to get a fit y(x), with output [P,y_hat,fit_measures]
% Launch with three parameter x,y,P_prev to get a fit y(x), with output
% [P,y_hat,fit_measures] starting from the priors P_prev
% (without gradient descent)

if nargin==0, launch_test(); return; end
if nargin==2, gx=estimate(x,P); return; end
if nargin==3, gx=estimate(x,P,u); return; end

%% transformation functions
beta = 50;
% posi_fun  = @(x)  log(1+exp(beta*x)) /beta + 1e-5;
    function out = posi_fun(x)
        out = log(1+exp(beta*x)) /beta + 1e-5;
        if isinf(out), out = x; end
    end
    function out = dposi_fun(x)
        out = exp(beta*x)./(exp(beta*x)+1);
        if isnan(out), out = 1; end
    end
shape_fun = @(x,c,p,s) (1-abs((x-c)/s).^p);
dshape_fun_dc = @(x,c,p,w) p*(x - c)  .* abs((x - c)/w).^(p - 2)/w^2;
dshape_fun_dw = @(x,c,p,w) p*(x - c).^2.* abs((x - c)/w).^(p - 2)/w^3;
dshape_fun_dp = @(x,c,p,w) -abs((x - c)/w).^p .* log(abs((x - c)/w));

%% unpacking parameters
c   = P(1);
A   = P(2);

w1  = posi_fun(P(3));
dw1 = dposi_fun(P(3));
w2  = posi_fun(P(6));
dw2 = dposi_fun(P(6));

p1  = exp(min(3,P(4)));
dp1 = exp(min(3,P(4)));
p2  = exp(min(3,P(7)));
dp2 = exp(min(3,P(7)));

b1  = P(5);
b2  = P(8);

gx = NaN(1,size(u,2));
dgdp = zeros(numel(P),size(u,2));

%% core function
idx_bf = find(u<c);
param = {u(idx_bf),c,p1,w1};
tmp_gx = shape_fun(param{:});
gx(idx_bf)  = (A-b1)*posi_fun(tmp_gx)+b1;
dgx = (A-b1)*dposi_fun(tmp_gx);

dgdp(1,idx_bf) = dshape_fun_dc(param{:}).*dgx;   % c
dgdp(2,idx_bf) = posi_fun(tmp_gx);              %A
dgdp(3,idx_bf) = dshape_fun_dw(param{:}).*dw1.*dgx;%w1
dgdp(4,idx_bf) = dshape_fun_dp(param{:}).*dp1.*dgx;%p1
dgdp(5,idx_bf) = -posi_fun(tmp_gx) + 1;% base level


idx_af = find(u>c);
param = {u(idx_af),c,p2,w2};
tmp_gx = shape_fun(param{:});
gx(idx_af)  = (A-b2)*posi_fun(tmp_gx)+b2;
dgx = (A-b2)*dposi_fun(tmp_gx);

dgdp(1,idx_af) = dshape_fun_dc(param{:}).*dgx; % c
dgdp(2,idx_af) = posi_fun(tmp_gx); % A
dgdp(6,idx_af) = dshape_fun_dw(param{:}).*dw2.*dgx; % w2
dgdp(7,idx_af) = dshape_fun_dp(param{:}).*dp2.*dgx; % p2
dgdp(8,idx_af) = -posi_fun(tmp_gx) + 1;% base level


idx_on = u==c;
param = {u(idx_on),c,p2,w2};
tmp_gx = shape_fun(param{:});
gx(idx_on)  = A*posi_fun(tmp_gx);
dgdx = [];
dgdp(2,idx_on) = 1;
end

function launch_test()
% recovery analysis
n_simu = 500;
n_trial = 256;
u_range = linspace(-0.10,3.10,n_trial);
SNR = 1/2;%1;%1/2;
rng(42)

%% simulate data
simu_param = NaN(8,n_simu);
randi_MinMax = @(m,M) randi(M-(m-1))+(m-1);
for ii=1:n_simu
    simu_param(2,ii) = randi_MinMax(-100,100)/100; %A
    simu_param(5,ii) = randi_MinMax(-100,100)/100; %b1
    simu_param(8,ii) = randi_MinMax(-100,100)/100; %b2
    simu_param(1,ii) = randi_MinMax(100,200)/100; %c
    simu_param(3,ii) = randi_MinMax(25,100)/100; %w1
    simu_param(6,ii) = randi_MinMax(25,100)/100; %w2
    simu_param(4,ii) = randi_MinMax(-100,100)/100; %p1
    simu_param(7,ii) = randi_MinMax(-100,100)/100; %p2
end

data = NaN(n_simu,n_trial);
dgdp = NaN(n_simu,n_trial,8);

%% fit
fitted_param = NaN(size(simu_param));
h=figure();
for ii_simu=1:n_simu
    % simulate data
    [data(ii_simu,:),~,tmp_dgdp] = g_shape_quick([],simu_param(:,ii_simu),u_range,[]);
    dgdp(ii_simu,:,:) = tmp_dgdp';
    % effective A
    [~,idx_c] = min(abs(u_range-simu_param(1,ii_simu)));
    simu_param(2,ii_simu) = data(ii_simu,idx_c);
    
    % fit data
    noise = 1/SNR*normrnd(0,std(data(ii_simu,:)),size(data(ii_simu,:)));
    y = data(ii_simu,:)+noise;
    
    P = estimate(u_range,y);
    fitted_param(:,ii_simu) = P;
    
    simu= NaN(size(y));
    for ii_t=1:numel(u_range)
        simu(ii_t) = g_shape_quick([],P',u_range(ii_t),[]);
    end
    R2 = 1-var(data(ii_simu,:)-simu)/var(data(ii_simu,:));
    
    if ii_simu>5 && mod(ii_simu,5)==0
        figure(h)
        set(h,'color','w')
        clf
        mat1 = simu_param(:,1:ii_simu)';
        mat2 = fitted_param(:,1:ii_simu)';
        
        subplot(2,1,1)
        imagesc( corr(mat1,mat2))
        title({'Correlations between generative and fitted parameters';...
            ['n_{simu}=',num2str(ii_simu), ...
            ', SNR=',num2str(SNR),...
            ', Sum(diag)=',num2str(sum(diag(corr(mat1,mat2))),3)]})
        ylabel('Generative')
        xlabel('Fitted')
        
        xticks(1:8)
        yticks(1:8)
        xticklabels({'c','A','w1','p_1','b_1','w2','p_2','b_2'})
        yticklabels({'c','A','w1','p_1','b_1','w2','p_2','b_2'})
        title(colorbar,'Correlation')
        colormap(cool)
        caxis([-1,1])
        axis tight
        axis square
        
        
        subplot(2,1,2)
        param_relevance = squeeze(mean(dgdp(1:ii_simu,:,:).^2,2));
        param_relevance = param_relevance./max(param_relevance,[],2);
        imagesc( corr(mat1.*param_relevance,mat2.*param_relevance))
        title({'Impact-weighted correlations';...
            ['n_{simu}=',num2str(ii_simu),...
            ', SNR=',num2str(SNR) ]})
        ylabel('Generative')
        xlabel('Fitted')
        
        xticks(1:8)
        yticks(1:8)
        xticklabels({'c','A','w1','p_1','b_1','w2','p_2','b_2'})
        yticklabels({'c','A','w1','p_1','b_1','w2','p_2','b_2'})
        title(colorbar,'Correlation')
        colormap(cool)
        caxis([-1,1])
        axis tight
        axis square
        
    end
end

end

function [P,y_hat] = estimate(x,y,P_prev)
n=numel(x);
x=reshape(x,1,n);
y=reshape(y,1,n);

fft_y = fft(y);
fft_y( abs(fft_y) < prctile(abs(fft_y),33)) = 0;
denoised_signal = ifft(fft_y);

if nargin==2
    [P,y_hat] = estimate_subfun(x,denoised_signal);
elseif nargin==3
    [P,y_hat] = estimate_subfun(x,denoised_signal,P_prev);
end
end

function [P,y_hat] = estimate_subfun(x,y,prev_P)
%% format
tmp = arrayfun(@(x){x},1:8);
[idx_c,idx_A,idx_w1,idx_p1,idx_b1,idx_w2,idx_p2,idx_b2] =deal(tmp{:});

n = numel(x);

%%  initialisation
current_param = NaN(8,1);

smooth_kernel = ones(round(numel(y)/10),1);
smooth_kernel = smooth_kernel / sum(smooth_kernel);
y_smoothed = conv(y,smooth_kernel,'valid');
y_smoothed = [NaN(1,floor((numel(y)-numel(y_smoothed))/2)),...
    y_smoothed,...
    NaN(1,ceil((numel(y)-numel(y_smoothed))/2),1)];

% b1 & b2
b1 = nanmean(y_smoothed(1:round(end/5)));
b2 = nanmean(y_smoothed(end-round(end/5)+1:end));
y_smoothed(isnan(y_smoothed(1:round(end/2)))) = b1;
y_smoothed(round(n/2)-1+find(isnan(y_smoothed(round(n/2):n)))) = b2;

% A & c
signal_above = y_smoothed - max(b1,b2);
signal_below = min(b1,b2) - y_smoothed;

above_strength = sum(signal_above(signal_above>0));
below_strength = sum(signal_below(signal_below>0));

[is_above,is_below,is_middle] = deal(false);
if above_strength>3*below_strength, is_above = true; end
if below_strength>3*above_strength, is_below = true; end

if ~is_above && ~is_below,          is_middle = true; end

if is_middle, ii_c = round(numel(y)/2); A = y_smoothed(ii_c); end

if is_above, [A,ii_c] = max(y_smoothed);  end
if is_below, [A,ii_c] = min(y_smoothed); end

c = x(ii_c);

% w1 & w2
mean_dev_b1 = std(y_smoothed(1:round(end/5)));
mean_dev_b2 = std(y_smoothed(end-round(end/5)+1:end));

dev_from_b1 = abs(y_smoothed-b1) > 3*mean_dev_b1;
dev_from_b2 = abs(y_smoothed-b2) > 3*mean_dev_b2;

n_sum = round(numel(y)/5);
neighboor_dev_from_b1 = conv(dev_from_b1,[zeros(n_sum,1);ones(n_sum,1)],'same');
[~,ii_w1] = min(fliplr(neighboor_dev_from_b1(1:ii_c))); % get the last minima before the peak (i.e. when it starts to deviate)
ii_w1 = ii_c - ii_w1;
w1 = c - x(ii_w1);

neighboor_dev_from_b2 = conv(dev_from_b2,[ones(n_sum,1);zeros(n_sum,1)],'same');
[~,ii_w2] = min(neighboor_dev_from_b2(ii_c:end)); % get the first minima after the peak (i.e. when it stops to deviate)
ii_w2 = ii_c + ii_w2;
w2 = x(ii_w2) - c;

% p1 & p2
p1 = 0;
p2 = 0;

% store param
current_param(idx_b1) = b1;
current_param(idx_b2) = b2;

current_param(idx_c) = c;
current_param(idx_A) = A;

current_param(idx_w1) = w1;
current_param(idx_w2) = w2;

current_param(idx_p1) = p1;
current_param(idx_p2) = p2;


%% linear approx
ii_count = 0;
ii_count_max = 10;
lambda = 1;
original_param = current_param;

% random sampling around the refe point and 1-step gradient descent
n_random_draw = 200;
random_var_P = [current_param,sign(current_param).*max(abs(current_param),0.5).*(1 + normrnd(0,1,size(current_param,1),n_random_draw))];

if exist('prev_P','var')
    random_var_P(:,2+ceil(n_random_draw/2):end)= sign(prev_P).*max(abs(prev_P),0.5).*(1 + normrnd(0,1,size(prev_P,1),ceil(n_random_draw/2)));
    random_var_P(:,end+1) = prev_P;
end

r2s_var = NaN(size(random_var_P,2),1);
for ii_var = 1:numel(r2s_var)
    [simu,~,dgdP] = g_shape_quick([],random_var_P(:,ii_var),x,[]);
    while ii_count < ii_count_max
        
        X_lin = [dgdP';lambda*eye(size(dgdP,1))];
        if ii_var>1+ceil( n_random_draw/2) || ~exist('prev_P','var')
            P_delta = pinv(X_lin'*X_lin)*X_lin'*([y-simu,lambda*(original_param-random_var_P(:,ii_var))'])';
        else
            P_delta = pinv(X_lin'*X_lin)*X_lin'*([y-simu,lambda*(prev_P-random_var_P(:,ii_var))'])';
        end
        
        [simu2,~,dgdP] = g_shape_quick([],random_var_P(:,ii_var) + P_delta,x,[]);
        if var(y-simu) > var(y-simu2)
            random_var_P(:,ii_var) = random_var_P(:,ii_var) + P_delta;
            simu = simu2;
            ii_count = ii_count + 1;
        else, break
        end
    end
    ii_count = 0;
    r2s_var(ii_var) = corr(y',g_shape_quick([],random_var_P(:,ii_var),x,[])' );
end
[~,best_idx] = max(r2s_var);

P = random_var_P(:,best_idx);
y_hat = g_shape_quick([],P,x,[]);

return

end


