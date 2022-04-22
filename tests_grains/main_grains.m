%% ============== Experiment 1.1 ====================
% deblurring with true radius with uncertainty in x_update

clear; close all force; clc;

% Add packages and folder paths
addpath('egrssMatlab')
input_folder = 'sharp_sample_images';

% Load and blur sharp image
image = 'grains.tif'; 
x_true = im2double(imread([input_folder '\' image]));

r_true = 3;                                       % radius of PSF
kernel = 'disk';                                  % also try: 'gaussian', 'motion'
PSF = fspecial(kernel, r_true);
p = (size(PSF, 1) - 1) / 2;
x_trueBC = padarray(x_true,[p p], 'symmetric');   % incorporate BC
b_blurred = conv2(x_trueBC, PSF, 'valid');        % perform convolution

% Add percentage relative Gaussian noise with noise level delta
eta = 1;                             % noise level in percentage
d = randn(size(b_blurred));          % generate i.i.d. noise from N(0,1)
d1 = d/norm(d)*norm(b_blurred);      % scale it to match the norm of data
b = b_blurred + (eta/100)*d1;        % scale with percentage noise and add to data

% Other Parameters
lambda_TV = 10;             % Regularization parameter for TV
Sx = 250;                   % Number of samples in x_update
delta_r = 0.3;              % Initial radius variance 
use_chol = 1;               % Include uncertainty in model

% use true noise standard deviation
sigma_e = 1/norm(d)*norm(b_blurred)*eta/100;    % V[k*X] = k^2*V[X]  <=>  std[k*X] = k*std[X]

% Deblur with true kernel and with taking uncertainty into account
x0 = zeros(size(b));
x_deblur = x_update(x0, r_true, delta_r, b, sigma_e, Sx, lambda_TV, use_chol);
rel_err = norm(x_deblur-x_true)/norm(x_true);

figure(1); 
subplot(2,2,1); imshow(x_true); title("Ground truth image"); 
subplot(2,2,2); imagesc(PSF); title("Point spread function"); axis square;
subplot(2,2,3); imshow(b); title("Blurred and noisy image"); 
subplot(2,2,4); imshow(x_deblur); 
title(['w. uncertainty, true r = ' num2str(r_true) ', rel\_err = ' num2str(rel_err)]); drawnow; 
drawnow;

% Save workspace to file
test_folder = ['C:\Users\Rainb\OneDrive\Documents\Skole\DTU\' ...
    '9. semester\HDC 2021 - Image deblurring project course\HDC paper\'...
    'HDC_2021_team_DTU_paper\tests_grains\deblur_with_true\with_uncertainty_xupdate\'];
save([test_folder 'test_workspace.mat'])

%% ============== Experiment 1.2 ====================
% deblurring with true radius without uncertainty in x_update

clear; close all force; clc;

% Add packages and folder paths
addpath('egrssMatlab')
input_folder = 'sharp_sample_images';

% Load and blur sharp image
image = 'grains.tif'; 
x_true = im2double(imread([input_folder '\' image]));

r_true = 3;                                       % radius of PSF
kernel = 'disk';                                  % also try: 'gaussian', 'motion'
PSF = fspecial(kernel, r_true);
p = (size(PSF, 1) - 1) / 2;
x_trueBC = padarray(x_true,[p p], 'symmetric');   % incorporate BC
b_blurred = conv2(x_trueBC, PSF, 'valid');        % perform convolution

% Add percentage relative Gaussian noise with noise level delta
eta = 1;                             % noise level in percentage
d = randn(size(b_blurred));          % generate i.i.d. noise from N(0,1)
d1 = d/norm(d)*norm(b_blurred);      % scale it to match the norm of data
b = b_blurred + (eta/100)*d1;        % scale with percentage noise and add to data

% Other Parameters
lambda_TV = 1e-2;             % Regularization parameter for TV
Sx = 250;                   % Number of samples in x_update
delta_r = 0.3;              % Initial radius variance 
use_chol = 0;               % Include uncertainty in model

% use true noise standard deviation
sigma_e = 1/norm(d)*norm(b_blurred)*eta/100;    % V[k*X] = k^2*V[X]  <=>  std[k*X] = k*std[X]

% Deblur with true kernel and without taking uncertainty into account
x0 = zeros(size(b));
x_deblur = x_update(x0, r_true, delta_r, b, sigma_e, Sx, lambda_TV, use_chol);
rel_err = norm(x_deblur-x_true)/norm(x_true);

figure(1); 
subplot(2,2,1); imshow(x_true); title("Ground truth image"); 
subplot(2,2,2); imagesc(PSF); title("Point spread function"); axis square;
subplot(2,2,3); imshow(b); title("Blurred and noisy image"); 
subplot(2,2,4); imshow(x_deblur); 
title(['wo. uncertainty, r = ' num2str(r_true) ', rel\_err = ' num2str(rel_err)]); drawnow; 
drawnow;

% Save workspace to file
test_folder = ['C:\Users\Rainb\OneDrive\Documents\Skole\DTU\' ...
    '9. semester\HDC 2021 - Image deblurring project course\HDC paper\'...
    'HDC_2021_team_DTU_paper\tests_grains\deblur_with_true\without_uncertainty_xupdate\'];
save([test_folder 'test_workspace.mat'])


%% ============== Experiment 2.1 ====================
% deblurring with a wrong radius +- 0.5 with uncertainty in x_update

clear; close all force; clc;

% Add packages and folder paths
addpath('egrssMatlab')
input_folder = 'sharp_sample_images';

% Load and blur sharp image
image = 'grains.tif'; 
x_true = im2double(imread([input_folder '\' image]));

r_true = 3;                                       % radius of PSF
kernel = 'disk';                                  % also try: 'gaussian', 'motion'
PSF = fspecial(kernel, r_true);
p = (size(PSF, 1) - 1) / 2;
x_trueBC = padarray(x_true,[p p], 'symmetric');   % incorporate BC
b_blurred = conv2(x_trueBC, PSF, 'valid');        % perform convolution

% Add percentage relative Gaussian noise with noise level delta
eta = 1;                             % noise level in percentage
d = randn(size(b_blurred));          % generate i.i.d. noise from N(0,1)
d1 = d/norm(d)*norm(b_blurred);      % scale it to match the norm of data
b = b_blurred + (eta/100)*d1;        % scale with percentage noise and add to data

% Other Parameters
lambda_TV = 10;             % Regularization parameter for TV
Sx = 250;                   % Number of samples in x_update
mu_r = r_true + 0.5;        % Initial radius estimate
delta_r = 0.3;              % Initial radius variance 
use_chol = 1;               % Include uncertainty in model

% use true noise standard deviation
sigma_e = 1/norm(d)*norm(b_blurred)*eta/100;    % V[k*X] = k^2*V[X]  <=>  std[k*X] = k*std[X]

% Deblur with wrong kernel and with taking uncertainty into account
x0 = zeros(size(b));
x_deblur = x_update(x0, mu_r, delta_r, b, sigma_e, Sx, lambda_TV, use_chol);
rel_err = norm(x_deblur-x_true)/norm(x_true);

figure(1); 
subplot(2,2,1); imshow(x_true); title("Ground truth image"); 
subplot(2,2,2); imagesc(PSF); title("Point spread function"); axis square;
subplot(2,2,3); imshow(b); title("Blurred and noisy image"); 
subplot(2,2,4); imshow(x_deblur); 
title(['w. uncertainty, r = ' num2str(mu_r) ', rel\_err = ' num2str(rel_err)]); drawnow; 
drawnow;

% Save workspace to file
test_folder = ['C:\Users\Rainb\OneDrive\Documents\Skole\DTU\' ...
    '9. semester\HDC 2021 - Image deblurring project course\HDC paper\'...
    'HDC_2021_team_DTU_paper\tests_grains\deblur_with_wrong\with_uncertainty_xupdate\'];
save([test_folder 'test_workspace.mat'])

%% ============== Experiment 2.2 ====================
% deblurring with a wrong radius +- 0.5 without uncertainty in x_update

clear; close all force; clc;

% Add packages and folder paths
addpath('egrssMatlab')
input_folder = 'sharp_sample_images';

% Load and blur sharp image
image = 'grains.tif'; 
x_true = im2double(imread([input_folder '\' image]));

r_true = 3;                                       % radius of PSF
kernel = 'disk';                                  % also try: 'gaussian', 'motion'
PSF = fspecial(kernel, r_true);
p = (size(PSF, 1) - 1) / 2;
x_trueBC = padarray(x_true,[p p], 'symmetric');   % incorporate BC
b_blurred = conv2(x_trueBC, PSF, 'valid');        % perform convolution

% Add percentage relative Gaussian noise with noise level delta
eta = 1;                             % noise level in percentage
d = randn(size(b_blurred));          % generate i.i.d. noise from N(0,1)
d1 = d/norm(d)*norm(b_blurred);      % scale it to match the norm of data
b = b_blurred + (eta/100)*d1;        % scale with percentage noise and add to data

% Other Parameters
lambda_TV = 1e-2;             % Regularization parameter for TV
Sx = 250;                   % Number of samples in x_update
mu_r = r_true - 0.5;        % Initial radius estimate
delta_r = 0.3;              % Initial radius variance 
use_chol = 0;               % Include uncertainty in model

% use true noise standard deviation
sigma_e = 1/norm(d)*norm(b_blurred)*eta/100;    % V[k*X] = k^2*V[X]  <=>  std[k*X] = k*std[X]

% Deblur with wrong kernel and without taking uncertainty into account
x0 = zeros(size(b));
x_deblur = x_update(x0, mu_r, delta_r, b, sigma_e, Sx, lambda_TV, use_chol);
rel_err = norm(x_deblur - x_true)/norm(x_true);

figure(1); 
subplot(2,2,1); imshow(x_true); title("Ground truth image"); 
subplot(2,2,2); imagesc(PSF); title("Point spread function"); axis square;
subplot(2,2,3); imshow(b); title("Blurred and noisy image"); 
subplot(2,2,4); imshow(x_deblur); 
title(['wo. uncertainty, r = ' num2str(mu_r) ', rel\_err = ' num2str(rel_err)]); drawnow; 
drawnow;

% Save workspace to file
test_folder = ['C:\Users\Rainb\OneDrive\Documents\Skole\DTU\' ...
    '9. semester\HDC 2021 - Image deblurring project course\HDC paper\'...
    'HDC_2021_team_DTU_paper\tests_grains\deblur_with_wrong\without_uncertainty_xupdate\'];
save([test_folder 'test_workspace.mat'])


%% ============== Experiment 3.1 ====================
% Estimating radius via r_update and deblurring without uncertainty

clear; close all force; clc;

% Add packages and folder paths
addpath('egrssMatlab')
input_folder = 'sharp_sample_images';

% Load and blur sharp image
image = 'grains.tif'; 
x_true = im2double(imread([input_folder '\' image]));

r_true = 3;                                       % radius of PSF
kernel = 'disk';                                  % also try: 'gaussian', 'motion'
PSF = fspecial(kernel, r_true);
p = (size(PSF, 1) - 1) / 2;
x_trueBC = padarray(x_true,[p p], 'symmetric');   % incorporate BC
b_blurred = conv2(x_trueBC, PSF, 'valid');        % perform convolution

% Add percentage relative Gaussian noise with noise level delta
eta = 1;                             % noise level in percentage
d = randn(size(b_blurred));          % generate i.i.d. noise from N(0,1)
d1 = d/norm(d)*norm(b_blurred);      % scale it to match the norm of data
b = b_blurred + (eta/100)*d1;        % scale with percentage noise and add to data

% Other Parameters
lambda_TV = 1e-2;           % Regularization parameter for TV
Sx = 250;                   % Number of samples in x_update
Sr = 250;                   % Number of samples in r_update
mu_r = r_true + 0.5;        % Initial radius estimate
delta_r = 0.3;              % Initial radius variance 
use_chol = 0;               % Include uncertainty in model
use_egrss = 0;              % Use egrss package for faster computations
alpha = 0.5;                % Relaxation parameter for variance update
n_iter = 5;                % Number of iterations for r estimate
progress = zeros(n_iter+1,5);

% use true noise standard deviation
sigma_e = 1/norm(d)*norm(b_blurred)*eta/100;    % V[k*X] = k^2*V[X]  <=>  std[k*X] = k*std[X]

x0 = zeros(size(b));
x = x0;
[mu_r,delta_r, abs(mu_r-r_true), 0]
progress(1,:) = [mu_r,delta_r,abs(mu_r-r_true),norm(x0-x_true)/norm(x_true),0];

% iterate
for k = 1:n_iter
    % Update x
    x = x_update(x, mu_r, delta_r, b, sigma_e, Sx, lambda_TV, use_chol);
    rel_err_k = norm(x-x_true)/norm(x_true);

    figure(1)
    if k == 1 
        initial_deblur = x; 
        init_r = mu_r;
        init_rel_err = norm(initial_deblur-x_true)/norm(x_true);
    end
    subplot(1,2,1)
    imshow(initial_deblur);
    title(['initial, est. r = ' num2str(init_r) ', rel\_err = ' num2str(init_rel_err)]);
    subplot(1,2,2)
    imshow(x); 
    title(['iter ' num2str(k) ', est. r = ' num2str(mu_r) ', rel\_err = ' num2str(rel_err_k)]);
    drawnow;
    
    % Update r
    %[mu_r, delta_r] = r_update(x, b, mu_r, delta_r, sigma_e, Sr, alpha, use_egrss);
    [mu_r, delta_r] = r_update(x_true, b, mu_r, delta_r, sigma_e, Sr, alpha, use_egrss);

    % Show result and save progress
    [mu_r,delta_r, abs(mu_r-r_true), k]
    progress(k+1,:) = [mu_r,delta_r, abs(mu_r-r_true),norm(x-x_true)/norm(x_true), k];
end

x = zeros(size(b));

% Deblur with final radius estimate 
x_estimate = x_update(x, mu_r, delta_r, b, sigma_e, Sx, lambda_TV, use_chol);
rel_err_final = norm(x_estimate-x_true)/norm(x_true);

figure(2); 
subplot(2,2,1); imshow(x_true); title("Ground truth image"); 
subplot(2,2,2); imagesc(PSF); title("Point spread function"); axis square;
subplot(2,2,3); imshow(b); title("Blurred and noisy image"); 
subplot(2,2,4); imshow(x_estimate); 
title(['wo. uncertainty, est. r = ' num2str(mu_r) ', rel\_err = ' num2str(rel_err_final)]); drawnow; 
drawnow;

figure(3);
hold on
plot(0:k, r_true*ones(1,k+1), '--k', 'linewidth',2)
plot(0:k, progress(:,1), '-b', 'linewidth',2)
patch([0:k, fliplr(0:k)], [(progress(:,1)+progress(:,2))', fliplr((progress(:,1)-progress(:,2))')],'b','EdgeColor','none','FaceAlpha',0.5)
plot(0:k, progress(:,3), '-r', 'linewidth',2)
xlabel('$k$','fontsize',14,'interpreter','latex')
ylabel('$\mu_r$','fontsize',14,'interpreter','latex')
ylim([0,mu_r+0.5])
hold off
legend('$r_{true}$','$\mu_r$','$\mu_r \pm \delta_r$','$|\mu_r-r_{true}|$','fontsize',14,'interpreter','latex','location','best')

% Save workspace to file
test_folder = ['C:\Users\Rainb\OneDrive\Documents\Skole\DTU\' ...
    '9. semester\HDC 2021 - Image deblurring project course\HDC paper\'...
    'HDC_2021_team_DTU_paper\tests_grains\estimate_without_uncertainty\'];
%save([test_folder 'test_workspace.mat'])

%% ============== Experiment 3.2 ====================
% Estimating radius via r_update and deblurring with uncertainty.
% W.I.P.