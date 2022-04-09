clear; close all force; clc;

% Add packages and folder paths
addpath('egrssMatlab')
input_folder = 'sharp_sample_images';
% =============== Load and blur sharp image ==============

images = {'barbara.tif' 'boat.tif' 'cameraman.tif' 'peppers.tif'}; 
im = 1;         % index for choosing between the four above
x_true = im2double(imread([input_folder '\' images{im}]));

% blur the image
r_true = 3;                                       % radius of PSF
PSF = fspecial('disk', r_true);
p = (size(PSF, 1) - 1) / 2;
x_trueBC = padarray(x_true,[p p], 'symmetric');   % incorporate BC
b_blurred = conv2(x_trueBC, PSF, 'valid');        % perform convolution

% Add percentage relative Gaussian noise with noise level delta
eta = 1;                            % noise level in percentage
d = randn(size(b_blurred));         % generate i.i.d. noise from N(0,1)
d1 = d/norm(d)*norm(b_blurred);      % scale it to match the norm of data
b = b_blurred + (eta/100)*d1;        % scale with percentage noise and add to data

figure(1); 
subplot(2,2,1); imshow(x_true); title("Ground truth image"); 
subplot(2,2,2); imagesc(PSF); title("Point spread function"); axis square;
subplot(2,2,3); imshow(b_blurred); title("Blurred image"); 
subplot(2,2,4); imshow(b); title("Blurred and noisy image"); 
drawnow;

% =============== Algorithm initialization ===============

% Options
save_deblur = 0;    % Save output deblurred image?
save_x_test = 1;    % Save output deblurred image as .mat file (used in testing)
use_egrss = 1;      % Use egrss package for r_update? If 0 only works on small-scale.
use_gpu = 0;        % Use gpu for faster computations? Requires Parallel computing toolbox
use_chol = 1;       % Take uncertainty into account in x_update
estimate_r = 1;     % binary flag for using the iterative radius estimation 

n_iter = 5;         % Number of iterations for r_update
Sr = 100;           % Number of samples for r_update
Sx = 100;           % Number of samples for x_update
alpha = 0.5;        % Relaxation parameter in varience est for r_update

% patch parameters for r estimation
use_patch = 0;
patch_width  = 150; % Width of patch for radius estimation
patch_height = 150; % Height of patch for radius estimation

% Other Parameters
lambda_TV = 1e-3;            % Regularization parameter for TV
mu_r0 = r_true - r_true/5;   % Initial radius, set to 10% smaller/larger than true one
delta_r0 = r_true/5;              % Initial variance 

% pre-allocate iteration progression
progress = zeros(n_iter+1,4);
x_iterates = zeros(length(b), width(b), n_iter+1);

% We plots N(mu_r0, delta_r0^2) to see if r_true is within +-3*delta_r0
x_grid = linspace(0,2*r_true);
figure(4); hold on; 
plot(x_grid, normpdf(x_grid, mu_r0, delta_r0), 'linewidth', 2);
plot(r_true*ones(1,10),linspace(0,normpdf(r_true, mu_r0, delta_r0),10), '-r', 'linewidth', 2)
legend('$\mathcal{N}(\mu_{r_0},\delta_{r_0}^2)$', '$r_{true}$','interpreter', 'latex', 'fontsize', 14, 'location', 'best')

%% =============== Algorithm start ===============

% Estimate noise standard deviation
%sigma_e = std2(b(1:50,1:50)); % estimate noise std from small corner patch
% use true noise standard deviation
sigma_e = 1/norm(d)*norm(b_blurred)*eta/100;    
% V[k*X] = k^2*V[X]  <=>  std[k*X] = k*std[X]

% initial iteration
mu_r = mu_r0;
delta_r = delta_r0;
x0 = zeros(size(b));
[mu_r,delta_r,norm(x0-x_true)/norm(x_true),0]
progress(1,:) = [mu_r,delta_r,norm(x0-x_true)/norm(x_true),0];
x_iterates(:,:,1) = x0;

% ==== Iteration for r estimation =====
if estimate_r       
    figure(2); 
    imshow(b);    
    if use_patch
        % ==== Prepare patches =====
        mid = floor(size(b)/2);
        hpatch_width = floor(patch_width/2);
        hpatch_height = floor(patch_height/2);
        b_patch = b(mid(1)-hpatch_height:mid(1)+hpatch_height, mid(2)-hpatch_width:mid(2)+hpatch_width);
        x = x0(mid(1)-hpatch_height:mid(1)+hpatch_height, mid(2)-hpatch_width:mid(2)+hpatch_width);
        title("Blurred and noisy image with patch"); 
        rectangle('Position',[mid(2)-hpatch_width,mid(1)-hpatch_height,patch_width,patch_height],'linewidth',3,'edgecolor','red')
    else
        b_patch = b;
        title("Blurred and noisy image"); 
        x = x0;
    end    
    drawnow; 
    % iterate
    for k = 1:n_iter
        % Update x
        x = x_update(x, mu_r, delta_r, b_patch, sigma_e, Sx, lambda_TV, use_chol);
        % save initial reconstruction explicitly
        if k == 1; initial_deblur = x; end 
        figure(2); imshow(x); title('Current deblurred image'); drawnow;

        % Update r
        %[mu_r, delta_r] = r_update(x, b_patch, mu_r, delta_r, sigma_e, Sr, alpha, use_egrss);
        [mu_r, delta_r] = r_update(x_true, b_patch, mu_r, delta_r, sigma_e, Sr, alpha, use_egrss);

        % Show result and save progress
        [mu_r,delta_r, norm(x-x_true)/norm(x_true), k]
        progress(k+1,:) = [mu_r,delta_r, norm(x-x_true)/norm(x_true), k];
        x_iterates(:,:,k+1) = x;
    end
end

if use_gpu == 1
    x = gpuArray(zeros(size(b)));
    b = gpuArray(b);
else
    x = zeros(size(b));
end

% ==== Deblur with radius estimate ====
x_estimate = x_update(x, mu_r, delta_r, b, sigma_e, 0, lambda_TV, use_chol);
if use_gpu == 1
    x_estimate = gather(x_estimate);
end

figure(3); 
subplot(1,2,1); imshow(b); title("Blurred and noisy image"); drawnow;
subplot(1,2,2); imshow(x_estimate); title("Deblurred image"); 

% Save to file
if save_deblur == 1
    % saves image to output folder
    test_folder = ['C:\Users\Rainow Slayer\OneDrive\Documents\Skole\DTU\' ...
        '9. semester\HDC 2021 - Image deblurring project course\HDC paper\'...
        'tests\r_true__' num2str(r_true) '\noise_level__' num2str(eta) '\lambdaTV__' num2str(lambda_TV)];
    % save to file in specific folder
    output_file = [output_folder '/' images{im}(1:end-4) '_r_' num2str(r_true)...
        '_nl_' num2str(eta) '_lambdaTV_' num2str(lambda_TV) '.png'];
    imwrite(x_estimate, output_file)
end
if save_x_test == 1
    % save parameters used for test in a struct
    test_params.true_im = x_true;
    test_params.blurred_im = b;
    test_params.r_true = r_true;
    test_params.rel_noise_lvl = eta;
    test_params.alpha_var = alpha;
    test_params.lambda = lambda_TV;
    test_params.r_init = mu_r0;
    test_params.dr_init = delta_r0;
    test_params.r_final = mu_r;
    test_params.dr_final = delta_r;
    test_params.initial_deblur = initial_deblur;
    test_params.deblurred_im = x_estimate;
    test_params.image_iterates = x_iterates;
    test_params.iterations = progress;
    % create a folder in the test folder system for this value of lambda_TV
    % (you can create the folder system by running the folder_generator.m 
    % file, just remember to change the folder paths to your own below and
    % in the folder_generator.m file)
    test_folder = ['C:\Users\Rainow Slayer\OneDrive\Documents\Skole\DTU\' ...
        '9. semester\HDC 2021 - Image deblurring project course\HDC paper\'...
        'tests\r_true__' num2str(r_true) '\noise_level__' num2str(eta) '\'];
    folder = [test_folder 'lambdaTV__' num2str(lambda_TV)];
    status = mkdir(folder);
    % save to file in specific folder
    filename = [folder '\' images{im}(1:end-4) '_test.mat'];
    save(filename, '-struct', 'test_params')
end