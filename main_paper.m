clear; close all force; clc;

% Add packages and folder paths
addpath('egrssMatlab')
input_folder = 'sharp_sample_images';

% =============== Load and blur sharp image ==============

images = {'barbara.tif' 'boat.tif' 'cameraman.tif' 'peppers.tif'}; 
im = 1;         % index for choosing between the four above
x_true = im2double(imread([input_folder '\' images{im}]));
%x_true = imresize(x_true, 0.5);                  % uncomment if we do not use egrss package


% blur the image
r_true = 3;                                       % radius of PSF
kernel = 'disk';                                  % also try: 'gaussian', 'motion'
PSF = fspecial(kernel, r_true);
p = (size(PSF, 1) - 1) / 2;
x_trueBC = padarray(x_true,[p p], 'symmetric');   % incorporate BC
b_blurred = conv2(x_trueBC, PSF, 'valid');        % perform convolution

% Add percentage relative Gaussian noise with noise level delta
eta = 5;                             % noise level in percentage
d = randn(size(b_blurred));          % generate i.i.d. noise from N(0,1)
d1 = d/norm(d)*norm(b_blurred);      % scale it to match the norm of data
b = b_blurred + (eta/100)*d1;        % scale with percentage noise and add to data

blurred_data = figure(1); 
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
lambda_TV = 5;            % Regularization parameter for TV
mu_r0 = r_true - r_true/5;   % Initial radius, set to 10% smaller/larger than true one
delta_r0 = r_true/5;         % Initial variance 

% pre-allocate iteration progression
progress = zeros(n_iter+1,5);
x_iterates = zeros(length(b), width(b), n_iter+1);

% We plot N(mu_r0, delta_r0^2) to see if r_true is within +-3*delta_r0
x_grid = linspace(0,2*r_true);
figure(2); hold on; 
plot(x_grid, normpdf(x_grid, mu_r0, delta_r0), 'linewidth', 2);
plot(r_true*ones(1,10),linspace(0,normpdf(r_true, mu_r0, delta_r0),10), '-r', 'linewidth', 2)
legend('$\mathcal{N}(\mu_{r_0},\delta_{r_0}^2)$', '$r_{true}$','interpreter', 'latex', 'fontsize', 14, 'location', 'best')

% Estimate noise standard deviation
%sigma_e = std2(b(1:50,1:50)); % estimate noise std from small corner patch
% use true noise standard deviation
sigma_e = 1/norm(d)*norm(b_blurred)*eta/100;    
% V[k*X] = k^2*V[X]  <=>  std[k*X] = k*std[X]

% initial iteration
mu_r = mu_r0;
delta_r = delta_r0;
x0 = zeros(size(b));
[mu_r,delta_r,abs(mu_r-r_true),0]
progress(1,:) = [mu_r,delta_r,abs(mu_r-r_true),norm(x0-x_true)/norm(x_true),0];
x_iterates(:,:,1) = x0;

% ============ Deblur with true kernel ===========
x_deblur_true = x_update(x0, r_true, delta_r, b, sigma_e, Sx, lambda_TV, use_chol);
rel_err_opt = norm(x_deblur_true-x_true)/norm(x_true);
% optimal result we can expect
true_kernel_deblur_fig = figure(3); 
subplot(1,2,1); imshow(b); title("Blurred and noisy image"); drawnow;
subplot(1,2,2); imshow(x_deblur_true); 
title(['Deblur w. true kernel, rel\_err = ' num2str(rel_err_opt)]); drawnow; 

%% =============== Algorithm start ===============
if estimate_r       
    if use_patch
        % ==== Prepare patches =====
        mid = floor(size(b)/2);
        hpatch_width = floor(patch_width/2);
        hpatch_height = floor(patch_height/2);
        b_patch = b(mid(1)-hpatch_height:mid(1)+hpatch_height, mid(2)-hpatch_width:mid(2)+hpatch_width);
        x = x0(mid(1)-hpatch_height:mid(1)+hpatch_height, mid(2)-hpatch_width:mid(2)+hpatch_width);
        figure(4); imshow(b); title("Blurred and noisy image with patch"); 
        rectangle('Position',[mid(2)-hpatch_width,mid(1)-hpatch_height,patch_width,patch_height],'linewidth',3,'edgecolor','red')
        drawnow; 
    else
        b_patch = b;
        x = x0;
    end    
    % iterate
    for k = 1:n_iter
        % Update x
        x = x_update(x, mu_r, delta_r, b_patch, sigma_e, Sx, lambda_TV, use_chol);
        % save initial reconstruction explicitly
        if k == 1 
            initial_deblur = x; 
            initial_deblur_fig = figure(5); imshow(initial_deblur); title('Initial deblurred image'); drawnow;
        else  
            figure(6); imshow(x); title('Current deblurred image'); drawnow;
        end
        
        % Update r
        %[mu_r, delta_r] = r_update(x, b_patch, mu_r, delta_r, sigma_e, Sr, alpha, use_egrss);
        [mu_r, delta_r] = r_update(x_true, b_patch, mu_r, delta_r, sigma_e, Sr, alpha, use_egrss);

        % Show result and save progress
        [mu_r,delta_r, abs(mu_r-r_true), k]
        progress(k+1,:) = [mu_r,delta_r, abs(mu_r-r_true),norm(x-x_true)/norm(x_true), k];
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
x_estimate = x_update(x, mu_r, delta_r, b, sigma_e, Sx, lambda_TV, use_chol);
rel_err_final = norm(x_estimate-x_true)/norm(x_true);

if use_gpu == 1
    x_estimate = gather(x_estimate);
end

%%
final_deblur_fig = figure(7); 
subplot(1,2,1); imshow(b); title("Blurred and noisy image"); drawnow;
subplot(1,2,2); imshow(x_estimate); title(['Deblur w. est. kernel, rel\_err = ' num2str(rel_err_final)]); drawnow; 

%%
r_conv_fig = figure(8);
hold on
plot(0:k, r_true*ones(1,k+1), '--k', 'linewidth',2)
plot(0:k, progress(:,1), '-b', 'linewidth',2)
patch([0:k, fliplr(0:k)], [(progress(:,1)+progress(:,2))', fliplr((progress(:,1)-progress(:,2))')],'b','EdgeColor','none','FaceAlpha',0.5)
plot(0:k, progress(:,3), '-r', 'linewidth',2)
xlabel('$k$','fontsize',14,'interpreter','latex')
ylabel('$\mu_r$','fontsize',14,'interpreter','latex')
ylim([0,r_true+0.5])
hold off
legend('$r_{true}$','$\mu_r$','$\mu_r \pm \delta_r$','$|\mu_r-r_{true}|$','fontsize',14,'interpreter','latex','location','best')

%% Save test data to file
if save_x_test == 1
    % save parameters used for test in a struct
    test_params.true_im = x_true;
    test_params.blurred_im = b;
    test_params.initial_deblur = initial_deblur;
    test_params.opt_deblur = x_deblur_true;
    test_params.opt_rel_err = rel_err_opt;
    test_params.final_deblur = x_estimate;
    test_params.final_rel_err = rel_err_final;
    test_params.r_update_image_iter = x_iterates;
    T = table(progress(:,5),progress(:,1),progress(:,2), progress(:,3), progress(:,4),...
          'VariableNames',{'k','mu_r','delta_r','|r-r_true|','||x-x_true||/||x||'});
    test_params.r_update_iter_info = T;
    test_params.kernel = kernel;
    test_params.r_true = r_true;
    test_params.true_noise_lvl = eta;
    test_params.alpha_var = alpha;
    test_params.lambdaTV = lambda_TV;
    test_params.r_init = mu_r0;
    test_params.dr_init = delta_r0;
    test_params.r_final = mu_r;
    test_params.dr_final = delta_r;
    % create a folder in the test folder system for this value of lambda_TV
    % (you can create the folder system by running the folder_generator.m 
    % file, just remember to change the folder paths to your own below and
    % in the folder_generator.m file)
    test_folder = ['C:\Users\Rainow Slayer\OneDrive\Documents\Skole\DTU\' ...
        '9. semester\HDC 2021 - Image deblurring project course\HDC paper\'...
        'tests_2\r_true__' num2str(r_true) '\noise_level__' num2str(eta) '\'];
    folder = [test_folder 'lambdaTV__' num2str(lambda_TV)];
    status = mkdir(folder);
    % save struct to file in specific folder
    filename = [folder '\' images{im}(1:end-4) '_test.mat'];
    save(filename, '-struct', 'test_params')

    % save images to files
    final_deblur_file = [folder '/' images{im}(1:end-4) '_final_deblur' '.png'];
    blurred_data_file = [folder '/' images{im}(1:end-4) '_blurred_data' '.png'];
    true_deblur_file = [folder '/' images{im}(1:end-4) '_true_kernel_deblur' '.png'];
    r_convergence_file = [folder '/' images{im}(1:end-4) '_kernel_radius_convergence' '.png'];
    init_deblur_file = [folder '/' images{im}(1:end-4) '_initial_deblur' '.png'];
    saveas(final_deblur_fig, final_deblur_file);
    saveas(blurred_data, blurred_data_file);
    saveas(true_kernel_deblur_fig, true_deblur_file);
    saveas(r_conv_fig, r_convergence_file);
    saveas(initial_deblur_fig, init_deblur_file);

end
