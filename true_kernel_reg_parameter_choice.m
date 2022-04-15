% Add packages and folder paths
addpath('egrssMatlab')
input_folder = 'sharp_sample_images';

% =============== Load and blur sharp image ==============

images = {'barbara.tif' 'boat.tif' 'cameraman.tif' 'peppers.tif'}; 
im = 1;         % index for choosing between the four above
x_true = im2double(imread([input_folder '\' images{im}]));

% Options
use_chol = 1;       % Take uncertainty into account in x_update
Sx = 100;

R = 1:4;
ETA = [1,5,10,20];
lambdas = logspace(-3,2,30);

for j = 3%:length(R)
    for k = 2%:length(ETA)

        % blur the image
        r_true = R(j);                                    % radius of PSF
        kernel = 'disk';                                  % also try: 'gaussian', 'motion'
        PSF = fspecial(kernel, r_true);
        p = (size(PSF, 1) - 1) / 2;
        x_trueBC = padarray(x_true,[p p], 'symmetric');   % incorporate BC
        b_blurred = conv2(x_trueBC, PSF, 'valid');        % perform convolution

        % Add percentage relative Gaussian noise with noise level delta
        eta = ETA(k);                        % noise level in percentage
        d = randn(size(b_blurred));          % generate i.i.d. noise from N(0,1)
        d1 = d/norm(d)*norm(b_blurred);      % scale it to match the norm of data
        b = b_blurred + (eta/100)*d1;        % scale with percentage noise and add to data

        % blurred_data = figure(1); 
        % subplot(2,2,1); imshow(x_true); title("Ground truth image"); 
        % subplot(2,2,2); imagesc(PSF); title("Point spread function"); axis square;
        % subplot(2,2,3); imshow(b_blurred); title("Blurred image"); 
        % subplot(2,2,4); imshow(b); title("Blurred and noisy image"); 
        % drawnow;
        
        % Other Parameters
        delta_r0 = r_true/5;         % Initial variance 

        % Estimate noise standard deviation
        % sigma_e = std2(b(1:50,1:50)); % estimate noise std from small corner patch
        % use true noise standard deviation
        sigma_e = 1/norm(d)*norm(b_blurred)*eta/100;    
        % V[k*X] = k^2*V[X]  <=>  std[k*X] = k*std[X]
        delta_r = delta_r0;
        x0 = zeros(size(b));

        rel_err_opt = zeros(1,length(lambdas));

        for i = 1:length(lambdas)
            % ============ Deblur with true kernel ===========
            lambda_TV = lambdas(i);            % Regularization parameter for TV
            x_deblur_true = x_update(x0, r_true, delta_r, b, sigma_e, Sx, lambda_TV, use_chol);
            rel_err_opt(i) = norm(x_deblur_true-x_true)/norm(x_true);
        end


        figure(1)
        semilogx(lambdas,rel_err_opt,'-b','linewidth',2)
        ylabel('relative error','interpreter','latex','fontsize',14); 
        xlabel('$\lambda_{TV}$','interpreter','latex','fontsize',14); 
        title(['Deblurring with true kernel, $r_{true}=$ ' num2str(r_true) ' and noise = ' num2str(eta) '\%'],'interpreter','latex','fontsize',14); 
        drawnow; 
        %
        test_folder = ['C:\Users\Rainb\OneDrive\Documents\Skole\DTU\' ...
                '9. semester\HDC 2021 - Image deblurring project course\HDC paper\'...
                'reg_parameter_choice\r_true__' num2str(r_true) '\noise_level__' num2str(eta) '\'];
        out.r_true = r_true;
        out.lambdas = lambdas;
        out.kernel = kernel;
        out.true_noise_lvl = eta;
        out.alpha_var = alpha;
        out.delta_r= delta_r;
        out.relative_errors = rel_err_opt;
        [~,i_opt] = min(rel_err_opt);
        out.lambda_opt = lambdas(i_opt);
        save([test_folder images{im}(1:end-4) '_rel_error_vs_lambdaTV.mat'], '-struct', 'out')
    end
end

