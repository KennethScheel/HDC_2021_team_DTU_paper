
fold1 = 'C:\Users\Rainow Slayer\OneDrive\Documents\Skole\DTU\9. semester\HDC 2021 - Image deblurring project course\HDC paper\tests\';
mkdir(fold1)

nR = 1:4;       
nL = 0:0.5:5;

for i = 1:length(nR)
    fold2 = [fold1 'r_true__' num2str(nR(i))];
    status = mkdir(fold2)
    for j = 1:length(nL)
        fold3 = [fold2 '\' 'noise_level__' num2str(nL(j))];
        status = mkdir(fold3)
    end
end