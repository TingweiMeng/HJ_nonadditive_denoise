function generate_data(case_num, t)
% load original image and add Poisson noise on it
% it will automatically save the original image x_ori, the noisy image
% x_noisy, the parameter t and the psnr in a mat file
%
% % INPUT: 
%   case_num (which image to load): 
%               1 for 4box, 2 for wfc3_uvis_full_field, 3 for abell_2744
%   t:          positive scalar, parameter in the noise (see our paper)
%
if case_num == 1 % 4box
    x_ori = zeros(256,256) + 0.2;
    x_ori(50:206, 50:206) = 0.4;
    x_ori(80:176, 80:176) = 0.6;
    x_ori(110:146, 110:146) = 0.8;
elseif case_num == 2 % wfc3_uvis_full_field
    imageData = imread('wfc3_uvis_full_field.jpg');
elseif case_num == 3 % abell_2744
    imageData = imread('Abell2744.jpg');
else
    fprintf('wrong case number\n');
end

if case_num > 1 % convert images to gray level and range [0,1]
    if size(imageData,3) == 3
        imageData = rgb2gray(imageData);
    end
    x_ori = im2double(imageData);
end
figure; imshow(x_ori); title('original image');
%% add noise
x_noisy = generate_poisson_noise(x_ori, t);
figure; imshow(x_noisy); title('noisy image');

psnr = 10* log(max(x_ori(:))^2 / mean((x_ori(:) - x_noisy(:)).^2))/ log(10);

folder_name = sprintf('./case%d',case_num);
if ~exist(folder_name, 'dir')
   mkdir(folder_name)
end
save(sprintf('%s/data.mat', folder_name), 't', 'x_ori', 'x_noisy', 'psnr');
end