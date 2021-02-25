barbara_img = double(imread('../HW2/barbara256.png')) + 2*randn(256,256);
lambda = 1;
U = kron(haarmtx(8).', haarmtx(8).');

H = haarmtx(8);

patch_size = 8;
num_rows = 256;
num_columns = 256;
reconstructed_frames = zeros(num_rows, num_columns);
all_weights = zeros(num_rows, num_columns);

phi = randn(32,64);

A = phi*U;
alpha = eigs(double(A.' * A), 1) + 1;

for patch_start_row = 1:(num_rows-patch_size+1) % Iterating over all possible patches
    for patch_start_column = 1:(num_columns-patch_size+1)
        current_patch = barbara_img(patch_start_row:patch_start_row+patch_size-1,patch_start_column:patch_start_column+patch_size-1);
        
        y = phi * double(reshape(current_patch.',[1, 64]).');
        theta = 255*rand(64, 1);
        
        % Implementing ISTA
        iter=1;
        norm_residual = 801;
        while norm_residual>800 && iter<100
            theta = wthresh(theta + (1/alpha) * A.' * (y - A*theta), 's', lambda/(2*alpha));
            residual = y - A*theta;
            norm_residual = norm(residual);
            iter = iter + 1;
        end
                
        x =  H.' * reshape(theta, [8,8]).' * H;
        
        % Rejoining the patches
        weights = all_weights(patch_start_row:patch_start_row+patch_size-1,patch_start_column:patch_start_column+patch_size-1);

        reconstructed_frames(patch_start_row:patch_start_row+patch_size-1,patch_start_column:patch_start_column+patch_size-1) ...
         = (reconstructed_frames(patch_start_row:patch_start_row+patch_size-1,patch_start_column:patch_start_column+patch_size-1) .* (weights./(weights+1))) + (x ./ (weights+1));

        all_weights(patch_start_row:patch_start_row+patch_size-1,patch_start_column:patch_start_column+patch_size-1) ...
            = (weights+1);
       
    end
end
imshow(uint8(reconstructed_frames))
original_img = imread('../HW2/barbara256.png');
% Displaying RMSE
sqrt(sum((uint8(reconstructed_frames) - uint8(original_img)).^2))/sqrt(sum(uint8(original_img).^2))
