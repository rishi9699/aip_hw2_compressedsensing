barbara_img = imread('../HW2/barbara256.png');
lambda = 1;
U = kron(dctmtx(8).', dctmtx(8).');

patch_size = 8;
num_rows = 256;
num_columns = 256;
reconstructed_frames = zeros(num_rows, num_columns);
all_weights = zeros(num_rows, num_columns);

for patch_start_row = 1:(num_rows-patch_size+1) % Iterating over all possible patches
    for patch_start_column = 1:(num_columns-patch_size+1)
        current_patch = barbara_img(patch_start_row:patch_start_row+patch_size-1,patch_start_column:patch_start_column+patch_size-1);
        phi = randn(32,64);
        y = phi * double(reshape(current_patch.',[1, 64]).');
        A = phi*U;
        alpha = eigs(double(A.' * A), 1) + 1;
        
        theta = 255*rand(64, 1);
        for iter=1:100
            theta = wthresh(theta + (1/alpha) * A.' * (y - A*theta), 's', lambda/(2*alpha));
        end
        
        % Converting DCT coeffs to image patch
        x = idct2( reshape(theta, [8,8]).');
        
        % Rejoining the patches
        weights = all_weights(patch_start_row:patch_start_row+patch_size-1,patch_start_column:patch_start_column+patch_size-1);

        reconstructed_frames(patch_start_row:patch_start_row+patch_size-1,patch_start_column:patch_start_column+patch_size-1) ...
         = (reconstructed_frames(patch_start_row:patch_start_row+patch_size-1,patch_start_column:patch_start_column+patch_size-1) .* (weights./(weights+1))) + (x ./ (weights+1));

        all_weights(patch_start_row:patch_start_row+patch_size-1,patch_start_column:patch_start_column+patch_size-1) ...
            = (weights+1);
    end
end
