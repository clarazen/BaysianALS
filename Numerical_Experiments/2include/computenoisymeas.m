function [Y,sigmae] = computenoisymeas(truth,I,SNR)

    %% Vectorization of corrupted tensor 
    noise_norm = norm(truth)/(10^(SNR/20));
    sigmae     = noise_norm^2/(numel(truth)-1);
    Y     = [];
    for i = 1:I
        noise      = randn(size(truth));
        noise      = noise/norm(noise)*noise_norm;
        y          = truth+noise;
        Y          = [Y y];
    end
end