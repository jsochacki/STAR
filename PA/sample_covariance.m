function [sample_covariance_matrix] = sample_covariance(observations)

[K N] = size(observations);
sample_mean_vector = sum(observations,2) / N;
A = observations - sample_mean_vector;
sample_covariance_matrix = (A * A.') / (N - 1);
end