clearvars;
rows_per_proc = 32; % Number of rows per processor
n_procs = 2; % Length of processor grid
N = n_procs * rows_per_proc; % Global number of rows
rng(0);

sparse_flag = 1;
nnz_per_row = ceil(N/16);

% Construct matrices, either sparse or dense
if (sparse_flag)
    A = RandSpMat(N, nnz_per_row);
    B = RandSpMat(N, nnz_per_row);
else
    A = randn(N);
    B = randn(N);
end

% Set the model data
model_data = struct('ib', 1000, ... % Inverse bandwidth
                    'send_oh', 10, ... % Send overhead
                    'recv_oh', 10, ... % Receive overhead
                    'use_parallel_comm', 1);
n_msgs_list = 2;
out_cycles = zeros(length(n_msgs_list), 1);
i = 1;
n_msgs_per_proc = 2;
[C, t] = MM_AsyncCannon(A, B, n_procs, rows_per_proc, n_msgs_per_proc, model_data);
M = C - A*B;
t
norm(M(:))

% Construct random matrix
% N: number of rows and columns
% nnz_per_row: number of non-zero values per row, i.e., the sparsity
function M = RandSpMat(N, nnz_per_row)
    rows = [];
    cols = [];
    vals = [];
    for i = 1:N
        rows = [rows; i*ones(nnz_per_row, 1)];
        cols = [cols; randperm(N, nnz_per_row)'];
        vals = [vals; randn(nnz_per_row, 1)];
    end
    
    M = sparse(rows, cols, vals, N, N);
end