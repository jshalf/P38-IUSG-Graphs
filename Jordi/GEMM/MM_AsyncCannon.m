function [C, t] = MM_AsyncCannon(A, B, n_procs, n_per_proc, n_msgs_per_proc, model_data)
    global tile_gemm_count;
    global tot_tile_gemms_count;
    
    [A_nrows, A_ncols] = size(A);
    [B_nrows, B_ncols] = size(B);
    if (A_nrows / n_procs ~= round(A_nrows / n_procs))
        fprintf('ERROR: number of processors must divide number of rows\n');
        C = nan;
        return;
    end
    if (n_per_proc / n_msgs_per_proc ~= round(n_per_proc / n_msgs_per_proc))
        fprintf('ERROR: number of messages per proc must divide number of rows per proc\n');
        C = nan;
        return;
    end
    
    ib = model_data.ib;
    send_oh = model_data.send_oh;
    recv_oh = model_data.recv_oh;
    use_parallel_comm = model_data.use_parallel_comm;

    p_nrows = A_nrows / n_procs;
    msg_nrows = p_nrows / n_msgs_per_proc;
    
    n_tile_gemms = zeros(n_procs);
    tile_gemm_count = zeros(n_procs);
    tot_tile_gemms = 0;
    for i = 1:n_procs
        for j = 1:n_procs
            for n = 1:n_msgs_per_proc
                for m = 1:n_msgs_per_proc
                    for k = 1:n_msgs_per_proc
                        n_tile_gemms(i, j) =  n_tile_gemms(i, j) + 1;
                        tot_tile_gemms = tot_tile_gemms + 1;
                    end
                end
            end
        end
    end
    tot_tile_gemms = n_procs * tot_tile_gemms;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Partition A and B
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A_cell = cell(n_procs, n_procs, n_msgs_per_proc, n_msgs_per_proc);
    B_cell = cell(n_procs, n_procs, n_msgs_per_proc, n_msgs_per_proc);
    for i = 1:n_procs
        for j = 1:n_procs
            for n = 1:n_msgs_per_proc
                for m = 1:n_msgs_per_proc            
                    ii = (i-1)*p_nrows + (n-1)*msg_nrows + (1:msg_nrows);
                    jj = (j-1)*p_nrows + (m-1)*msg_nrows + (1:msg_nrows);
                    A_cell{i, j, n, m} = A(ii, jj);
                    B_cell{i, j, n, m} = B(ii, jj);
                end
            end
        end
    end
    A = A_cell;
    B = B_cell;
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup data that is local (non-shared) to each processor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    c = cell(n_procs, n_procs, n_msgs_per_proc, n_msgs_per_proc);
    a = cell(n_procs, n_procs, n_msgs_per_proc, n_msgs_per_proc);
    b = cell(n_procs, n_procs, n_msgs_per_proc, n_msgs_per_proc);
    iter = ones(n_procs);
    for i = 1:n_procs
        for j = 1:n_procs
            for n = 1:n_msgs_per_proc
                for m = 1:n_msgs_per_proc
                    c{i, j, n, m} = [];
                end
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup comp/comm flags
    %%%%%%%%%%%%%%%%%%%%%%%%%
    A_recv = cell(n_procs, n_procs, n_procs, n_msgs_per_proc, n_msgs_per_proc);
    B_recv = cell(n_procs, n_procs, n_procs, n_msgs_per_proc, n_msgs_per_proc);
    A_send = cell(n_procs, n_procs, n_procs, n_msgs_per_proc, n_msgs_per_proc);
    B_send = cell(n_procs, n_procs, n_procs, n_msgs_per_proc, n_msgs_per_proc);
    
    A_recv_flag = zeros(n_procs, n_procs, n_procs, n_msgs_per_proc, n_msgs_per_proc);
    B_recv_flag = zeros(n_procs, n_procs, n_procs, n_msgs_per_proc, n_msgs_per_proc);
    A_send_flag = zeros(n_procs, n_procs, n_procs, n_msgs_per_proc, n_msgs_per_proc);
    B_send_flag = zeros(n_procs, n_procs, n_procs, n_msgs_per_proc, n_msgs_per_proc);
    
    A_send_buff_flag = zeros(n_procs, n_procs, n_procs, n_msgs_per_proc, n_msgs_per_proc);
    B_send_buff_flag = zeros(n_procs, n_procs, n_procs, n_msgs_per_proc, n_msgs_per_proc);
    
    A_SubGEMM_flag = zeros(n_procs, n_procs, n_procs, n_msgs_per_proc, n_msgs_per_proc);
    B_SubGEMM_flag = zeros(n_procs, n_procs, n_procs, n_msgs_per_proc, n_msgs_per_proc);
    A_SubGEMM_not_done = ones(n_procs, n_procs, n_procs, n_msgs_per_proc, n_msgs_per_proc);
    B_SubGEMM_not_done = ones(n_procs, n_procs, n_procs, n_msgs_per_proc, n_msgs_per_proc);
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Setup simulation data
    %%%%%%%%%%%%%%%%%%%%%%%%
    proc_busy = zeros(n_procs);
    proc_tot_time = zeros(n_procs);
    if (use_parallel_comm)
        A_comm_actor_busy = zeros(n_procs, n_procs, n_procs, n_msgs_per_proc, n_msgs_per_proc);
        B_comm_actor_busy = zeros(n_procs, n_procs, n_procs, n_msgs_per_proc, n_msgs_per_proc);
    else
        A_comm_actor_busy = zeros(n_procs, n_procs, n_procs);
        B_comm_actor_busy = zeros(n_procs, n_procs, n_procs);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Post first sends
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:n_procs
        for j = 1:n_procs
            for n = 1:n_msgs_per_proc
                for m = 1:n_msgs_per_proc
                    A_send{1, i, j, n, m} = A{i, j, n, m};
                    B_send{1, i, j, n, m} = B{i, j, n, m};
                    A_send_buff_flag(1, i, j, n, m) = 1;
                    B_send_buff_flag(1, i, j, n, m) = 1;
                end
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Start async Cannon's
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %for l = 1:n_procs
    t = 0;
    tot_tile_gemms_count = 0;
    comm_actor_wtime = 0;
    proc_wtime = 0;
    % Loop until all tiles have been multiplied
    while (tot_tile_gemms_count < tot_tile_gemms)
        %fprintf('%d %d\n', tot_tile_gemms_count, tot_tile_gemms);
        proc_wtime_start = tic;
        % Loop over all processors and move to the next iteration
        % if the number of MMs has been completed for this iteration
        for i = 1:n_procs
            for j = 1:n_procs
                l = iter(i, j);
                if (l <= n_procs)
                    %if (tile_gemm_count(i, j) == n_tile_gemms(i, j))
                    if (sum(A_SubGEMM_not_done(l, i, j, :)) == 0 && ...
                        sum(B_SubGEMM_not_done(l, i, j, :)) == 0)
                        for n = 1:n_msgs_per_proc
                            for m = 1:n_msgs_per_proc
                                a{i, j, n, m} = [];
                                b{i, j, n, m} = [];
                            end
                        end
                        tile_gemm_count(i, j) = 0;
                        iter(i, j) = iter(i, j) + 1;
                        %fprintf('%d %d %d\n', i, j, iter(i, j));
                    end
                end
            end
        end
        proc_wtime = proc_wtime + toc(proc_wtime_start);
        
        % Communication loop
        comm_actor_wtime_start = tic;
        for l = 1:n_procs
        for i = 0:(n_procs-1)
        for j = 0:(n_procs-1)
            if (l == 1)
                i_A = i+1;
                j_A = mod(n_procs+j-i, n_procs)+1;                            
                i_B = mod(n_procs+i-j, n_procs)+1;
                j_B = j+1;
            else
                i_A = i+1;
                j_A = mod(j+n_procs-1, n_procs)+1;
                i_B = mod(i+n_procs-1, n_procs)+1;
                j_B = j+1;
            end
            if (use_parallel_comm)
                for n = 1:n_msgs_per_proc
                for m = 1:n_msgs_per_proc
                    if (~A_comm_actor_busy(l, i+1, j+1, n, m))
                        if (~isempty(A_send{l, i+1, j+1, n, m}) && isempty(A_recv{l, i_A, j_A, n, m}))
                            if (A_send_flag(l, i+1, j+1, n, m))
                                A_recv{l, i_A, j_A, n, m} = A_send{l, i+1, j+1, n, m};
                                A_send{l, i+1, j+1, n, m} = [];
                            else
                                A_comm_actor_busy(l, i+1, j+1, n, m) = ...
                                    A_comm_actor_busy(l, i+1, j+1, n, m) + ...
                                        ib * nnz(A_send{l, i+1, j+1, n, m});
                                A_send_flag(l, i+1, j+1, n, m) = 1;
                            end
                        end
                    end

                    if (~B_comm_actor_busy(l, i+1, j+1, n, m))
                        if (~isempty(B_send{l, i+1, j+1, n, m}) && isempty(B_recv{l, i_B, j_B, n, m}))
                            if (B_send_flag(l, i+1, j+1, n, m))
                                B_recv{l, i_B, j_B, n, m} = B_send{l, i+1, j+1, n, m};
                                B_send{l, i+1, j+1, n, m} = [];
                            else
                                B_comm_actor_busy(l, i+1, j+1, n, m) = ...
                                    B_comm_actor_busy(l, i+1, j+1, n, m) + ...
                                        ib * nnz(B_send{l, i+1, j+1, n, m});
                                B_send_flag(l, i+1, j+1, n, m) = 1;                                 
                            end
                        end
                    end
                end
                end
            else
                if (~A_comm_actor_busy(l, i+1, j+1))
                    for n = 1:n_msgs_per_proc
                        for m = 1:n_msgs_per_proc
                            if (A_send_buff_flag(l, i+1, j+1, n, m))
                                if (A_send_flag(l, i+1, j+1, n, m))
                                    %fprintf('At time %d: comm actor (%d, %d, %d) deliver A\n', t, l, i, j);
                                    A_recv{l, i_A, j_A, n, m} = A_send{l, i+1, j+1, n, m};
                                    A_send{l, i+1, j+1, n, m} = [];
                                    A_send_buff_flag(l, i+1, j+1, n, m) = 0;
                                else
                                    %fprintf('At time %d: comm actor (%d, %d, %d) send A\n', t, l, i, j);
                                    A_comm_actor_busy(l, i+1, j+1) = ...
                                        A_comm_actor_busy(l, i+1, j+1) + ...
                                            ib * nnz(A_send{l, i+1, j+1, n, m});
                                    A_send_flag(l, i+1, j+1, n, m) = 1;
                                end
                                break;
                            end
                        end
                        if (A_comm_actor_busy(l, i+1, j+1))
                            break;
                        end
                    end
                end
                
                if (~B_comm_actor_busy(l, i+1, j+1))
                    for n = 1:n_msgs_per_proc
                        for m = 1:n_msgs_per_proc
                            if (B_send_buff_flag(l, i+1, j+1, n, m))
                                if (B_send_flag(l, i+1, j+1, n, m))
                                    %fprintf('At time %d: comm actor (%d, %d, %d) deliver B\n', t, l, i, j);
                                    B_recv{l, i_B, j_B, n, m} = B_send{l, i+1, j+1, n, m};
                                    B_send{l, i+1, j+1, n, m} = [];
                                    B_send_buff_flag(l, i+1, j+1, n, m) = 0;
                                else
                                    %fprintf('At time %d: comm actor (%d, %d, %d) send B\n', t, l, i, j);
                                    B_comm_actor_busy(l, i+1, j+1) = ...
                                        B_comm_actor_busy(l, i+1, j+1) + ...
                                            ib * nnz(B_send{l, i+1, j+1, n, m});
                                    B_send_flag(l, i+1, j+1, n, m) = 1;                                
                                end
                                break;
                            end
                        end
                        if (B_comm_actor_busy(l, i+1, j+1))
                            break;
                        end
                    end
                end
            end
        end
        end
        end
        comm_actor_wtime = comm_actor_wtime + toc(comm_actor_wtime_start);
        
        proc_wtime_start = tic;
        % Computation loop
        for i = 1:n_procs
        for j = 1:n_procs
            l = iter(i, j);
            if (l <= n_procs && proc_busy(i, j) == 0)
                for n = 1:n_msgs_per_proc
                    for m = 1:n_msgs_per_proc
                        if (~isempty(A_recv{l, i, j, n, m}))
                            if (A_recv_flag(l, i, j, n, m))
                                if (~A_SubGEMM_flag(l, i, j, n, m))
                                    %fprintf('At time %d: proc (%d, %d) A SubGEMM\n', t, i, j);
                                    A_SubGEMM_flag(l, i, j, n, m) = 1;
                                    a{i, j, n, m} = A_recv{l, i, j, n, m};
                                    [c, n_ops] = ASubGEMM(c, a{i, j, n, m}, b, i, j, n, m);
                                    proc_busy(i, j) = proc_busy(i, j) + n_ops;
                                    if (l < n_procs)
                                        proc_busy(i, j) = proc_busy(i, j) + send_oh;
                                    end
                                else
                                    if (l < n_procs)
                                        A_send{l+1, i, j, n, m} = a{i, j, n, m};
                                        A_send_buff_flag(l+1, i, j, n, m) = 1;
                                    end
                                    A_recv{l, i, j, n, m} = [];
                                    A_SubGEMM_not_done(l, i, j, n, m) = 0;
                                end
                            else
                                %fprintf('%d %d %d %d %d\n', l, i, j, n, m);
                                A_recv_flag(l, i, j, n, m) = 1;
                                proc_busy(i, j) = proc_busy(i, j) + recv_oh;
                            end
                            break;
                        end
                        %fprintf('hello %d\n', isempty(B_recv{l, i, j, n, m}));
                        if (~isempty(B_recv{l, i, j, n, m}))
                            if (B_recv_flag(l, i, j, n, m))
                                if (~B_SubGEMM_flag(l, i, j, n, m))
                                    %fprintf('At time %d: proc (%d, %d) B SubGEMM\n', t, i, j);
                                    B_SubGEMM_flag(l, i, j, n, m) = 1;
                                    b{i, j, n, m} = B_recv{l, i, j, n, m};
                                    [c, n_ops] = BSubGEMM(c, a, b{i, j, n, m}, i, j, n, m);
                                    proc_busy(i, j) = proc_busy(i, j) + n_ops;
                                    if (l < n_procs)
                                        proc_busy(i, j) = proc_busy(i, j) + send_oh;
                                    end
                                else
                                    %fprintf('At time %d: proc (%d, %d) send B\n', t, i, j);
                                    if (l < n_procs)
                                        B_send{l+1, i, j, n, m} = b{i, j, n, m};
                                        B_send_buff_flag(l+1, i, j, n, m) = 1;
                                    end
                                    B_recv{l, i, j, n, m} = [];
                                    B_SubGEMM_not_done(l, i, j, n, m) = 0;
                                end
                            else
                                %fprintf('At time %d: proc (%d, %d) recv B\n', t, i, j);
                                B_recv_flag(l, i, j, n, m) = 1;
                                proc_busy(i, j) = proc_busy(i, j) + recv_oh;
                            end
                            break;
                        end
                    end
                    if (proc_busy(i,j))
                        break;
                    end
                end 
                proc_tot_time(i, j) = proc_tot_time(i, j) + proc_busy(i, j);
            end
        end
        end
        proc_wtime = proc_wtime + toc(proc_wtime_start);
        
        proc_busy = max(0, proc_busy-1);
        A_comm_actor_busy = max(0, A_comm_actor_busy-1);
        B_comm_actor_busy = max(0, B_comm_actor_busy-1);
        
        
        t = t + 1;
        %fprintf('%d %d %d\n', t, tot_tile_gemms_count, tot_tile_gemms);
    end
    %fprintf('%d %f %f\n', t, proc_wtime, comm_actor_wtime);

    C = sparse(A_nrows, B_ncols);
    for i = 1:n_procs
        for j = 1:n_procs
            for n = 1:n_msgs_per_proc
                for m = 1:n_msgs_per_proc
                    ii = (i-1)*p_nrows + (n-1)*msg_nrows + (1:msg_nrows);
                    jj = (j-1)*p_nrows + (m-1)*msg_nrows + (1:msg_nrows);
                    C(ii, jj) = c{i, j, n, m};
                end
            end
        end
    end
end

function [c, n_ops] = ASubGEMM(c, a_ijnm, b, i, j, n, m)
    global tile_gemm_count;
    global tot_tile_gemms_count;

    n_ops = 0;
    [~, ~, n_msgs_per_proc, ~] = size(c);
    for k = 1:n_msgs_per_proc
        if (~isempty(b{i, j, m, k}))
            if (isempty(c{i, j, n, k}))
                c{i, j, n, k} = a_ijnm * b{i, j, m, k};
            else
                c{i, j, n, k} = c{i, j, n, k} + a_ijnm * b{i, j, m, k};
            end
            tile_gemm_count(i, j) = tile_gemm_count(i, j) + 1;
            tot_tile_gemms_count = tot_tile_gemms_count + 1;
            
            [n_rows_a, n_cols_a] = size(a_ijnm);
            [n_rows_b, n_cols_b] = size(b{i, j, m, k});
            n_ops = n_ops + nnz(b{i, j, m, k})*n_rows_a + nnz(a_ijnm)*n_cols_b;%n_rows_a * n_cols_a * n_cols_b;
        end
    end
end

function [c, n_ops] = BSubGEMM(c, a, b_ijnm, i, j, n, m)
    global tile_gemm_count;
    global tot_tile_gemms_count;

    n_ops = 0;
    [~, ~, n_msgs_per_proc, ~] = size(c);
    for k = 1:n_msgs_per_proc
        if (~isempty(a{i, j, k, n}))
            if (isempty(c{i, j, k, m}))
                c{i, j, k, m} = a{i, j, k, n} * b_ijnm;
            else
                c{i, j, k, m} = c{i, j, k, m} + a{i, j, k, n} * b_ijnm;
            end
            tile_gemm_count(i, j) = tile_gemm_count(i, j) + 1;
            tot_tile_gemms_count = tot_tile_gemms_count + 1;
            
            [n_rows_a, n_cols_a] = size(a{i, j, k, n});
            [n_rows_b, n_cols_b] = size(b_ijnm);
            n_ops = n_ops + nnz(b_ijnm)*n_rows_a + nnz(a{i, j, k, n})*n_cols_b;%n_rows_a * n_cols_a * n_cols_b;
        end
    end
end