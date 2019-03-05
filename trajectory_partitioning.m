%% import file
tic
[filename, filepath] = uigetfile ...
    ('*.*', 'Please select a trajectory file');

if filename == 0
    return
end

%% pre-processing
T_all = dlmread([filepath, filename]);
dim_tra = T_all(1, 1);
num_tra = T_all(2, 1);
num_vec_tra = T_all(3:end, 2);
T_mat = T_all(3:end, 3:end);
T = cell(1, size(T_mat, 1));

for i = 1:size(T_mat, 1)
    T{i} = reshape(T_mat(i, 1:num_vec_tra(i) * dim_tra), dim_tra, []);
    T{i} = T{i}';
end

%% algorithm
% Approximate Trajectory Partitioning
global point_list;
point_list = T{6};
global partion_list;
partion_list = [];
% add p1
partion_list = horzcat(partion_list, point_list(1, :));
start_idx = 1;
length = 1;
partition_idx = [];
partition_idx = horzcat(partition_idx, 1);
len = size(point_list, 1);

while (start_idx + length < len)
    current_idx = start_idx + length;
    cost_par = MDL_par(start_idx, current_idx);
    cost_nonpar = MDL_nonpar(start_idx, current_idx);

    if (cost_par > cost_nonpar)
        partion_list = vertcat(partion_list, point_list(current_idx - 1, :));
        partition_idx = horzcat(partition_idx, current_idx);
        start_idx = current_idx - 1;
        length = 1;
    else
        length = length + 1;
    end

end

partion_list = vertcat(partion_list, point_list(len, :));
partition_idx = horzcat(partition_idx, len);

% plot
plot(point_list(:, 1), point_list(:, 2), '-b*');
hold on
plot(partion_list(:, 1), partion_list(:, 2), '-bo');
axis square

eva = evaluate(partion_list, partition_idx)
time = toc

%% MDL par
function y = MDL_par(start_idx, current_idx)
    global point_list;
    path = point_list(start_idx, :) - point_list(current_idx, :);
    L = log2_(norm(path));
    L_DH = vertical_func(start_idx, current_idx) + theta_func(start_idx, current_idx);
    y = L + L_DH;
end

function MDL_no = MDL_nonpar(start_idx, current_idx)%%original path
    global point_list;
    sum_ = 0;

    for i = start_idx:current_idx - 1
        path = point_list(i + 1, :) - point_list(i, :);
        sum_ = sum_ + log2_(norm(path));
    end

    MDL_no = sum_;
end

%% vertical function
function y = vertical_func(start_idx, current_idx)
    sum = 0;
    global point_list;
    par_vec = point_list(current_idx, :) - point_list(start_idx, :); % p1-->p4
    par_vec = horzcat(par_vec, 0);
    vert_dist_list = []; % Euclidean distance list: l1, l2, l3.....

        for j = start_idx:current_idx
        path_tem_ver = point_list(start_idx, :) - point_list(j, :); % p1-->p2, p1-->p3
        path_tem_ver = horzcat(path_tem_ver, 0);

        cross_product_ver = cross(path_tem_ver, par_vec);
        verti_dist = norm(cross_product_ver) / norm(par_vec); % l1, l2, l3
        vert_dist_list = horzcat(vert_dist_list, verti_dist); % add to list
    end

    for i = 1:size(vert_dist_list, 2) - 1
        l1 = vert_dist_list(1, i)^2 + vert_dist_list(1, i + 1)^2;
        l2 = vert_dist_list(1, i) + vert_dist_list(1, i + 1);
        % common case: l1 + l2 == 0, divide by 0
        if l1 + l2 == 0
            sum = 0;
        else
            sum = sum + log2_(l1 / l2);
        end

    end

    y = sum;
end

%% theta function
function y = theta_func(start_idx, current_idx)
    global point_list;
    sum = 0;
    partion_path = point_list(current_idx, :) - point_list(start_idx, :);
    partion_path = horzcat(partion_path, 0);

    for j = start_idx:current_idx - 1
        path_tem = point_list(j + 1, :) - point_list(j, :);
        path_tem = horzcat(path_tem, 0);
        cross_product = cross(path_tem, partion_path);
        theta_dist = norm(cross_product) / norm(partion_path);
        sum = sum + log2_(theta_dist);
    end

    y = sum;
end

%% evaluate the cost
function y = evaluate(partition, partition_idx_arr)
    L = 0;
    row = size(partition, 1);
    L_DH = 0;

    for j = 1:row - 1
        sub_path = partition(j + 1, :) - partition(j, :);
        L = L + log2_(norm(sub_path));

        sub_start_idx = partition_idx_arr(1, j);
        sub_end_idx = partition_idx_arr(1, j + 1);
        vert_dist = vertical_func(sub_start_idx, sub_end_idx);
        theta_dist = theta_func(sub_start_idx, sub_end_idx);
        L_DH = L_DH + vert_dist + theta_dist;
    end

    y = L + L_DH
end

%% log2
function y = log2_(input)
    if (input < 1)
        y = 0;
        % y=log2(input);
    else
        y = log2(input);
    end
    % y = ceil(y);
end
