%% import file
tic
[filename, filepath] = uigetfile ...
    ('*.*', 'Please select a trajectory file');

if filename == 0
    return
end

%% pre processing
T_all = dlmread([filepath, filename]);
dim_tra = T_all(1, 1); % the dimension
num_tra = T_all(2, 1); % the number of trajectories
num_vec_tra = T_all(3:end, 2);
T_mat = T_all(3:end, 3:end);
T = cell(1, size(T_mat, 1));

for i = 1:size(T_mat, 1)
    T{i} = reshape(T_mat(i, 1 : num_vec_tra(i)*dim_tra), dim_tra, []);
    T{i} = T{i}';
end

TT = []
for i = 1:200
    One = ones(size(T{i}, 1), 1) * i % add the index of trajectory
    T{i} = horzcat(T{i}, One) % append index
    TT = [TT; T{i}]; % append T{i}
end
dlmwrite('hurricane1950_2006_trans300.XYI', TT, 'delimiter', ' ');
