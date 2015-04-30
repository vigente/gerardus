% takes the trajectories data ('trajs.txt') from the Smoldyn output and
% rearranges it to match the format for my model
sim_index = ['A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';'P';'Q';'R';'S';'T';'U';'V';'W';'X';'Y';'Z'];
sim_index2 = ['AA';'AB';'AC';'AD';'AE';'AF';'AG';'AH';'AI';'AJ';'AK';'AL';'AM';'AN';'AO';'AP';'AQ';'AR';'AS';'AT';'AU';'AV';'AW';'AX';'AY';'AZ';'BA';'BB';'BC';'BD'];

for k = 1:size(sim_index, 1)
    j = sim_index(k);
    for i = 1:4
        filetoimport = sprintf('trajs%s%d.txt', j, i);
        trajs = dlmread(filetoimport);

        % remove the unneeded columns
        trajs = trajs(:,[1 4:6]);

        % calculate number of timesteps and number of molecules
        timesteps = trajs(:,1);
        no_of_timesteps = max(timesteps)-1; % to allow for missing data
        no_of_molecules = nnz(timesteps ==1);

        % rearrange matrix
        positions = zeros(no_of_timesteps, no_of_molecules, 3);
        for j = 1:no_of_timesteps
             timestep_data = trajs((no_of_molecules*(j-1)+1):no_of_molecules*j,:);
             positions(j,:,:) = timestep_data(:,2:4);
        end

        % order to be (molecules, xyz, timestep) as for my model
        positions = permute(positions,[2 3 1]);
        eval(sprintf('positions%d = positions', i));
        filetosave = sprintf('positions_%s%d.mat', j, i);

        save(filetosave, sprintf('positions%d', i), '-v7.3');
        clearvars -except i
    end
end
for k = 1:size(sim_index2, 1)
    j = sim_index2(k,:);
    for i = 1:4
        filetoimport = sprintf('trajs%s%d.txt', j, i);
        trajs = dlmread(filetoimport);

        % remove the unneeded columns
        trajs = trajs(:,[1 4:6]);

        % calculate number of timesteps and number of molecules
        timesteps = trajs(:,1);
        no_of_timesteps = max(timesteps)-1; % to allow for missing data
        no_of_molecules = nnz(timesteps ==1);

        % rearrange matrix
        positions = zeros(no_of_timesteps, no_of_molecules, 3);
        for j = 1:no_of_timesteps
             timestep_data = trajs((no_of_molecules*(j-1)+1):no_of_molecules*j,:);
             positions(j,:,:) = timestep_data(:,2:4);
        end

        % order to be (molecules, xyz, timestep) as for my model
        positions = permute(positions,[2 3 1]);
        eval(sprintf('positions%d = positions', i));
        filetosave = sprintf('positions_%s%d.mat', j, i);

        save(filetosave, sprintf('positions%d', i), '-v7.3');
        clearvars -except i
    end
end
