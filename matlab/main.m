

clear;  clc;
% addpath(genpath('data'));
% addpath('utility');
close all
%% 
% Raw folder on Cardinal
% folder_path = "data/Cardinal_skin_states";
% folder_name = "20230201_1613";
%%
folder_path ="/Users/raymond/Documents/repo/raymond_dev/body-detection/matlab/data_demo";
folder_name = "20230201_1613_ca_skin";

% Read configuration of data folders
data_config = readtable('/Users/raymond/Documents/repo/raymond_dev/body-detection/matlab/data_demo/data_configuration.txt');

% notify the resolution for 1g of the data
[is_folder_noted,folder_id] = max(contains(data_config.date_time,folder_name));
if is_folder_noted == 1
    disp("The resolution for 1g of the data: " + data_config.resolution_1g(folder_id))
else
    disp("There is no information of 1g-resolution")
end

% initial setting
resolution_1g = 2^8;     
initial_skin_status = 0; % off-skin

% group data in seconds, then detect
skin = OnBodyDetection(initial_skin_status,resolution_1g);

info_data_error = 1;
error = skin.group_into_second_data(folder_path,folder_name,info_data_error);
IR_SAMPLING_RATE = 32; 
IR_SAMPLING_DURATION = 2;
skin.set_ir_sampling_parameter(IR_SAMPLING_RATE,IR_SAMPLING_DURATION); 

skin.determine_skin_state(skin.second_data);
skin.plot_skin_states(data_config);

