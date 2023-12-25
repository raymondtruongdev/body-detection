% second_data = skin.second_data;
% skin_states = skin.skin_states;

path_save = '/Users/raymond/Documents/repo/raymond_dev/body-detection/matlab/data_demo/20230201_1613_ca_skin_forC';

fileID_accel100 = fopen(fullfile(path_save, 'accel100hz.csv'), 'w');

fileID_ir = fopen(fullfile(path_save, 'ppg_ir.csv'), 'w');
fileID_accel100_filter = fopen(fullfile(path_save, 'accel100hz_filter.csv'), 'w');
fileID_accel12h5 = fopen(fullfile(path_save, 'accel12h5.csv'), 'w');
fileID_skin = fopen(fullfile(path_save, 'skin.csv'), 'w');

for k = 1:length(second_data)

    data = second_data(k);

    write_array_csv(fileID_accel100, k);
    write_array_csv(fileID_accel100, round(data.axlX));
    write_array_csv(fileID_accel100, round(data.axlY));
    write_array_csv(fileID_accel100, round(data.axlZ));

    write_array_csv(fileID_accel100_filter, k);
    write_array_csv(fileID_accel100_filter, round(data.accel100hFiltered(:,1)));
    write_array_csv(fileID_accel100_filter, round(data.accel100hFiltered(:,2)));
    write_array_csv(fileID_accel100_filter, round(data.accel100hFiltered(:,3)));

    write_array_csv(fileID_accel12h5, k);
    write_array_csv(fileID_accel12h5, round(data.accel12h5(:,1)));
    write_array_csv(fileID_accel12h5, round(data.accel12h5(:,2)));
    write_array_csv(fileID_accel12h5, round(data.accel12h5(:,3)));

    write_array_csv(fileID_ir, k);
    write_array_csv(fileID_ir, round(data.ir'));

    write_array_csv(fileID_skin, k);
    write_array_csv(fileID_skin, round(skin_states));


end
fclose(fileID_accel100);
fclose(fileID_ir);
fclose(fileID_skin);
fclose(fileID_accel12h5);

%% Filtering by feed each second data block
% second_data = skin.second_data;

load('second_data_and_skin.mat');

a = 1;
b = [22504, 24378, 26033, 27436, 28559, 29378, 29877, 30044, 29877, 29378, 28559, 27436, 26033, 24378, 22504];

myModuleFir = ModuleFIR(a, b);
myModuleFir.resetBuffer();
myModuleDownSampling = ModuleDownSampling(100/12.5);

for k = 1:length(second_data)
    accel_raw = [second_data(k).axlX, second_data(k).axlY,second_data(k).axlZ];
    second_data(k).accel100hFiltered = myModuleFir.process(accel_raw);
    second_data(k).accel12h5 = myModuleDownSampling.process(second_data(k).accel100hFiltered);

end
accel100hFiltered = vertcat(second_data(:).accel100hFiltered);

%% Filtering with full data array
data_in = vertcat(second_data(:).axlX);
accel100hBuildIn = filter_data(data_in);

%%
figure(1);
clf("reset");
hold on;
plot(data_in, 'r');
plot(accel100hBuildIn, '*g');
plot(accel100hFiltered(:, 1), 'blue');
plot(accel100hFiltered(:, 1)-accel100hBuildIn, 'black');

%%
accel12h5X = vertcat(second_data(:).axlX_12_5);
accel12h5 = vertcat(second_data(:).accel12h5);
accel12h5XN = accel12h5(:, 1);
figure(1);
clf("reset");
hold on;
plot(accel12h5X(:), 'r');
plot(accel12h5XN, 'green');
plot(accel12h5X-accel12h5XN, 'black');
