function fData = filter_data(data)
% collected data are in 100 Hz
% downgrade them to 12.5 Hz

b = [22504, 24378, 26033, 27436, 28559, 29378, 29877, 30044, 29877, 29378, 28559, 27436, 26033, 24378, 22504];
a = 1;

fData = filter(b, a, data);
fData = bitshift(fData, -19, 'int64');

end