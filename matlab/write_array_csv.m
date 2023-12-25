function write_array_csv(fileID, array_data)

if (isempty(array_data))
    fprintf(fileID, '\n');
else
    if length(array_data) == 1
        fprintf(fileID, '%d\n', array_data(1));
    end
    if length(array_data) > 1
        fprintf(fileID, '%d', array_data(1));
        fprintf(fileID, ', %d', array_data(2:end));
        fprintf(fileID, '\n');
    end
end
end