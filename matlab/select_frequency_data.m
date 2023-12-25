function [newdata] = select_frequency_data(data,freq_mode)

switch freq_mode

    case 'mode_data_00' % raw data
        %%% Current filter
        %         H = design(fdesign.lowpass('N,Fp,Fst',14,.05,.0625),'ALLFIR');
        %         b = H(2).Numerator;
        %         a = 1;

        b = [22504	24378	26033	27436	28559	29378	29877	30044	29877	29378	28559	27436	26033	24378	22504];
        %         b = bitshift(b,-19);
        a = 1;

        fData = filter(b, a, data);
        fData = fData(8:8:end,:); % change to 12.5hz
        fData = bitshift(fData,-19,'int64');

    case 'mode_data_01' % data 12.5 Hz
        fData = data;

    case 'mode_data_02' % from 100hz then convert to 'duplicate' 100Hz
        data1= zeros(size(data));
        for k=1:size(data,1)
            if mod(k,2)==1
                data1(k,:) = data(k,:);
            else
                data1(k,:) = data(k-1,:);
            end
        end

        b= [22504	24378	26033	27436	28559	29378	29877	30044	29877	29378	28559	27436	26033	24378	22504]/(2^19);
        a = 1;
        fData = filter(b, a, data1);
        fData = fData(1:8:end,:); % change to 12.5hz

    case 'mode_data_03' % from 50hz android duplicate to'duplicate' 100Hz
        data1= zeros(size(data,1)*2,size(data,2));
        for k=1:size(data,1)
            data1((k-1)*2+1,:) = data(k,:);
            data1(k*2,:) = data(k,:);

        end

        b= [22504	24378	26033	27436	28559	29378	29877	30044	29877	29378	28559	27436	26033	24378	22504]/(2^19);
        a = 1;
        fData = filter(b, a, data1);
        fData = fData(1:8:end,:); % change to 12.5hz

    case 'mode_data_04' % from 50hz android
        % calculate for interger in C
        data2= data;
        a2 = 1;
        b2 =  [48756	54872	58756	60088	58756	54872	48756];
        fData2 = int32(filter(b2, a2, data2));
        fData2 = bitshift(fData2,-19);
        fData2 = fData2(4:4:end,:); % change to 12.5hz

        % return result
        fData = fData2;

end

newdata = fData;

end