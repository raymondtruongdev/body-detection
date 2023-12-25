classdef ModuleFIR < handle
    %  This class is used to filter the data in vector signal_in with the
    %  filter described by vectors A and B to create the filtered
    %  data signal_out.  The filter is a "Direct Form II Transposed"
    %  implementation of the standard difference equation:
    %    a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
    %                         - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
    %
    %   If a(1) is not equal to 1, FILTER normalizes the filter
    %   coefficients by a(1).
    %
    % Using coefficients with this class:
    % a = 1;
    % b = [22504	24378	26033	27436	28559	29378	29877	30044	29877	29378	28559	27436	26033	24378	22504];
    % data_in = [-563, -580, -616, -637, -631, -602, -573, -544, -515, -483, -456, -427, -395, -361, -326, ...
    %            -288, -248, -215, -184, -161, -138, -113, -93, -81, -85, -106, -128, -152, -170, -181, ...
    %            -194, -205, -222, -246, -282, -328, -359, -351, -327, -319, -308, -281, -261, -271, ...
    %            -300, -317, -301, -242, -182, -176, -189, -179, -144, -102, -68, -39, -12, 17, 58, 100, ...
    %             142, 195, 238, 252, 220, 171, 136, 135, 167, 196, 203, 177, 119, 50, -42, -201, -374, -471, ...
    %            -519, -554, -547, -517, -506, -493, -461, -433, -406, -378, -343, -311, -284, -266, -252, ...
    %            -234, -225, -206, -175, -137, -103, -86, -69];
    % my_iir_x = ModuleFIR(a,b);
    % my_iir_x.resetBuffer();
    % fData = my_iir_x.ApplyFilter(data_in);

    % Theory:
    % We can use Matlab design filter to create the coefficients :
    % H = design(fdesign.lowpass('N,Fp,Fst', 14, .05, .0625), 'ALLFIR');
    % b = H(2).Numerator;
    % a = 1;

    % Using coefficients with Matlab build-in function:
    % a = 1;
    % b = [22504, 24378, 26033, 27436, 28559, 29378, 29877, 30044, 29877, 29378, 28559, 27436, 26033, 24378, 22504];
    % % b = bitshift(b,-19);
    % fData = filter(b, a, data_in);
    % fData = bitshift(fData, -19, 'int64');


    properties
        a = [];
        b = [];
        buffer_in = [];
        buffer_out = [];
    end

    properties (SetAccess = private)
        MAX_BUFFER_LENGTH = 1;
    end

    methods
        function obj = ModuleFIR(a, b)
            obj.a = a;
            obj.b = b;
            obj.MAX_BUFFER_LENGTH = max([length(a), length(b)]);
            obj.resetBuffer();
        end

        function signal_out = process(obj, signal_in)
            coeff_a = obj.a;
            coeff_b = obj.b;
            % using previous data
            if isempty(obj.buffer_in)
                obj.buffer_in = zeros(obj.MAX_BUFFER_LENGTH, size(signal_in, 2));
                obj.buffer_out = obj.buffer_in;
            end
            data_in = [obj.buffer_in; signal_in];
            data_out = [obj.buffer_out; zeros(size(signal_in, 1), size(signal_in, 2))];

            % Get the index of 1st sample of new data in array
            start_index = obj.MAX_BUFFER_LENGTH + 1;
            % Only apply filter for new data
            for n = start_index:length(data_in)
                filterOutput = zeros(1, size(signal_in, 2));

                for i = 1:length(coeff_b)
                    if (n >= i+1)
                        filterOutput = filterOutput + coeff_b(i) * data_in(n - i+1,:);
                    end
                end

                for i = 2:length(coeff_a)
                    if (n >= i+1)
                        filterOutput = filterOutput - (coeff_a(i) * data_out(n - i+1,:));
                    end
                end
                filterOutput = bitshift(filterOutput, -19, 'int64');
                data_out(n, :) = filterOutput;

            end
            % output the filtered data
            signal_out = data_out(start_index:end, :);
            % get buffer for next calculation
            obj.buffer_in = data_in(end-obj.MAX_BUFFER_LENGTH+1:end, :);
            obj.buffer_out = data_out(end-obj.MAX_BUFFER_LENGTH+1:end, :);

        end
        function obj = resetBuffer(obj)
            obj.buffer_in = [];
            obj.buffer_out = [];
        end
    end
end