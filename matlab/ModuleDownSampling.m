classdef ModuleDownSampling < handle
    % y = downsample(x,n) decreases the sample rate of x by keeping the
    % first sample and then every nth sample after the first. If x is a
    % matrix, the function treats each column as a separate sequence.

    properties (SetAccess = public)
        N = 1;
        buffer = [];
    end

    properties (SetAccess = private)

    end

    methods (Access = public)
        function obj = ModuleDownSampling(n)
            if (rem(n, 1) ~= 0) && (n <= 0)
                error('N must be an integer and > 0')
            else
                obj.N = n;
            end
            obj.resetBuffer();
        end

        function out_data = process(obj, signal)
            % Decreases the sample rate of signal
            obj.buffer = [obj.buffer; signal];
            out_data = [];
            LEN_BUFFER = length(obj.buffer);
            idx = 0; 
            while (idx + obj.N) <= LEN_BUFFER
                idx = idx + obj.N;
                out_data = [out_data; obj.buffer(idx, :)];
            end
            if (idx < LEN_BUFFER)
                obj.buffer = obj.buffer(idx+1:end, :);
            else
                obj.buffer = [];
            end
        end
        function obj = resetBuffer(obj)
            obj.buffer = [];
        end

    end


end
