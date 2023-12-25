classdef OnBodyDetection < handle
    % % This class is used to implement on-off-wrist algorithm
    % % To use the class:
    %    - skin = OnBodyDetection(1,9): create a class of SkinDetection that
    %                                 has the initial skin state is on skin
    %                                 (1), and the resolution of 9 bits for
    %                                 1g acceleration
    %                                 
    %    - skin.determine_skin_state(second_data): determine skin states every seconds
    %
    %    - skin.plot_skin_states: to display the plot of state transition. There are
    %                              3 states in switching on-wrist and off-wrist:
    %                              On-skin, Off-skin, and Transition
    
    properties (SetAccess = private)
        initial_status % 0 = OFF_SKIN ;    1 = ON_SKIN;        2 = TRANSITION 
        second_data       
        resolution_1g
        error_output 
        data_folder
        IR_SAMPLING_DURATION % seconds
        IR_SAMPLING_RATE % 

        % Threshold constants
        ON_RELATIVE_FLIPPING_INDEX = 0.42;
        OFF_RELATIVE_FLIPPING_INDEX = 0.6;
        SPIKE_RELATIVE_FLIPPING_INDEX = 0.05;
        MINIMAL_PPG_IR = 2.2*10^6;
        MAX_IR_SD = 4; % max IR sampling duration in seconds, should be > IR_SD
        PPG_TIMEOUT_DURING_ON = 60; % seconds
        THRES_ACCEL_RANGE_SECOND_ACTIVE = 200;
        PPG_IR_STD = 770;
    end
    
    properties (SetAccess = protected)
        skin_states % for the whole session
        heart_rates
        timestamps
        debug_data
        raw_accel
        raw_ppg
        raw_amb
        raw_hr
        raw_var
        raw_ppg_no_dc

        OFF_SKIN = 0;
        ON_SKIN = 1;
        TRANSITION = 2;
    end
    
    methods
        function obj = OnBodyDetection(initial_skin_status, resolution_1g)
            arguments
                initial_skin_status = 0;  
                resolution_1g = 2^8;
            end

            obj.initial_status = initial_skin_status;
            obj.resolution_1g = resolution_1g;

        end
        
        %%
        function set_ir_sampling_parameter(obj, sampling_rate, sampling_duration)
            
            obj.IR_SAMPLING_DURATION = sampling_duration;
            obj.IR_SAMPLING_RATE = sampling_rate;
        end
        
        %%
        function error = group_into_second_data(obj,folder_path,folder_name, info_data_error)
        % Read raw data of acceleration and PPG. 
        % Then, group them into data batches for every second in the same
        % duration of accel and PPG.
        % The error has ouput 1 when there is any PPG batch missed in a
        % second.

            buffer_data = struct;
            obj.data_folder = folder_name;
            error = 0;
            csv_files= dir(fullfile(folder_path,folder_name,'*.csv'));
            acxl = struct;   
            accel_12_5Hz = struct;
            amb = struct;     
            ppg = struct;      
            hr = struct;
            
            LIST_NAME_AXL       = {'ACCEL_100HZ.CSV'  , 'RAW_ACC.CSV'} ;
            LIST_NAME_PPG       = {'PPG_25HZ.CSV'     , 'RAW_PPG.CSV'}   ;
            LIST_NAME_PPG_NO_DC = {'RAW_PPG_DRIVER.CSV'} ;
            LIST_NAME_AMB       = {'AMB_25HZ.CSV'     , 'AMB.CSV'}   ;
            LIST_NAME_HR        = {'HR_1HZ.CSV'       , 'HR.CSV',  'HRPERHZ_STEP.CSV'} ;
           
            for csvId=1:length(csv_files)
                csv_filename =  upper(csv_files(csvId).name);
                switch csv_filename
                    case LIST_NAME_AXL

                        data_table     =  readtable(fullfile( csv_files(csvId).folder,csv_filename ));
                        acxl.timestamp = data_table.Timestamp  ;
                        acxl.axlX     = data_table.X ;
                        acxl.axlY     = data_table.Y ;
                        acxl.axlZ     = data_table.Z ;

                        obj.raw_accel = acxl;
                        
                        
                    case LIST_NAME_AMB
                        data_table     =  readtable(fullfile( csv_files(csvId).folder,csv_filename ));
                        amb.timestamp  = data_table.Timestamp ;
                        amb.amb       = data_table.AMB ;    
                        
                    case LIST_NAME_PPG
                        data_table     =  readtable(fullfile( csv_files(csvId).folder,csv_filename ));
                        table_field = data_table.Properties.VariableNames;
                        
                        % extract the right field of optical sensor
                        ppg.timestamp  = data_table.Timestamp;
                        LIST_FIELD_PPG = {'PPG', 'Green'};      
                        LIST_FIELD_Red = {'Red'};
                        LIST_FIELD_Ir = {'Ir'};

                        [is_green, ppg_id] = max(contains(table_field,LIST_FIELD_PPG));
                        [is_red, red_id] = max(contains(table_field,LIST_FIELD_Red));
                        [is_ir, ir_id] = max(contains(table_field,LIST_FIELD_Ir));

                        if is_green
                            field_name = data_table.Properties.VariableNames{ppg_id};
                            ppg.green = data_table.(field_name);
                        end
                        if is_red
                            field_name = data_table.Properties.VariableNames{red_id};
                            ppg.red = data_table.(field_name);
                        end
                        if is_ir
                            field_name = data_table.Properties.VariableNames{ir_id};
                            ppg.ir = data_table.(field_name);
                        end

                        obj.raw_ppg = ppg;
                   
                        % extract the right field of AMB
                        LIST_FIELD_AMB = {'Amb'};
                        [is_amb, amb_id] = max(contains(table_field,LIST_FIELD_AMB));
                        
                        if is_amb
                            amb.timestamp  = data_table.Timestamp;
                            field_name = data_table.Properties.VariableNames{amb_id};
                            amb.amb = data_table.(field_name);

                            obj.raw_amb = amb;
                        end

                    case LIST_NAME_PPG_NO_DC
                        data_table     =  readtable(fullfile( csv_files(csvId).folder,csv_filename ));
                        table_field = data_table.Properties.VariableNames;
                        
                        % extract the right field of optical sensor
                        ppg_no_dc.timestamp  = data_table.Timestamp;
                        LIST_FIELD_PPG = {'PPG', 'Green'};      
                        LIST_FIELD_PPG_DC_OFFSET = {'Green_Dc_Offset'};
                        LIST_FIELD_Red = {'Red'};
                        LIST_FIELD_Red_DC_OFFSET = {'Red_Dc_Offset'};
                        LIST_FIELD_Ir = {'Ir'};
                        LIST_FIELD_Ir_DC_OFFSET = {'Ir_Dc_Offset'};


                        [is_green, ppg_id] = max(contains(table_field,LIST_FIELD_PPG));
                        [is_red, red_id] = max(contains(table_field,LIST_FIELD_Red));
                        [is_ir, ir_id] = max(contains(table_field,LIST_FIELD_Ir));
                        [is_green_offset, ppg_offset_id] = max(contains(table_field,LIST_FIELD_PPG_DC_OFFSET));
                        [is_red_offset, red_offset_id] = max(contains(table_field,LIST_FIELD_Red_DC_OFFSET));
                        [is_ir_offset, ir_offset_id] = max(contains(table_field,LIST_FIELD_Ir_DC_OFFSET));

                        if is_green
                            field_name = data_table.Properties.VariableNames{ppg_id};
                            ppg_no_dc.green = data_table.(field_name);

                            if is_green_offset
                                field_name = data_table.Properties.VariableNames{ppg_offset_id};
                                ppg_no_dc.green_dc_offset = data_table.(field_name);
                            end
                        end
                        if is_red
                            field_name = data_table.Properties.VariableNames{red_id};
                            ppg_no_dc.red = data_table.(field_name);

                            if is_red_offset
                                field_name = data_table.Properties.VariableNames{red_offset_id};
                                ppg_no_dc.red_dc_offset = data_table.(field_name);
                            end
                        end
                        if is_ir
                            field_name = data_table.Properties.VariableNames{ir_id};
                            ppg_no_dc.ir = data_table.(field_name);

                            if is_ir_offset
                                field_name = data_table.Properties.VariableNames{ir_offset_id};
                                ppg_no_dc.ir_dc_offset = data_table.(field_name);
                            end
                        end

                        obj.raw_ppg_no_dc = ppg_no_dc;

                    case LIST_NAME_HR
                        data_table     =  readtable(fullfile( csv_files(csvId).folder,csv_filename ));
                        hr.timestamp   = data_table.Timestamp ;
                        table_field = data_table.Properties.VariableNames;
                        
                        % extract the right field of HR
                        LIST_FIELD_HR = {'HR', 'heartrate'};
                        
                        [is_hr, hr_id] = max(contains(table_field,LIST_FIELD_HR));
                        
                        if is_hr
                            field_name = data_table.Properties.VariableNames{hr_id};
                            hr.hr = data_table.(field_name);
                        end

                        % extract the right field of HR quality
                        LIST_FIELD_HRQ = {'quality', 'heartrate_quality'};
                        [is_hrq, hrq_id] = max(contains(table_field,LIST_FIELD_HRQ));
                        
                        if is_hrq
                            field_name = data_table.Properties.VariableNames{hrq_id};
                            hr.quality = data_table.(field_name);
                        end

                        obj.raw_hr = hr;
                end
            end

            if  isempty(fieldnames(acxl)) && isempty(fieldnames(ppg))
                return;
            end
            
            % convert to data on Cardinal
            % +/- 8g data using 12 bits and 12.5Hz accel data -> 1g = 2^8
            raw_8g_data = floor([acxl.axlX, acxl.axlY, acxl.axlZ]*(2^8)/(obj.resolution_1g));
            fData = obj.filter_data(raw_8g_data);
            fData = round(fData);
            accel_12_5Hz.timestamp = acxl.timestamp(8:8:end);
            accel_12_5Hz.axlX = fData(8:8:end,1);
            accel_12_5Hz.axlY = fData(8:8:end,2);
            accel_12_5Hz.axlZ = fData(8:8:end,3);

            % check timestamp consistence among data sets
            [~,ts_accel] = groupcounts(round(acxl.timestamp/1000));
            [~,ts_ppg] = groupcounts(round(ppg.timestamp/1000));

            % extract common timestamp
            ts = (max(ts_accel(1),ts_ppg(1)):min(ts_accel(end),ts_ppg(end)));
            
            for k=1:length(ts)
                buffer_data(k).timestamp              = ts(k) ;
               
                % get heartrate
                if ~isempty(obj.raw_hr)
                    hr_rows_the_same_second = (round(hr.timestamp/1000) == ts(k));
                    hr_value = hr.hr(hr_rows_the_same_second);
                    hrq_value = hr.quality(hr_rows_the_same_second);
                    
                    if isempty(hr_value)
                        buffer_data(k).hr                     = 0 ;
                        buffer_data(k).hr_quality             = 0 ;
                        buffer_data(k).isNewHR                = 0 ;
                    else
                        buffer_data(k).hr                     = round(mean(hr_value));
                        buffer_data(k).hr_quality             = round(mean(hrq_value));
    
                        if k == 1
                            buffer_data(k).isNewHR            = 1;
                        elseif buffer_data(k).hr == buffer_data(k-1).hr 
                            buffer_data(k).isNewHR            = 0;
                        else 
                            buffer_data(k).isNewHR            = 1;       
                        end
                    end
                else
                    buffer_data(k).hr                     = 0 ;
                    buffer_data(k).hr_quality             = 0 ;
                    buffer_data(k).isNewHR                = 0 ;
                end

                % get ppg data
                ppg_rows_the_same_second = (round(ppg.timestamp/1000) == ts(k));
                buffer_data(k).green = ppg.green(ppg_rows_the_same_second);
                buffer_data(k).red = ppg.red(ppg_rows_the_same_second);
                buffer_data(k).ir = ppg.ir(ppg_rows_the_same_second);

                % buffer_data(k).ac_green = decode_unsigned_decimal(ppg_no_dc.green(ppg_rows_the_same_second),'total_number_of_bit',24,'number_of_sign_bit',3,'within_adc_range',false); 
                % buffer_data(k).ac_red = decode_unsigned_decimal(ppg_no_dc.red(ppg_rows_the_same_second),'total_number_of_bit',24,'number_of_sign_bit',3,'within_adc_range',false);
                % buffer_data(k).ac_ir = decode_unsigned_decimal(ppg_no_dc.ir(ppg_rows_the_same_second),'total_number_of_bit',24,'number_of_sign_bit',3,'within_adc_range',false);

                % notify empty data
                if isempty(ppg.ir(ppg_rows_the_same_second))
                    if info_data_error == 1
                        disp("Empty PPG at (seconds): " + k);
                    end
                    error = 1;
                end
                
                % get accelerometer data 100Hz
                accel_rows_the_same_second = (round(acxl.timestamp/1000) == ts(k));
                buffer_data(k).axlX = acxl.axlX(accel_rows_the_same_second);
                buffer_data(k).axlY = acxl.axlY(accel_rows_the_same_second);
                buffer_data(k).axlZ = acxl.axlZ(accel_rows_the_same_second);

                 % get accelerometer data 12.5Hz
                accel_rows_the_same_second = (round(accel_12_5Hz.timestamp/1000) == ts(k));
                buffer_data(k).axlX_12_5 = accel_12_5Hz.axlX(accel_rows_the_same_second);
                buffer_data(k).axlY_12_5 = accel_12_5Hz.axlY(accel_rows_the_same_second);
                buffer_data(k).axlZ_12_5 = accel_12_5Hz.axlZ(accel_rows_the_same_second);                                   
            end
          
            obj.second_data = buffer_data;
        end
        %%
        function determine_skin_state(obj,second_data)
        % Loop through data batches in every seconds to detect any skin state transition. 
        % Motion is based on 12.5Hz accel data.
        % Metrics of relative flipping index and standard deviation of a 2-second PPG-IR segment are used to switch between skin states   

            PPG_SR = obj.IR_SAMPLING_RATE;
            
            % take parameters of IR-sampling
            if isempty(obj.IR_SAMPLING_DURATION)
                IR_SD = 1;
            else
                IR_SD = obj.IR_SAMPLING_DURATION;
            end
            
            debug = struct;
            % for logging
            ir_segment = zeros(PPG_SR*IR_SD,1); % buffer IR_SD seconds of PPG_IR 
            
            % for algo on device
            ppg_segment = zeros(PPG_SR*IR_SD,1); % buffer IR_SD seconds of PPG_IR
            is_buffering_ppg = 0;
            
            skin_status = obj.initial_status; % take the initial skin state
            ppg_second_counter = 0;
            
            total_sec = length(second_data);
            obj.skin_states = zeros(total_sec,1);
            obj.timestamps = zeros(total_sec,1);
            
            for id_sec = 1:total_sec
                
                sensor = obj.get_specific_data(second_data,id_sec);
                obj.timestamps(id_sec) = sensor.timestamp;

                % buffering PPG if needed in the previous loop
                if is_buffering_ppg == 1
                    ppg_segment = obj.buffer_ppg_segment(ppg_segment, sensor.ppg_ir);
                    ppg_second_counter = ppg_second_counter + 1;
                end
                
                % for logging
                % Use circular buffer to store required length of data
                ir_segment(1:PPG_SR*(IR_SD-1)) = ir_segment(PPG_SR+1:PPG_SR*IR_SD);
                if ~isempty(sensor.ppg_ir)
                    ir_segment(1+PPG_SR*(IR_SD-1):PPG_SR*IR_SD) = sensor.ppg_ir(end-PPG_SR+1:end);
                else
                    ir_segment(1+PPG_SR*(IR_SD-1):PPG_SR*IR_SD) = 0;
                end
                
                motion_range = obj.calculate_feature(sensor);
                
                switch skin_status
                    case obj.OFF_SKIN
                        if (motion_range >= obj.THRES_ACCEL_RANGE_SECOND_ACTIVE)                
                            skin_status        = obj.TRANSITION;   % change to state 'TRANSITION'
                            is_buffering_ppg = 1;
                        end
                        
                    case obj.TRANSITION              
                        % utilize the PPG in the next 4 seconds (at most)
                        if ppg_second_counter >= IR_SD
                            
                            rfi = obj.get_flipping_index(ppg_segment)/(PPG_SR*IR_SD);
                            if (rfi > 0) && ...
                               (rfi <= obj.ON_RELATIVE_FLIPPING_INDEX) && ...
                               (min(ppg_segment) >= obj.MINIMAL_PPG_IR)
                                
                                skin_status   = obj.ON_SKIN;
                                
                                ppg_checking_timer = 1;
                                is_buffering_ppg = 0;
                                ppg_second_counter = 0;

                            elseif ppg_second_counter == obj.MAX_IR_SD
                                skin_status    =  obj.OFF_SKIN;
                                is_buffering_ppg = 0;
                                ppg_second_counter = 0;
                            end
                        end
 
                    case obj.ON_SKIN
                        % using PPG timer to separate times of checking
                        % off-skin
                        
                        if ppg_checking_timer == 1

                            % utilize the PPG in the next 4 seconds (at most)
                            if (motion_range <= obj.THRES_ACCEL_RANGE_SECOND_ACTIVE) || (is_buffering_ppg == 1)
                                is_buffering_ppg = 1;
                            
                                if ppg_second_counter >= IR_SD
                                    
                                    rfi = obj.get_flipping_index(ppg_segment)/(PPG_SR*IR_SD);
                                    if ((rfi <= obj.SPIKE_RELATIVE_FLIPPING_INDEX) && (iqr(ppg_segment) == 0)) || ... % (if any) few positive spikes of PPG-IR 
                                       ((rfi > obj.OFF_RELATIVE_FLIPPING_INDEX) && (std(ppg_segment) <= obj.PPG_IR_STD)) || ... % noise with small deviation
                                       ((rfi > obj.OFF_RELATIVE_FLIPPING_INDEX) && (max(ppg_segment) <= obj.MINIMAL_PPG_IR)) % noise without proper lower PPG-IR threshold      
                                            
                                        skin_status   = obj.OFF_SKIN;
                                        is_buffering_ppg = 0;
                                        ppg_second_counter = 0;

                                    elseif ppg_second_counter == obj.MAX_IR_SD
                                        % reset all if the last check is not satisfied
                                        ppg_checking_timer = 0;
                                        ppg_delay_counter = 0;
                                        is_buffering_ppg = 0;
                                        ppg_second_counter = 0;
                                       
                                    end
                                end                                     
                            end
                        else
                            ppg_delay_counter = ppg_delay_counter + 1;
                            if ppg_delay_counter >= obj.PPG_TIMEOUT_DURING_ON
                                ppg_checking_timer = 1;
                            end
                        end
                end
                
                % log data for plotting
                obj.skin_states(id_sec) = skin_status;
                
                % for debugging
                if isempty(motion_range)
                    debug.axl_range(id_sec) = 0;
                else
                    debug.axl_range(id_sec) = motion_range;
                end
                debug.ir_flipping_index(id_sec) = obj.get_flipping_index(ir_segment)/(PPG_SR*IR_SD);
                debug.std_ir_segment(id_sec) = std(ir_segment);
                debug.ppg_ir(id_sec) = {sensor.ppg_ir'};
                debug.skin_state(id_sec) = skin_status;

            end
            
            obj.debug_data = debug;
        end
        
        %% Plot skin states
        function plot_skin_states(obj,data_config)
            
            % process title
            [is_folder_noted,folder_id] = max(contains(data_config.date_time,obj.data_folder));
            if is_folder_noted == 1
                folder_title = data_config.surface(folder_id) + " on from " + data_config.time_on_skin(folder_id) + " to " + data_config.time_off_skin(folder_id);
            else
                folder_title = obj.data_folder;
            end

            skin_state = obj.skin_states;
            ppg_ir = obj.raw_ppg.ir;

            figure
            ax1 = subplot(3,1,1);
            plot(skin_state,'LineWidth',2);
            set(gca,'XLim',[1 length(skin_state)]);
            set(gca,'YTICK',(0:2), 'yticklabels',["OFF_SKIN"; "ON_SKIN";"TRANSITION"],'TickLabelInterpreter', 'none');
            title(folder_title,'Interpreter','none')

            if is_folder_noted == 1
                hold on
                    xline(data_config.time_on_skin(folder_id),'Color','g','LineWidth',2);
                    xline(data_config.time_off_skin(folder_id),'Color','r','LineWidth',2)
                hold off
            end
            
            ax2 = subplot(3,1,2);
            yyaxis left
            plot(obj.debug_data.ir_flipping_index)
            set(gca,'XLim',[1 length(obj.debug_data.ir_flipping_index)]);
            ylabel('RFI')

            yyaxis right
            plot(obj.debug_data.axl_range)
            set(gca,'XLim',[1 length(obj.debug_data.axl_range)]);
            ylim([0 1000]);
            ylabel('Accel range')
            xlabel('Time (seconds)');

            linkaxes([ax1, ax2],'x');

            subplot(3,1,3);
            plot(ppg_ir)
            set(gca,'XLim',[1 length(ppg_ir)]);
            title('PPG-IR');
            xlabel('Sample')
            
        end

        %% Plot raw PPG data
        function plot_restored_dc_ppg(obj,data_config)
        % Plot the restored-dc values of PPG 
            
            % process title
            [is_folder_noted,folder_id] = max(contains(data_config.date_time,obj.data_folder));
            if is_folder_noted == 1
                folder_title = data_config.surface(folder_id) + ", restored dc PPG" + ", on from " + data_config.time_on_skin(folder_id)*32 + " to " + data_config.time_off_skin(folder_id)*32;
            else
                folder_title = obj.data_folder + ", restored dc PPG";
            end

            ppg = obj.raw_ppg;

            figure
            ax1 = subplot(3,1,1);
            plot(ppg.green,'LineWidth',2,'Color','g');
            ylabel('Green')
            title(folder_title,'Interpreter','none')
            
            ax2 = subplot(3,1,2);
            plot(ppg.red,'LineWidth',2,'Color','r','LineStyle',':');
            ylabel('Red')

            ax3 = subplot(3,1,3);
            plot(ppg.ir,'LineWidth',2,'Color','r');
            ylabel('Infrared')

            linkaxes([ax1, ax2, ax3],'x');                    
        end

        function plot_ac_ppg(obj,data_config)
        % Plot the ac values of PPG to see whether they are in the range

            % process title
            [is_folder_noted,folder_id] = max(contains(data_config.date_time,obj.data_folder));
            if is_folder_noted == 1
                folder_title = data_config.surface(folder_id) + ", ac PPG"+ ", on from " + data_config.time_on_skin(folder_id)*32 + " to " + data_config.time_off_skin(folder_id)*32;
            else
                folder_title = obj.data_folder + ", ac PPG";
            end

            ppg_no_dc = obj.raw_ppg_no_dc;

            decoded_green = decode_unsigned_decimal(ppg_no_dc.green,'total_number_of_bit',24,'number_of_sign_bit',3,'within_adc_range',false);
            decoded_red = decode_unsigned_decimal(ppg_no_dc.red,'total_number_of_bit',24,'number_of_sign_bit',3,'within_adc_range',false);
            decoded_ir = decode_unsigned_decimal(ppg_no_dc.ir,'total_number_of_bit',24,'number_of_sign_bit',3,'within_adc_range',false);

            figure
            ax1 = subplot(3,1,1);
            plot(decoded_green,'LineWidth',2,'Color','g');
            yline([-2^13*171 2^13*171],'LineWidth',2);
            ylabel('Green')
            title(folder_title,'Interpreter','none')
            
            ax2 = subplot(3,1,2);
            plot(decoded_red,'LineWidth',2,'Color','r','LineStyle',':');
            yline([-2^13*171 2^13*171],'LineWidth',2);
            ylabel('Red')

            ax3 = subplot(3,1,3);
            plot(decoded_ir,'LineWidth',2,'Color','r');
            yline([-2^13*171 2^13*171],'LineWidth',2);
            ylabel('Infrared')

            linkaxes([ax1, ax2, ax3],'x');                    
        end

        function plot_dc_offset_ppg(obj,data_config)
        % Plot the dc offset values of PPG    
            % process title
            [is_folder_noted,folder_id] = max(contains(data_config.date_time,obj.data_folder));
            if is_folder_noted == 1
                folder_title = data_config.surface(folder_id) + ", dc offset PPG" + ", on from " + data_config.time_on_skin(folder_id)*32 + " to " + data_config.time_off_skin(folder_id)*32;
            else
                folder_title = obj.data_folder + ", dc offset PPG";
            end

            ppg = obj.raw_ppg_no_dc;

            figure
            ax1 = subplot(3,1,1);
            plot(ppg.green_dc_offset,'LineWidth',2,'Color','g');
            ylabel('Green')
            title(folder_title,'Interpreter','none')
            
            ax2 = subplot(3,1,2);
            plot(ppg.red_dc_offset,'LineWidth',2,'Color','r','LineStyle',':');
            ylabel('Red')

            ax3 = subplot(3,1,3);
            plot(ppg.ir_dc_offset,'LineWidth',2,'Color','r');
            ylabel('Infrared')

            linkaxes([ax1, ax2, ax3],'x');                    
        end
    end
    
    methods (Access = protected)
        %%
        function sensor = get_specific_data(~, second_data, id_sec)
        % Use the 12.5 Hz accel and PPG-IR in the second      
            sensor = struct;
            sensor.timestamp    = second_data(id_sec).timestamp;
            sensor.id           = id_sec;
            sensor.ppg_ir       = second_data(id_sec).ir;
            
            % Use 12.5Hz accel data
            sensor.axlX         = second_data(id_sec).axlX_12_5;
            sensor.axlY         = second_data(id_sec).axlY_12_5;
            sensor.axlZ         = second_data(id_sec).axlZ_12_5;
        end
        
        %%
        function axl_range = calculate_feature(~, sensor)
        % return the sum of maximum acceleration difference in all 3 axes

            maxX = max(sensor.axlX);
            maxY = max(sensor.axlY);
            maxZ = max(sensor.axlZ);
            
            minX = min(sensor.axlX);
            minY = min(sensor.axlY);
            minZ = min(sensor.axlZ);
            
            axl_range = (maxX-minX)+(maxY-minY)+(maxZ-minZ);  
        end

        %%
        function out_ppg_segment = buffer_ppg_segment(obj, ppg_segment, ppg_ir)
        % add in the lastest batch of PPG and output the lastest buffered PPG segment    
            IR_SD = obj.IR_SAMPLING_DURATION;
            PPG_SR = obj.IR_SAMPLING_RATE;

            out_ppg_segment = ppg_segment;
            out_ppg_segment(1:PPG_SR*(IR_SD-1)) = ppg_segment(PPG_SR+1:PPG_SR*IR_SD);

            if ~isempty(ppg_ir)
                out_ppg_segment(1+PPG_SR*(IR_SD-1):PPG_SR*IR_SD) = ppg_ir(end-PPG_SR+1:end);
            else
                % NOTE: this is to handle a situation where ppg_ir is empty in collected
                % data. This is probably not something that can happen when running in
                % device
                out_ppg_segment(1+PPG_SR*(IR_SD-1):PPG_SR*IR_SD) = 0;
            end
            
        end
         %%
        function fData = filter_data(~, data)
        % collected data are in 100 Hz
        % downgrade them to 12.5 Hz

            b = [22504	24378	26033	27436	28559	29378	29877	30044	29877	29378	28559	27436	26033	24378	22504];
            a = 1;
    
            fData = filter(b, a, data);
            fData = bitshift(fData,-19,'int64');

        end
        function flip_index = get_flipping_index(~,data)
            % Compute the quatity of switching side (up / down) of column vectors
            % The input data should have at least 3 values otherwise flip_index = 0

            if size(data,1) >= 3
                data_diff = diff(data);
                flip_index = zeros(1,size(data_diff,2));
                sample_nunber = size(data_diff,1);

                % go to each column
                for col = 1:size(data_diff,2)
                    data_diff_vec = data_diff(:,col);

                    % jump to the first different from zero
                    ids = 1;

                    while  data_diff_vec(ids) == 0
                        ids = ids + 1;
                        if ids == sample_nunber
                            break
                        end
                    end

                    if ids < sample_nunber
                        current_val = data_diff_vec(ids);
                        for idx = ids+1:sample_nunber
                            if current_val*data_diff_vec(idx) < 0
                                flip_index(col) = flip_index(col) + 1;
                                current_val = data_diff_vec(idx);
                            elseif current_val*data_diff_vec(idx) > 0
                                current_val = data_diff_vec(idx);
                            end
                        end
                    end
                end
            else
                flip_index = 0;
            end
        end
    end
end

