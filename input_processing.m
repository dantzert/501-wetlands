%% depth time series processing

% output_filename is the input to the state space wetland model

% input_filename is the results of the SWMM simulation


%process("FWS_A_u.csv","FWS_A.csv",0,400) % if we go off the end it will just stop reading (can set to Inf)
%process("FWS_B_u.csv","FWS_B.csv",0,400) 
%process("FWS_C_u.csv","FWS_C.csv",0,400) 
process("nitrate_removal_u.csv","nitrate_removal.csv",0,400) 
process("quantity_control_u.csv","quantity_control.csv",0,400)
%process("uncontrolled_u.csv","uncontrolled.csv", 0, 400)

function process(output_filename, input_filename, start_day, duration_days)
    %% make $u$, which is the wetland_depth timeseries and frequency analysis of that signal
    %
    % in the future, could incorporate water quality parameters as well as
    % these would be relevant to plant health / which plants do well
    opts = detectImportOptions(input_filename);
    opts.SelectedVariableNames = ["Var1","wetland_depth_m"]; % Var1 is the time vector
    opts.DataLines = [2 + 3600*start_day*24,2 + 3600*start_day*24 + 3600*24*duration_days];
    % column 1 is time in seconds (dt for pyswmm script is 1 second)
    % column 2 is wetland depth in meters
    wetland_depth = readmatrix(input_filename,opts);
    
    %% frequency domain analysis
        
    % oscillation frequencies
    hourly = 1/3600; % 1/ seconds in an hour
    daily = 1/(3600*24); % 1/seconds in a day
    biweekly = 1/(24*3600*14); % 1/seconds in two weeks
    
    for i=1:(3600*24):length(wetland_depth) % only update frequency estimates every day
        dt = 1;
        Fs = dt;
        L = length(wetland_depth(1:i,2));
        f = Fs*(0:floor(L /2) )/ L;
        
        Y = fft(wetland_depth(1:i,2)); % only what we've seen so far
        P2 = abs(Y/L);
        P1 = P2(1:floor(L/2)+1);
        P1(2:end-1) = 2*P1(2:end-1);
        
        %plot(f,P1)
        %title("single sided amplitude spectrum of wetland depth")
        %xlabel("f (Hz)")
        %ylabel("|P1(f)|")
        
        slow_oscillation_strength = sum((f <= biweekly)' .* P1);
        wetland_depth(i,3) = slow_oscillation_strength;
        medium_oscillation_strength = sum( (f > biweekly)' .*(f <= daily)' .* P1);
        wetland_depth(i,4) = medium_oscillation_strength;
        fast_oscillation_strength = sum( (f>daily)' .* (f <= hourly)' .* P1);
        wetland_depth(i,5) = fast_oscillation_strength;
        rapid_oscillation_strength = sum( (f > hourly)' .* P1);
        wetland_depth(i,6) = rapid_oscillation_strength;
    
    end
    
    % interpolate frequency columns (using previous nonzero entry)
    for speed = 3:1:6
        prev_nonzero_entry = 0; 
        for i=1:1:length(wetland_depth)
            if (wetland_depth(i,speed) > 0)
                prev_nonzero_entry = wetland_depth(i,speed);
            else
                wetland_depth(i,speed) = prev_nonzero_entry;
            end
        end
    end


    
    % wetland_depth now has 4 more columns
    % column 3: slow oscillation strength
    % up to 6 which is the fastest oscillation strength
    
    % now finished building the input (u)
    % it's the immediate water level plus four bins of the strength of
    % different frequencies of oscillation
    
    writematrix(wetland_depth,output_filename);
    return;

end
