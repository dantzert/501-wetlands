%% wetland state space model



%% make $u$, which is the wetland_depth timeseries and frequency analysis of that signal
%
% in the future, could incorporate water quality parameters as well as
% these would be relevant to plant health / which plants do well
opts = detectImportOptions("FWS_C.csv");
opts.SelectedVariableNames = ["Var1","wetland_depth_m"];
opts.DataLines = [3600*90,3600*90 + 30000];
% column 1 is time in seconds (dt for pyswmm script is 1 second)
% column 2 is wetland depth in meters
wetland_depth = readmatrix("FWS_C.csv",opts);

%% frequency domain analysis
    
% oscillation frequencies
hourly = 1/3600; % 1/ seconds in an hour
daily = 1/(3600*24); % 1/seconds in a day
biweekly = 1/(24*3600*14); % 1/seconds in two weeks

for i=1:1:length(wetland_depth)
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


% wetland_depth now has 4 more columns
% column 3: slow oscillation strength
% up to 6 which is the fastest oscillation strength

% now finished building the input (u)
% it's the immediate water level plus four bins of the strength of
% different frequencies of oscillation

%% state variables $x$
%
% I'll consider vegetation, invertebrates, waterfowl (dabblers and divers), 
% and mammals in a simplified linear interaction of a wetland.
%
% *Vegetation*
%
% I'll consider two types, desirable and undesirable. Desirable plants
% create good habitat and food for invertebrates, waterfowl, and mammals.
% They respond well to slow variations in water level and poorly to quick
% drawdowns. 
%
% Undesirable vegetation doesn't provide habitat or food for animals and
% thrives in rapidly changing water levels.
%
% The vegetation is eaten down by all the animals present.
%
% *Invertebrates*
%
% Invertebrates are lumped into one state. They also do better under slow
% drawdowns and benefit from the presence of desirable vegetation. They are
% predated on by the other animals.
%
% *Waterfowl*
%
% The effect of water level on birds is actually nonlinear as they have
% preferred ranges of depths. The way I'll represent this is by divers
% having a greater positive impact from increased water depth (than dabblers) 
% with a negative feedback from the presence of dabblers. This means that
% divers would be outcompeted by dabblers at low water levels, which isn't
% exactly the dynamic, but should give the desired dynamics within the
% constraints of a linear model. 
%
% The waterfowl benefit from the presence of desirable vegetation and
% invertebrates. The frequency of oscillation of water level doesn't affect
% them in this model. They are predated on by mammals.
%
% *Mammals*
%
% Unaffected by water levels. They eat the vegetation and other animals in
% the wetland.
%
% x_1 = desirable veg
% 
% x_2 = undesirable veg
% 
% x_3 = invertebrates
% 
% x_4 = dabblers
%
% x_5 = divers
%
% x_6 = mammals
%

state_names = ["desirable veg","undesirable veg", "invertebrates", "dabblers", "divers", "mammals"];

%% interaction matrix ($A$)
% qualitative magnitude of effect
strong = 1*10^-4;
weak = 0.5*10^-4;

A = [0 -strong -weak -weak -weak -weak;
    -weak 0 0 0 0 0;
    strong weak 0 -weak -weak -weak; 
    weak 0 strong 0 0 -weak ;
    weak 0 strong -weak 0 -weak;
    weak/2 0 weak/2 weak/2 weak/2 0] - eye(6)*0.001

% include self-inhibitition to get stability


%% input matrix ($B$)
B = [0 10*strong 0 -weak -strong; 0 0 0 10*weak 10*strong; ...
    weak weak 0 -weak -strong; weak 0 0 0 0; ...
    strong 0 0 0 0; 0 0 0 0 0]

%% observation matrix ($C$)
%
% Suppose that the only method of monitoring the wetland relies on
% automated processing of audio from a microphone in the wetland. Then we
% can only observe the birds and mammals as the aquatic invertebrates
% aren't picked up and vegetation is silent. Suppose also that we can
% distinguish each type of audible animal perfectly.
%
C = [zeros(3,3) eye(3,3)];


%% define state space model 
u = wetland_depth(:,2:6);
t = wetland_depth(:,1);
wetland = ss(A,B,C,0,'StateName',state_names,'InputName',["level","slow","med","fast","rapid"])
%% simulate wetland response
x_o = [1;1;1;1;1;1]; % initial condition
[y,t,x] = lsim(wetland,u,t,x_o);

figure()
title("wetland response")
plot(t,x)
legend(state_names,"Location","East")
grid on