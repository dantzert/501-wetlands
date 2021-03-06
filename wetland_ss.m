%% wetland state space model

%wetland_response("FWS_A_wetland_response.jpg","FWS_A_u.csv")
close all
%wetland_response("FWS_B_wetland_response.jpg","FWS_B_u.csv")
close all
%wetland_response("FWS_C_wetland_response.jpg","FWS_C_u.csv")
close all
%wetland_response("nitrate_removal_wetland_response.jpg","nitrate_removal_u.csv")
close all
%wetland_response("uncontrolled_wetland_response.jpg","uncontrolled_u.csv")
close all
%wetland_response("quantity_control_wetland_response.jpg","quantity_control_u.csv")
close all

function wetland_response(output,u_file)
    wetland_depth = readmatrix(u_file);
    
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
    speed = 10^-6;
    strong = 1*speed;
    weak = 0.3*speed;
    
    A = [0 -weak 0 -weak -weak -weak;
        -weak 0 0 0 0 0;
        strong 0 0 -weak -weak -weak; 
        strong 0 strong 0 0 -weak ;
        weak 0 strong -weak 0 -weak;
        weak/4 0 weak/4 weak/4 weak/4 0] - eye(6)*strong;
    
    % include self-inhibitition to get stability
    eig(A);
    
    %% input matrix ($B$)
    B = [0 3*strong weak -weak -weak; 0 0 0 weak weak; ...
        strong weak 0 -weak -weak; weak 0 0 0 0; ...
        strong 0 0 0 0; 0 0 0 0 0];
    
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
    wetland = ss(A,B,C,0,'StateName',state_names,'InputName',["level","slow","med","fast","rapid"]);
    %% simulate wetland response
    x_o = [1;1;1;1;1;1]; % initial condition
    [y,t,x] = lsim(wetland,u,t,x_o);
    
    figure(1)
    subplot(2,1,1)
    plot(t,wetland_depth(:,2:6))
    legend(["current level (m)","slow","med","fast","rapid"],"Location","northwest")
    title("input characteristics")
    grid on
    
    
    subplot(2,1,2)
    plot(t,x)
    legend(state_names,"Location","northwest")
    title("wetland response")
    grid on
    set(gcf,'Position',get(0,'Screensize'));
    exportgraphics(gcf,output,"Resolution",300)
end
