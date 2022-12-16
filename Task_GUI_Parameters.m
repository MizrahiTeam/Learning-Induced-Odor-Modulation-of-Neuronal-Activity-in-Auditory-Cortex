classdef Task_GUI_Parameters < handle
    
    % ============ Fields =================================
    
    properties
        
        is_save_data    % whether to save data at end of task
        data_dir        % directory for data saving
        data_file       % data file name
        mouse           % mouse name     

        ITI                 % inter trial interval , in sec
        ITI_range           % maximum random deviation for ITI , in sec
        pGo                 % probability of go stimulus
        response_window     % time for mouse to get reinforcement after stimulus , in sec
        task_dur            % duration of entire task , in min
        n_licks             % number of licks to reinforcement

        Go          % frequency of go stimulus , in KHz
        NoGo        % frequency of NoGo stimulus , in KHz
        tone_dur    % duration of tone stimuli , in ms
        atten       % attenuation of tone from max_voltage , in dB

        is_auto_reward      % whether to automatically deliver reward after Go
        is_delay_punishment % whether to punish a false alarm by canceling next stimulus
        is_noise_punishment % whether to punish a false alarm with white noise
        is_prelick_delay    % whether to delay stim on early licks 
        is_repeat_incorrect % whether to repeat NoGo when mouse does FA
        
        is_catch            % whether to present catch trials
        n_freqs             % number of total frequencies (minimum of 2)
        freqs               % a list of presented frequencies , in KHz
        freq_probs          % an array for appearance probabilities for all frequencies
        pNoGo               % probability of NoGo stimulus in case of intermediate frequencies
        catch_odor1_prob    % probability of odor1 before catch stimulus
        catch_odor2_prob    % probability of odor2 before catch stimulus
        catch_odor3_prob    % probability of odor3 before catch stimulus

        is_odor             % whether to present an odor before stimuli
        Go_odor1_prob       % probability of odor1 before Go stimulus
        NoGo_odor1_prob     % probability of odor1 before NoGo stimulus
        Go_odor2_prob       % probability of odor2 before Go stimulus
        NoGo_odor2_prob     % probability of odor2 before NoGo stimulus
        Go_odor3_prob       % probability of odor3 before Go stimulus
        NoGo_odor3_prob     % probability of odor3 before NoGo stimulus
        odor_dur            % duration of odor , in ms
        odor_onset          % onset of odor before stimulus , in ms
        is_odor_alone       % are there odor alone trials
        p_odor_alone        % probability of odor alone trials
        
        reward_dur           % duration of reward, in ms
        reward_delay         % delay of reward from end of stim, in seconds
        noise_punishment_dur % duration of white noise as punishment , in sec
        delay_punishment_dur % duration of delay as punishment , in sec
        
    end
    
    properties (Constant = true)
        
        % Channels
        daq_dev = 'Dev1' ;          % data acquisition device
        stim_ch = 'ao0' ;           % analog output channel for auditory stimuli
        stimOn_ch = 'Port0/Line2' ; % digital output channel for stimulus onset
        lick_ch = 'Port0/Line4' ;   % digital input channel for detecting licks
        reward_ch = 'Port0/Line3' ; % digital output channel for delivering rewards
        odor1_ch = 'Port0/Line7' ;  % digital output channel for delivering odor1
        odor2_ch = 'Port0/Line6' ;  % digital output channel for delivering odor2
        odor3_ch = 'Port0/Line5' ;  % digital output channel for delivering odor3
        dummy_input_ch = 'ai0' ;    % dummy analog input channel to please MATLAB
        MFC_ch = 'ao1' ;            % analog output channel to control MFC
        
        % gui parameters 
        disp_max_time_span = 60 ;   % maximum time span for display , seconds
        update_rate = 20 ;          % the update rate of display, lick detection, and outputs
        trials_per_anal = 10 ;      % trials per online anlss
        video_sample_rate = 3 ;     % frame rate for recording video
        video_run_rate = 5 ;        % frame rate for running video
        scatter_size = 100 ;        % size of markers in scatter plots
        acq_sample_rate = 500 ;     % data acquisition sample rate , KHz
        out_sample_rate = 500 ;     % output sample rate for real mouse , KHz
        prot_delay = 3 ;            % delay from pressing start to delivering stimuli , seconds

        % outputs parameters
        min_freq = 4 ;          % minimum frequency for pure tones , KHz
        max_freq = 32 ;         % maximum frequency for pure tones , KHz
        max_dist_in_oct = 300 ; % maximum distance between Go and NoGo in percetage-octaves
        min_reward_dur = 5  ;   % duration of reward pulse , in ms
        max_reward_dur = 200 ;  % duration of reward pulse , in ms
        base_tone_voltage = 9 ; % voltage for unattenuated tone , volts
        max_atten = 90 ;        % maximum attenuation for pure tones , in dB 
        ramp_ratio = 0.05 ;     % ratio of ramp duration to entire tone duration
        n_pulses = 2 ;           % the number of pulses detected to be considered a lick
        max_n_freqs = 20 ;      % the maximal number of frequencies in case of 
                                % catch trials
        MFC_volt = 0.5 ;        % voltage to MFC which controls odor flow
        
        cal_freqs = 4*2.^(0:0.25:2.5) ; % frequecies for attenuation calibration
        cal_attens = [3.1239 7.1118 9.2751 1.0912 0 0.5606 1.2290 14 4.5716 11.4580 1] ; % decrease in dB for cal_freqs
        
        % task parameters
        
        prelick_delay = 3 ;     % delay when mouse licks prematurely
        
    end
    
    properties (Access = private)
        
        trial_count = 0 ;       % Counter of concluded trials
        
    end
    
    % ============ Events =================================
    
    events
        
        progress_anlss % event for online analysis progression
        
    end
    
    % ============ Methods ================================
    
    methods
        
        function obj = Task_GUI_Parameters()
            % Constructor Task_GUI_Parameters initializes a new Task_GUI_Parameters
            % object that is used for storing internal GUI parameters that
            % are not controlled by user and to count concluded task trials
            %
                       
        end
        
        function add_user_parameters(obj,handles)
            % function add_user_parameters adds all user controlled
            % parameters to the object
            %
            % input:
            % handles - task gui handles object
            %
            
            obj.is_save_data = handles.data_checkbox.Value ;
            obj.data_dir = handles.data_dir_button.UserData ;
            obj.data_file = handles.data_name.String ;
            obj.mouse = handles.mouse.String ;
            
            obj.ITI = str2double(handles.ITI.String) ;
            obj.ITI_range = str2double(handles.ITI_range.String) ;
            obj.pGo = str2double(handles.pGo.String) ;
            obj.response_window = str2double(handles.response_window.String) ;
            obj.task_dur = str2double(handles.task_dur.String) ;
            obj.n_licks = str2double(handles.n_licks.String) ;
            
            obj.Go = str2double(handles.Go.String) ;
            obj.NoGo = str2double(handles.NoGo.String) ;
            obj.tone_dur = str2double(handles.tone_dur.String) ;
            obj.atten = str2double(handles.atten.String) ;
            
            obj.is_auto_reward = handles.auto_reward_checkbox.Value ;
            obj.is_delay_punishment = handles.delay_punishment_checkbox.Value ;
            obj.is_noise_punishment = handles.noise_punishment_checkbox.Value ;
            obj.is_prelick_delay = handles.prelick_checkbox.Value ;
            obj.is_repeat_incorrect = handles.repeat_inc_checkbox.Value ;
            
            obj.is_catch = handles.catch_checkbox.Value ;
            obj.n_freqs = str2double(handles.n_freqs.String) ;
            obj.pNoGo = str2double(handles.pNoGo.String) ;
            obj.freqs = logspace(log10(obj.Go),log10(obj.NoGo),obj.n_freqs) ;
            
            obj.freq_probs = (1 - obj.pGo - obj.pNoGo) * ones(1,obj.n_freqs) / (obj.n_freqs - 2) ;
            obj.freq_probs(1) = obj.pGo ;
            obj.freq_probs(end) = obj.pNoGo ;
            
            obj.is_odor = handles.odor_checkbox.Value ;
            obj.Go_odor1_prob = str2double(handles.Go_odor1_prob.String) ;
            obj.NoGo_odor1_prob = str2double(handles.NoGo_odor1_prob.String) ;
            obj.catch_odor1_prob = str2double(handles.catch_odor1_prob.String) ;
            obj.Go_odor2_prob = str2double(handles.Go_odor2_prob.String) ;
            obj.NoGo_odor2_prob = str2double(handles.NoGo_odor2_prob.String) ;
            obj.catch_odor2_prob = str2double(handles.catch_odor2_prob.String) ;
            obj.Go_odor3_prob = str2double(handles.Go_odor3_prob.String) ;
            obj.NoGo_odor3_prob = str2double(handles.NoGo_odor3_prob.String) ;
            obj.catch_odor3_prob = str2double(handles.catch_odor3_prob.String) ;
            obj.odor_dur = str2double(handles.odor_dur.String) ;
            obj.odor_onset = str2double(handles.odor_onset.String) ;
            obj.is_odor_alone = handles.odor_alone_checkbox.Value ;
            obj.p_odor_alone = str2double(handles.p_odor_alone.String) ;
            
            obj.reward_dur = str2double(handles.reward_dur.String) ;
            obj.reward_delay = str2double(handles.reward_delay.String) ;
            obj.noise_punishment_dur = str2double(handles.noise_punishment_dur.String) ;
            obj.delay_punishment_dur = str2double(handles.delay_punishment_dur.String) ;
                
        end
        
        function val = get_trial_count(obj)
            % function get_trial_count gives the current trial_count as
            % registered by this object
            % 
            % output:
            % val - the trial count
            %
            
            val = obj.trial_count ;
            
        end
        
        function increase_trial_count(obj)
            % Function increase_trial_count increases the concluded trials
            % counter by one and if enough trials were concluded it
            % notifies to progress the online anlss
            %
  
            obj.trial_count = obj.trial_count + 1 ;
            
            if mod(obj.trial_count,obj.trials_per_anal) == 0
                notify(obj , 'progress_anlss')
            end
        
        end
        
        function zero_trial_count(obj)
            % Function zero_trial_count zeros the concluded trials counter 
            %
            
            obj.trial_count = 0 ;
            
        end  
        
    end
    
end