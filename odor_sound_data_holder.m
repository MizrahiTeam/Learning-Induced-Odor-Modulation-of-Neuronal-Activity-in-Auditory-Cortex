classdef odor_sound_data_holder < handle & matlab.mixin.Heterogeneous
    
        % ============ Fields =================================
    
    properties
        
        Parameters
        Go_times
        NoGo_times
        Lick_times
        Reward_times
        Punishment_times
        Go_odors
        NoGo_odors
        
        odor_order % odors ordered according to probabilities
        
        trial_types
        trial_times
        trial_odors
        is_early_punish
        
        hit_rate
        FA_rate
        d_prime
        bias
        d_prime_exp
        odor_d_prime
        mean_dDI
        
        BOTs
        n_neurons
        field
        tSOS
        tSOF
        tSOL
        tSOR
        mask = logical.empty([0 0 0 ]) ;
        Time
        
        ft
        lt
        ft_all
        lt_all
        
        Neurons
        
    end
    
    properties (Constant = true)
        
        frames_ch = 1 ;
        StimOn_ch = 2 ;
        reward_ch = 3 ;
        lick_ch = 4 ;
        n_channels = 4 ;
        
        frames_before_stim = 16 ;   % Number of frames to be taken before stimulus
        frames_after_stim = 37  ;   % Number of frames to be taken after stimulus
        baseline_frames = 7     ;   % Number of frames to be taken as baseline (starting from first frame) 
        frames_for_amp = 7     ;    % Number of frames to be taken into account for ampltude calculation
        smooth_fact = 1 ;           % Smoothing for dFF traces
        
        WS_downsample = 100     ;   % Downsample factor for WS vectors
        
        LICK_RATE_CRIT = 10 ;
        LICK_RATE_CRIT_TRIALS = 50 ;
        
        EARLY_PUNISH_CRIT = 0.8 ;
                
    end
    
    % ============ Methods ================================
    
    methods
        
        function obj = odor_sound_data_holder()
        end 
        
        function organize_data(obj)
            
            si = mode(diff(obj.tSOF)) ;
            obj.Time = ((1-obj.frames_before_stim):obj.frames_after_stim) * si ;
            
            obj.classify_trials()
            obj.choose_first_last_trials()
            obj.create_neurons()
            obj.calc_statistics()
            
        end
        
        function classify_trials(obj)
            
            p = obj.Parameters ; 
                    
            if isempty(p.n_licks)
                p.n_licks = 1 ;
            end   
            
            temp = [obj.Go_times obj.NoGo_times] ;
            temp(end) = [] ;
            obj.trial_types = zeros(size(temp)) ;
            obj.is_early_punish = zeros(size(temp)) ;
            obj.trial_odors = zeros(size(temp)) ;
            [obj.trial_times,order] = sort(temp) ;
            Sound = [ones(size(obj.Go_times)) zeros(size(obj.NoGo_times))] ;
            Odor = [obj.Go_odors obj.NoGo_odors] ;
            Sound = Sound(order) ;
            Odor = Odor(order) ;
            tSOL = obj.tSOL ;
            tSOS = obj.tSOS ;
            tone_dur = p.tone_dur*1e-3 ;
            odor_dur = p.odor_dur*1e-3 ;
            response_window = p.response_window ;
            for trial = 1 : length(temp)            
                obj.trial_odors(trial) = find(obj.odor_order == Odor(trial)) ;
                Odor_start = tSOS(trial) - odor_dur ;
                Start = tSOS(trial) + tone_dur ;
                End = Start + response_window ;
                prelick = logical(sum(tSOL > Odor_start & tSOL < Start)) ;
                is_lick = sum(tSOL > Start & tSOL < End) >= p.n_licks ;
                if prelick
                    obj.trial_types(trial) = 7 ;
                elseif (Sound(trial) == 1) && is_lick
                    obj.trial_types(trial) = 1 ;
                elseif (Sound(trial) == 0) && ~is_lick
                    obj.trial_types(trial) = 2 ;
                elseif (Sound(trial) == 1) && ~is_lick
                    obj.trial_types(trial) = 3 ;
                elseif (Sound(trial) == 0) && is_lick
                    obj.trial_types(trial) = 4 ;
                    punish_time = obj.Punishment_times(find(obj.Punishment_times > obj.trial_times(trial) , 1 , 'first')) ;
                    if ~isempty(punish_time)
                        obj.is_early_punish(trial) = (punish_time - obj.trial_times(trial)) < obj.EARLY_PUNISH_CRIT ;
                    end
                end
            end            
           
            if isempty(obj.ft)
                obj.ft = 1 ;
            end
            
            if isempty(obj.lt)
                obj.lt = length(obj.trial_types) ;
            end
                                   
        end
                    
        function choose_first_last_trials(obj)

            xx = obj.LICK_RATE_CRIT_TRIALS ;

            lick_rate = [] ;

            for i = (xx+1):(length(obj.trial_types))
                hit = sum(obj.trial_types((i-xx+1) : i)==1) ;
                CR = sum(obj.trial_types((i-xx+1) : i)==2) ;
                miss = sum(obj.trial_types((i-xx+1) : i)==3) ;
                FA = sum(obj.trial_types((i-xx+1) : i)==4) ;
                lick_rate(end+1) = 100*(hit + FA)/xx ;
            end
            
            Start = 1 ;
            Stop = find(lick_rate<=obj.LICK_RATE_CRIT , 1 ,'first') ;
            if isempty(Stop)
                Stop = length(obj.trial_types) ;
            end
            
            obj.calc_first_last_trials(Start,Stop)
    
        end
        
        function calc_first_last_trials(obj,ft,lt)
            
            obj.ft = ft ;
            obj.lt = lt ;
            obj.ft_all = zeros(3,4) ;
            obj.lt_all = zeros(3,4) ;
                        
            for odor = 1 : 3
                for type = 1 : 4
                    obj.ft_all(odor,type) = sum((obj.trial_odors(1:(ft-1)) == odor) &...
                        (obj.trial_types(1:(ft-1)) == type)) + 1 ;
                    obj.lt_all(odor,type) = sum((obj.trial_odors(1:lt) == odor) &...
                        (obj.trial_types(1:lt) == type)) ;                  
                end
            end
                   
        end
        
        function create_neurons(obj)
            
            obj.n_neurons = size(obj.BOTs,2) ;
            obj.Neurons = odor_sound_neuron.empty ;
            
            for neur = 1 : obj.n_neurons
                BOT = obj.BOTs(:,neur) ;
                obj.Neurons(neur) = odor_sound_neuron(BOT,obj.trial_types,...
                    obj.trial_odors,obj.tSOF,obj.tSOS,obj.ft_all,obj.lt_all,obj.Parameters,obj.is_early_punish) ;
            end
            
        end
               
        function calc_statistics(obj)
            
            types = obj.trial_types(obj.ft:obj.lt) ;
            odors = obj.trial_odors(obj.ft:obj.lt) ;

            for odor = 1:3

                hit(odor) = sum(types == 1 & odors == odor) ;
                CR(odor) = sum(types == 2 & odors == odor) ;
                miss(odor) = sum(types == 3 & odors == odor) ;
                FA(odor) = sum(types == 4 & odors == odor) ;
                n_Go(odor) = hit(odor)+miss(odor) ;
                n_NoGo(odor) = CR(odor)+FA(odor) ;

                obj.hit_rate(odor) = 100*hit(odor)/n_Go(odor) ;
                obj.FA_rate(odor) = 100*FA(odor)/n_NoGo(odor) ;
                obj.d_prime(odor) = (norminv(min([obj.hit_rate(odor)/100 (n_Go(odor)-0.5)/n_Go(odor)])) -...
                    norminv(max([obj.FA_rate(odor)/100 0.5/n_NoGo(odor)]))) ;
                obj.bias(odor) = (norminv(min([obj.hit_rate(odor)/100 (n_Go(odor)-0.5)/n_Go(odor)])) +...
                    norminv(max([obj.FA_rate(odor)/100 0.5/n_NoGo(odor)]))) ;;

            end
            for pair = 1:3
                obj.d_prime_exp(pair) = (norminv(min([obj.hit_rate(pair)/100 (n_Go(pair)-0.5)/n_Go(pair)])) -...
                    norminv(max([obj.FA_rate(4-pair)/100 0.5/n_NoGo(4-pair)]))) ;
            end
            
            n_Go(4) = sum(hit)+sum(miss) ;
            n_NoGo(4) = sum(FA)+sum(CR) ;
            obj.hit_rate(4) = 100*sum(hit)/n_Go(4) ;
            obj.FA_rate(4) = 100*sum(FA)/n_NoGo(4) ;          
            obj.d_prime(4) = (norminv(min([obj.hit_rate(4)/100 (n_Go(4)-0.5)/n_Go(4)])) -...
                    norminv(max([obj.FA_rate(4)/100 0.5/n_NoGo(4)]))) ;
            
            n_odor1 = sum(odors==1) ;
            n_odor3 = sum(odors==3) ;
            odor1_lick_rate = 100*(hit(1)+FA(1))/n_odor1 ;
            odor3_lick_rate = 100*(hit(3)+FA(3))/n_odor3 ;
            
            obj.odor_d_prime = (norminv(min([odor1_lick_rate/100 (n_odor1-0.5)/n_odor1])) -...
                    norminv(max([odor3_lick_rate/100 0.5/n_odor3]))) ;
                
            dDI = [] ;

            for neur = 1 : numel(obj.Neurons)
                n = obj.Neurons(neur) ;
                if ~n.invalid && n.responsive
                    dDI(end+1) = n.exp_DI_diff ;
                end
            end
                
            obj.mean_dDI = mean(dDI) ;
            
        end

    end
    
end