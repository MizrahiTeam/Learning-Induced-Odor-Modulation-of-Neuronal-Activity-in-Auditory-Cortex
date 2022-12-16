% ========================== Opening / closing functions =================

function varargout = odor_sound_anlss(varargin)

    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @odor_sound_anlss_OpeningFcn, ...
                       'gui_OutputFcn',  @odor_sound_anlss_OutputFcn, ...
                       'gui_LayoutFcn',  [] , ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
    
end

function odor_sound_anlss_OpeningFcn(hObject, eventdata, handles, varargin)

    handles.output = hObject;
    
    handles.main_folder = 'C:\Users\OWNER\Desktop\Thesis\Paper - Odor mediated expectation modulates neuronal responses in auditory cortex\Odor sound analysis' ;
    
    set(handles.figure1,'windowscrollWheelFcn',@scroll_func)
        
    handles.phys_mice = [1:3 5] ;
    
    handles.prob_colors = {[0.474 0.623 0.796],[0.5 0.5 0.5],[0.976 0.4 0.369]} ;
    handles.exp_colors = {[0 0.659 0.643],[0.5 0.5 0.5],[0.651 0.549 0.89]} ;
  
    handles.data_holder = {} ;                       
    
    guidata(hObject,handles) ;
    
end

function varargout = odor_sound_anlss_OutputFcn(hObject, eventdata, handles) 

    varargout{1} = handles.output;
    
end

% ========================== Main callbacks & functions ==================

function all_mice_anlss_menu_Callback(hObject, eventdata, handles)

    handles = guidata(hObject) ;
    
    Ax = handles.anlss_ax ;

    handles.all_mice_anlss_menu.UserData = 1 ;
    handles.one_mouse_anlss_menu.UserData = 0 ;
    menu = handles.all_mice_anlss_menu ;
    cla(Ax , 'reset')
    axis(Ax , 'xy')
    axis(Ax , 'auto')
    hold(Ax , 'on')
    set(Ax,'Visible','on','Box','off')
    set(Ax, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')

    switch menu.String{menu.Value}
        case 'Bias progress'
            plot_odor_behav_prog(Ax,hObject,'bias')
        case 'FA rate progress'
            plot_odor_behav_prog(Ax,hObject,'FA rate')
        case 'Hit rate progress'
            plot_odor_behav_prog(Ax,hObject,'hit rate')
        case "d' expectation progress"
            plot_odor_behav_prog(Ax,hObject,"d' expectation")
        case 'Effective days bias'
            plot_effective_days_behav(Ax,hObject,'bias')
        case "Effective days d' expectation"
            plot_effective_days_behav(Ax,hObject,"d' expectation")
        case 'Effective days HR'
            plot_effective_days_behav(Ax,hObject,'HR')
        case 'Effective days FAR'
            plot_effective_days_behav(Ax,hObject,'FAR')
        case 'sound vs. odor'
            plot_odor_sound_behav_prog(Ax,hObject)
        case 'Expectation DI (Go neurons)'
            all_neurons_DI(Ax,hObject,'Go neurons','box')
        case 'Expectation DI (NoGo neurons)'
            all_neurons_DI(Ax,hObject,'NoGo neurons','box')
        case 'DI progress'
            plot_DI_progress(Ax,hObject)
        case 'DI diff. progress'
            plot_DI_diff_progress(Ax,hObject)
        case "dDI to d'"
            plot_dDI_to_dprime(Ax,hObject)
        case 'Population mean (baseline)'
            plot_all_population_means(hObject,'baseline')
        case 'Population mean (bias)'
            plot_all_population_means(hObject,'bias')

    end
        

end

function one_mouse_anlss_menu_Callback(hObject, eventdata, handles)

    handles = guidata(hObject) ;
    
    Ax = handles.anlss_ax ;
    
    handles.all_mice_anlss_menu.UserData = 0 ;
    handles.one_mouse_anlss_menu.UserData = 1 ;


    menu = handles.one_mouse_anlss_menu ;
    cla(Ax,'reset')
    axis(Ax , 'xy')
    axis(Ax , 'auto')
    hold(Ax , 'on')
    set(Ax,'Visible','on','Box','off')
    set(Ax, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')

    switch menu.String{menu.Value}
        case 'Neuronal field'
            display_field(Ax,hObject)
        case 'Single session behavior'
            plot_single_sess_behav_all(Ax,hObject)
        case 'Single session HR'
            plot_single_sess_behav_parts(Ax,hObject,'HR')
        case 'Single session FAR'
            plot_single_sess_behav_parts(Ax,hObject,'FAR')
        case 'Single session bias'
            plot_single_sess_behav_parts(Ax,hObject,'bias')
        case "Single session d'"
            plot_single_sess_behav_parts(Ax,hObject,"d' expectation")
        case 'All sessions bias'
            plot_all_sess_behav(Ax,hObject,'bias')
        case 'All sessions FA rate'
            plot_all_sess_behav(Ax,hObject,'FA rate')                    
        case 'All sessions hit rate'
            plot_all_sess_behav(Ax,hObject,'hit rate')
        case "All sessions d'"
            plot_all_sess_behav(Ax,hObject,'latency')

    end
    
end

function update_anlss_menu(hObject)

    handles = guidata(hObject) ;

    handles.all_mice_anlss_menu.String = {'Bias progress',...
        'FA rate progress','Hit rate progress',...
        "d' expectation progress",'sound vs. odor'...
        'Effective days bias',"Effective days d' expectation",...
        'Effective days HR','Effective days FAR',...
        'Expectation DI (Go neurons)','Expectation DI (NoGo neurons)',...
        'DI progress','DI diff. progress',"dDI to d'",...
        'Population mean (bias)','Population mean (baseline)'} ;

    handles.one_mouse_anlss_menu.String = {'Neuronal field' , 'Single session behavior','Single session HR',...
        'Single session FAR','Single session bias',"Single session d'",'All sessions bias',...
        'All sessions FA rate','All sessions hit rate'} ;
    
    handles.all_mice_anlss_menu.Value = 1 ;
    handles.one_mouse_anlss_menu.Value = 1 ;
        
end

function neuron_filter_menu_Callback(hObject, eventdata, handles)

    handles = guidata(hObject) ;

    mouse = handles.mouse_menu.Value ;
    sess = handles.session_menu.Value ;
    h = handles.data_holder{mouse}(sess) ;

    handles.neuron_menu.String = {} ;
    for neur = 1  : h.n_neurons
        neuron_crit = 0 ; % Filter criteria
        n = h.Neurons(neur) ;
        switch handles.neuron_filter_menu.String{handles.neuron_filter_menu.Value}
            case 'Responsive neurons'
                if ~n.invalid && n.responsive
                    neuron_crit = 1 ;
                end
            case 'Go neurons'
                if ~n.invalid && n.responsive && n.pref_sound==1
                    neuron_crit = 1 ;
                end
            case 'NoGo neurons'
                if ~n.invalid && n.responsive && n.pref_sound==0
                    neuron_crit = 1 ;
                end

        end
        if neuron_crit
            handles.neuron_menu.String{end+1} = ['Neuron ' num2str(neur)] ;
        end
    end
    if isempty(handles.neuron_menu.String)
        handles.neuron_menu.String = {'No neurons'} ;
    end
    
    handles.neuron_menu.Value = 1 ;
    
    neuron_menu_Callback(hObject,0,0)

end

function update_neuron_filter_menu(hObject)

    handles = guidata(hObject) ;

    handles.neuron_filter_menu.String = {'Responsive neurons','Go neurons','NoGo neurons'} ;
    handles.neuron_filter_menu.Value = 1 ;
    
end

% ========================== Single mouse analysis =======================

function display_field(Ax,hObject)

    handles = guidata(hObject) ;
    sess = handles.session_menu.Value ;
    mouse = handles.mouse_menu.Value ;
    h = handles.data_holder{mouse}(sess) ;
    
    if ~isempty(h.field)
        legend(Ax,'hide')

        set(Ax,'YLimMode','auto','XLimMode','auto')
        xlabel(Ax,'')
        ylabel(Ax,'')

        col_hmask=sum(h.mask,3); % Collapsed hmask
        ilu_hmask=[zeros(size(col_hmask,1),1) diff(col_hmask,2,2) zeros(size(col_hmask,1),1)]+[zeros(1,size(col_hmask,2)) ; diff(col_hmask,2,1) ; zeros(1,size(col_hmask,2))]; % illustrated hmask

        field_with_mask = h.field;
        field_with_mask(:,:,3)=ilu_hmask;
        imshow(field_with_mask,'Parent',Ax)
        axis(Ax,'square')

        for i = 1:size(h.mask,3)
            temp = regionprops(h.mask(:,:,i),'centroid') ;
            text(Ax,temp.Centroid(1),temp.Centroid(2),num2str(i),'fontsize',14,'color','blue')
        end
    end
    
end

function plot_all_sess_behav(Ax,hObject,type)

    handles = guidata(hObject) ;
    mouse = handles.mouse_menu.Value ;
    hh = handles.data_holder{mouse} ;
    
    Colors = {[0 0 1],[0.5 0.5 0.5],[1 0 0]} ;
    days = {'Baseline'} ;
    for i = 2 : numel(hh)
       days{i} = ['Bias ' num2str(i-1)] ;
    end
    
    for j = 1 : numel(hh)
        h = hh(j) ;
        trial_types = h.trial_types(h.ft:h.lt) ;
        trial_odors = h.trial_odors(h.ft:h.lt) ;

        for i = 1:3
            hit = sum(trial_types == 1 & trial_odors == i) ;
            CR = sum(trial_types == 2 & trial_odors == i) ;
            miss = sum(trial_types == 3 & trial_odors == i) ;
            FA = sum(trial_types == 4 & trial_odors == i) ;
            n_Go = hit+miss ;
            n_NoGo = CR+FA ;

            accuracy(i,j) = 100*(hit + CR)/(n_Go + n_NoGo) ;
            hit_rate(i,j) = 100*hit/n_Go ;
            FA_rate(i,j) = 100*FA/n_NoGo ;
            d_prime(i,j) = (norminv(min([hit_rate(i,j)/100 (n_Go-0.5)/n_Go])) - norminv(max([FA_rate(i,j)/100 0.5/n_NoGo]))) ;
            bias(i,j) = (norminv(min([hit_rate(i,j)/100 (n_Go-0.5)/n_Go])) + norminv(max([FA_rate(i,j)/100 0.5/n_NoGo]))) ;
            latency(i,j) = mean(h.trial_latencies(trial_types == 1 & trial_odors == i)) ;

        end
    end
    switch type
        case 'bias'
            behav = bias ;
            Title = 'Bias' ;
            YLIM = [-1 2] ;
        case 'FA rate'
            behav = FA_rate ;
            Title = 'FA rate' ;
            YLIM = [-50 50] ;
        case 'hit rate'
            behav = hit_rate ;
            Title = 'hit rate' ;
            YLIM = [-50 50] ;
    end
    plot(Ax,[0 (numel(hh)+1)],[0 0],'k')
    for i = [1 3]
        plot(Ax,behav(i,:)-behav(2,:),'Color',Colors{i})
    end
%     plot(Ax,bias(1,:)-bias(3,:),'k')    

    ylim(Ax , YLIM) 
    set(Ax,'XTick',1:numel(hh),'XTickLabel',days)
    ylabel(Ax,Title,'fontsize',16)
    
end

function plot_single_sess_behav_all(Ax,hObject)

    handles = guidata(hObject) ;
    legend(Ax,'hide')
    sess = handles.session_menu.Value ;
    mouse = handles.mouse_menu.Value ;
    h = handles.data_holder{mouse}(sess) ; 
    
    hold(handles.anlss_ax , 'on')
    
    ylim2 = [-0.3 3] ;
    ylim1 = ylim2*110/ylim2(2) ;
    colororder(Ax,{'k','k'})
    
    hit_rate = h.hit_rate ;
    FA_rate = h.FA_rate ;
    bias = h.bias ;
    d_prime = h.d_prime_exp ;
    
    for odor = 1:3        
        bar(Ax, [0 4] + odor ,[hit_rate(odor) FA_rate(odor)],'FaceColor' , handles.prob_colors{odor},'BarWidth',0.1)
    end 
    
    ylim(Ax , ylim1)
    ylabel(Ax,'Rate (%)','fontsize',16)
    
    plot(Ax,[8 8],ylim1,'--k','Linewidth',2)
    
    yyaxis(Ax,'right')
        
    for odor = 1:3
        
        bar(Ax, 8 + odor ,bias(odor),'FaceColor' , handles.prob_colors{odor},'BarWidth',0.5)
        bar(Ax, 12 + odor ,d_prime(odor),'FaceColor' , handles.exp_colors{odor},'BarWidth',0.5)
    end 

    ylim(Ax , ylim2)
    ylabel(Ax,"d'/Bias (AU)",'fontsize',16)
     
    set(Ax,'XTick',(0:4:16) + 2,'XTickLabel',{'Hit rate','FA rate','Bias',"d'"})
    
end

function plot_single_sess_behav_parts(Ax,hObject,type)

    handles = guidata(hObject) ;
    legend(Ax,'hide')
    sess = handles.session_menu.Value ;
    mouse = handles.mouse_menu.Value ;
    h = handles.data_holder{mouse}(sess) ; 
    
    hold(handles.anlss_ax , 'on')
        
    for odor = 1 : 3
        switch type
            case 'HR'
                hit_rate = h.hit_rate ;
                bar(Ax, odor ,hit_rate(odor),'FaceColor' , handles.prob_colors{odor},'BarWidth',0.8)
                YLIM = [0 100] ;
                Title = 'Hit rate (%)' ;
                yticks = [0 50 100] ;
                ylabels = {'0','50','100'} ;
            case 'FAR'
                FA_rate = h.FA_rate ;
                bar(Ax, odor ,FA_rate(odor),'FaceColor' , handles.prob_colors{odor},'BarWidth',0.8)
                YLIM = [0 100] ;
                Title = 'FA Rate (%)' ;
                yticks = [0 50 100] ;
                ylabels = {'0','50','100'} ;
            case 'bias'
                bias = h.bias ;
                bar(Ax, odor ,bias(odor),'FaceColor' , handles.prob_colors{odor},'BarWidth',0.8)
                YLIM = [0 2] ;
                Title = 'Bias' ;
                yticks = [0 1 2] ;
                ylabels = {'0','001','002'} ;
            case "d' expectation"
                d_prime = h.d_prime_exp ;
                bar(Ax, odor ,d_prime(odor),'FaceColor' , handles.exp_colors{odor},'BarWidth',0.8)
                YLIM = [0 2] ;
                Title = "d'" ;
                yticks = [0 1 2] ;
                ylabels = {'0','001','002'} ;
        end
    end
        
    
    axis(Ax , [0.3 3.7 YLIM])
%     ylabel(Ax,Title,'fontsize',16)     
    set(Ax,'XTick',[],'YTick',yticks,'YTickLabel',ylabels)
    
end

% ========================== All mice analysis ===========================

function all_neurons_DI(Ax,hObject,neur_type,view)

    handles = guidata(hObject) ;
    h = handles.data_holder ;
    
    y_values = [] ;
    x_values = [] ;
    ids = [] ;
    sig_signs = {'n.s.','*','**','***','****','*****'} ;
    pp = [inf,0.05,0.005,0.0005,0.00005,0.000005] ;
    
    mice = handles.phys_mice ;
    sessions = 3:5 ;
    
    switch neur_type
        case 'Go neurons'
            pref_sound = 1 ;
        case 'NoGo neurons'
            pref_sound = 0 ;
    end
    
    for mouse = mice

        neurons = 1 : h{mouse}(1).n_neurons ;
        for sess = sessions
            hh = h{mouse}(sess) ;
            for neur = neurons
                n = hh.Neurons(neur) ;
                if ~n.invalid && n.responsive && sum(n.pref_sound == pref_sound)
                    neurons(neurons == neur) = [] ;
                    y_values = [y_values ; n.exp_DI] ;
                    xloc = 1:3 ;
                    x_values = [x_values ; xloc] ;
                    YLIM = [0.4 1.2] ;
                    yticks = 0.5:0.1:1 ;
                    Title = 'DI' ;
                    XLabels = {'Expected','Neutral','Unexpected'} ;
                    ids(end+1,:) = [mouse sess neur] ;
                end
            end
        end
    end
    clc
    handles.h_scat = {} ;
    for neur = 1 : size(y_values,1)
        handles.h_scat{neur} = scatter(Ax,...
            x_values(neur,:)-0.2+0.4*rand(size(x_values(neur,:))),y_values(neur,:),20,...
            'MarkerEdgeColor','k','MarkerFaceColor','k') ;
    end

    for i = 1 : 3
        Data{i} = y_values(x_values==i) ;
        plot(Ax,i,mean(y_values(x_values==i)),'+r','markersize',20)
    end
    
    p(1) = signrank(Data{1},Data{2}) ;
    text(Ax,1.2,1.05,sig_signs{sum(3*p(1)<=pp)},'fontsize',16)
    plot(Ax,[1. 1.9],[1 1],'k')
    p(2) = signrank(Data{2},Data{3}) ;
    text(Ax,2.2,1.05,sig_signs{sum(3*p(2)<=pp)},'fontsize',16)
    plot(Ax,[2.1 2.9],[1 1],'k')
    p(3) = signrank(Data{1},Data{3}) ;
    text(Ax,1.7,1.15,sig_signs{sum(3*p(3)<=pp)},'fontsize',16)
    plot(Ax,[1.1 2.9],[1.1 1.1],'k')

    axis(Ax , [0 4 YLIM]) 
    set(Ax,'XTick',xloc,'XTickLabel',XLabels,'YTick',yticks,'fontsize',16)
    ylabel(Ax,Title,'fontsize',16)
    
    if strcmp(view,'box')
        figure
        boxplot(y_values)
        hold on
        text(gca,1.2,1.05,sig_signs{sum(3*p(1)<=pp)},'fontsize',16)
        plot(gca,[1. 1.9],[1 1],'k')
        text(gca,2.2,1.05,sig_signs{sum(3*p(2)<=pp)},'fontsize',16)
        plot(gca,[2.1 2.9],[1 1],'k')
        text(gca,1.7,1.15,sig_signs{sum(3*p(3)<=pp)},'fontsize',16)
        plot(gca,[1.1 2.9],[1.1 1.1],'k')
        axis(gca , [0 4 YLIM]) 
        set(gca,'XTick',[],'YTick',yticks,'fontsize',16)
%         ylabel(gca,Title,'fontsize',16)
    end
    
    guidata(hObject,handles)

end

function plot_odor_behav_prog(Ax,hObject,type)

    handles = guidata(hObject) ;
    
    Colors = handles.prob_colors ;
    sig_signs = {'n.s.','*','**','***','****','*****'} ;
    pp = [1,0.05,0.005,0.0005,0.00005,0.000005] ;
    shapes = {'+','s','d','v','^'} ;
    
    for mouse = 1 : numel(handles.data_holder)
        n_days(mouse) = numel(handles.data_holder{mouse}) ;
    end
    
    behav = cell(1,max(n_days)) ;
    
    for mouse = 1 : numel(handles.data_holder)
        hh = handles.data_holder{mouse} ;
        for day = 1 : n_days(mouse)
            h = hh(day) ;               
            switch type
                case 'bias'
                    behav{day}(:,end+1) = h.bias ;
                    Title = '\DeltaBias' ;
                    YLIM = [-1.2 1.5] ;
                    yticks = [-1 0 1] ;
%                     ylabels = {'-01','0','01'} ;
                case 'latency'
                    behav{day}(:,end+1) = h.hit_latency ;
                    Title = '\DeltaLatency (sec)' ;
                    YLIM = [-0.1 0.4] ;
                case 'FA rate'
                    behav{day}(:,end+1) = h.FA_rate(1:3) ;
                    Title = '\DeltaFA rate' ;
                    YLIM = [-45 45] ;
                    yticks = [-40 -20 0 20 40] ;
%                     ylabels = {num2str(YLIM(1)),'0',num2str(YLIM(2))} ;
                case 'hit rate'
                    behav{day}(:,end+1) = h.hit_rate(1:3) ;
                    Title = '\DeltaHit rate' ;
                    YLIM = [-30 30] ;
                    yticks = [YLIM(1) 0 YLIM(2)] ;
%                     ylabels = {num2str(YLIM(1)),'0',num2str(YLIM(2))} ;
                case "d' expectation"
                    behav{day}(:,end+1) = h.d_prime_exp ;
                    Title = "\Deltad'" ;
                    YLIM = [-1.5 1.5] ;
                    yticks = [-1 0 1] ;
%                     ylabels = {'-01','0','01'} ;
                    Colors = handles.exp_colors ;

            end
        end
    end
    
    for day = 1 : numel(behav)
        d_behav{day} = behav{day} - repmat(behav{day}(2,:),3,1) ;
        m_behav(:,day) = mean(d_behav{day},2) ;
        s_behav(:,day) = std(d_behav{day},[],2)/sqrt(size(behav{day},2)) ;
    end
    
    for mouse = 1 : numel(handles.data_holder)
        for odor = [1 3]
            trace = [] ;
           for day = 1 : numel(behav)
               trace(end+1) = d_behav{day}(odor,mouse) ;
           end
           plot(Ax,(1:length(trace)) - 0.1,trace,shapes{mouse},'Color',...
               Colors{odor} + 0.5*(ones(1,3)-Colors{odor}),'LineWidth',0.2,'MarkerSize',4)
        end
    end
    
    for odor = [1 3]
        errorbar(Ax,m_behav(odor,:),s_behav(odor,:),'Color',Colors{odor},'LineWidth',1.5)
    end
    
    plot(Ax,0:numel(behav)+1,zeros(1,numel(behav)+2),'k')
    
    for day = 1 : numel(d_behav)
        odor1 = d_behav{day}(1,:) ;
        odor2 = d_behav{day}(2,:) ;
        odor3 = d_behav{day}(3,:) ;
        
        p = signrank(odor1,odor3,'Tail','Right') ;
        if ~ isnan(p)
            text(Ax,day,0.9*YLIM(2),sig_signs{sum(p<=pp)},'fontsize',16)
        end
        
    end
    

    days = {'PBS'} ;
    for i = 1 : numel(behav)
        days{i+1} = ['BS ' num2str(i)] ;
    end
    axis(Ax , [0.5 7.5 YLIM]) 
    set(Ax,'XTick',1:numel(behav),'XTickLabel',{},'YTick',yticks)
%     ylabel(Ax,Title,'fontsize',16)
    
end

function plot_odor_sound_behav_prog(Ax,hObject)

    handles = guidata(hObject) ;
    
    sig_signs = {'n.s.','*','**','***','****','*****'} ;
    pp = [1,0.05,0.005,0.0005,0.00005,0.000005] ;
    YLIM = [-1 2.5] ;
       
    for mouse = 1 : numel(handles.data_holder)
        n_days(mouse) = numel(handles.data_holder{mouse}) ;
    end
       
    odor_dprime = cell(1,max(n_days)) ;
    sound_dprime = cell(1,max(n_days)) ;    
        
    for mouse = 1 : numel(handles.data_holder)
        hh = handles.data_holder{mouse} ;
        for day = 1 : numel(hh)           
            h = hh(day) ;

            sound_dprime{day}(end+1) = h.d_prime(4) ;            
            odor_dprime{day}(end+1) = h.odor_d_prime ;
        end
    end
    
 
    for day = 1 : numel(sound_dprime)
        m_sound_behav(day) = mean(sound_dprime{day}) ;
        s_sound_behav(day) = std(sound_dprime{day})/sqrt(numel(sound_dprime{day})) ;
        m_odor_behav(day) = mean(odor_dprime{day}) ;
        s_odor_behav(day) = std(odor_dprime{day})/sqrt(numel(odor_dprime{day})) ;
    end
   
    errorbar(Ax,m_sound_behav,s_sound_behav,'k') 
    errorbar(Ax,m_odor_behav,s_odor_behav,'--k')
    plot(Ax,[0 numel(m_sound_behav)+1],[0 0],'k')
    

    
    for day = 1 : numel(m_sound_behav)
        
        p(day) = signrank(sound_dprime{day},odor_dprime{day},'Tail','Right') ;
        if ~ isnan(p(day))
            text(Ax,day,0.9*YLIM(2),sig_signs{sum(p(day)<=pp)},'fontsize',16)
        end
    end
    

    days = {'PBS'} ;
    for i = 1 : numel(m_sound_behav)
        days{i+1} = ['BS ' num2str(i)] ;
    end
    axis(Ax , [0 max(n_days)+1 YLIM]) 
    set(Ax,'XTick',1:numel(m_sound_behav),'XTickLabel',days)
    set(Ax,'XTick',1:numel(m_sound_behav),'XTickLabel',{})
    ylabel(Ax,"d'",'fontsize',16)
    
end

function plot_effective_days_behav(Ax,hObject,type)

    handles = guidata(hObject) ;
    
    first_day = 3 ;
    last_day = 5 ;
    
    
    Colors = handles.prob_colors ;
    sig_signs = {'n.s.','*','**','***','****','*****'} ;
    pp = [inf,0.05,0.005,0.0005,0.00005,0.000005] ;
    
    for mouse = 1 : numel(handles.data_holder)
        n_days(mouse) = numel(handles.data_holder{mouse}) ;
    end
    
    behav = cell(1,max(n_days)) ;
    
    for mouse = 1 : numel(handles.data_holder)
        hh = handles.data_holder{mouse} ;
        for day = 1 : n_days(mouse)
            h = hh(day) ;

            switch type
                case 'bias'
                    behav{day}(:,end+1) = h.bias(1:3) ;
                    Title = 'Bias' ;
                    YLIM = [-2 5] ;
                    yticks = -2:2:4 ;
                    ylabels = {'-02','000','002','004'} ;
                case 'latency'
                    behav{day}(:,end+1) = h.hit_latency(1:3) ;
                    Title = 'Latency (sec)' ;
                    YLIM = [-0.1 0.4] ;
                case 'FAR'
                    behav{day}(:,end+1) = h.FA_rate(1:3) ;
                    Title = 'FA rate' ;
                    YLIM = [0 110] ;
                    yticks = [0 50 100] ;
                    ylabels = {'0','50','100'} ;
                case 'HR'
                    behav{day}(:,end+1) = h.hit_rate(1:3) ;
                    Title = 'Hit rate' ;
                    YLIM = [0 110] ;
                    yticks = [0 50 100] ;
                    ylabels = {'0','50','100'} ;
                case "d' expectation"
                    behav{day}(:,end+1) = h.d_prime_exp(1:3) ;
                    Title = "d'" ;
                    YLIM = [0 4] ;
                    yticks = 0:2:4 ;
                    ylabels = {'000','002','004'} ;
                    Colors = handles.exp_colors ;

            end
        end
    end
       
    d_behav = [] ;
    for day = first_day : last_day
        d_behav = [d_behav behav{day}] ;
    end
    m_behav = mean(d_behav,2) ;
    s_behav = std(d_behav,[],2)/sqrt(size(behav,2)) ;
    
    for odor = 1:3
        scatter(Ax,odor*ones(1,size(d_behav,2)),d_behav(odor,:),100,Colors{odor},'filled')
%         errorbar(Ax,odor,m_behav(odor),s_behav(odor),'.','Color','k','LineWidth',7,'CapSize',15)
    end
    for sess = 1 : size(d_behav,2)
        plot(Ax,1:3,d_behav(:,sess),'k')
    end
    
    [~,p,~,stats] = ttest(d_behav(1,:),d_behav(3,:),'Tail','Right') ;
    if ~ isnan(p)
        text(Ax,2,0.9*YLIM(2),sig_signs{sum(3*p<=pp)},'fontsize',16)
        plot(Ax,[1 3],0.8*YLIM(2)*[1 1],'k')
    end


    [~,p,~,stats] = ttest(d_behav(1,:),d_behav(2,:),'Tail','Right') ;
    if ~ isnan(p)
        text(Ax,1.5,0.7*YLIM(2), sig_signs{sum(3*p<=pp)},'fontsize',16)
        plot(Ax,[1.1 1.9],0.6*YLIM(2)*[1 1],'k')
    end

    
    [~,p,~,stats] = ttest(d_behav(3,:),d_behav(2,:),'Tail','Left') ;
    if ~ isnan(p)
        text(Ax,2.5,0.7*YLIM(2), sig_signs{sum(3*p<=pp)},'fontsize',16)
        plot(Ax,[2.1 2.9],0.6*YLIM(2)*[1 1],'k')
    end


    axis(Ax , [0.5 3.5 YLIM]) 
    set(Ax,'XTick',[],'YTick',yticks,'YTickLabel',ylabels)
    
end

function plot_DI_progress(Ax,hObject)

    handles = guidata(hObject) ;
    h = handles.data_holder ;
    
    sig_signs = {'n.s.','*','**','***','****','*****'} ;
    pp = [inf,0.05,0.005,0.0005,0.00005,0.000005] ;
    
    mice = handles.phys_mice ;
    sessions = 1 : 7 ;
    DI = cell(1,length(sessions)) ;
    
    for mouse = mice
        for sess = sessions
            hh = h{mouse}(sess) ; 
            for neur = 1 : hh.n_neurons
                n = hh.Neurons(neur) ;
                if ~n.invalid && n.responsive
                    DI{sess} = [DI{sess} ; n.exp_DI([1 3])] ;
                end
            end
        end
    end
    
    for sess = sessions
        m_DI(sess,:) = mean(DI{sess},1) ;
        s_DI(sess,:) = std(DI{sess},[],1)/sqrt(size(DI{sess},1)) ;
        p = signrank(DI{sess}(:,1),DI{sess}(:,2)) ;
        p_text{sess} = sig_signs{sum(p<=pp)} ;
    end
    errorbar(Ax,m_DI(:,1),s_DI(:,1),'Color',handles.exp_colors{1},'LineWidth',2)
    errorbar(Ax,m_DI(:,2),s_DI(:,2),'Color',handles.exp_colors{3},'LineWidth',2)
    text(Ax,sessions,0.67*ones(size(sessions)),p_text,'fontsize',16)
    
    axis(Ax,[0 length(sessions)+1 0.55 0.7])
    days = {'PBS'} ;
    for i = 1 : numel(sessions)
        days{i+1} = ['BS ' num2str(i)] ;
    end
%     set(Ax,'Xtick',1:length(sessions),'XTickLabel',days)
    set(Ax,'Xtick',1:length(sessions),'XTickLabel',{},'YTick',[0.55 0.6 0.65 0.7])
%     ylabel(Ax,'DI','Fontsize',16)
        
        
    guidata(hObject,handles)
    
end

function plot_DI_diff_progress(Ax,hObject)

    handles = guidata(hObject) ;
    h = handles.data_holder ;
    
    sig_signs = {'n.s.','*','**','***','****','*****'} ;
    pp = [inf,0.05,0.005,0.0005,0.00005,0.000005] ;
    
    mice = handles.phys_mice ;
    sessions = 1 : 7 ;
    dDI = cell(1,length(sessions)) ;
    
    for mouse = mice
        for sess = sessions
            hh = h{mouse}(sess) ; 
            for neur = 1 : hh.n_neurons
                n = hh.Neurons(neur) ;
                if ~n.invalid && n.responsive
                    dDI{sess}(end+1) = -diff(n.exp_DI([1 3])) ;
                end
            end
        end
    end
    
    for sess = sessions
        m_dDI(sess,:) = mean(dDI{sess}) ;
        s_dDI(sess,:) = std(dDI{sess})/sqrt(length(dDI{sess})) ;
        p = signrank(dDI{sess},0) ;
        p_text{sess} = sig_signs{sum(p<=pp)} ;
    end
    plot(Ax,[0 8],[0 0],'--k')
    errorbar(Ax,m_dDI,s_dDI,'Color','k','LineWidth',2)
    text(Ax,sessions,0.04*ones(size(sessions)),p_text,'fontsize',16)
    
    axis(Ax,[0 length(sessions)+1 -0.05 0.02])
    days = {'PBS'} ;
    for i = 1 : numel(sessions)
        days{i+1} = ['BS ' num2str(i)] ;
    end
%     set(Ax,'Xtick',1:length(sessions),'XTickLabel',days)
    set(Ax,'Xtick',1:length(sessions),'XTickLabel',{},'YTick',[-0.05 -0.025 0])
%     ylabel(Ax,'DI','Fontsize',16)
        
        
    guidata(hObject,handles)
    
end

function plot_all_population_means(hObject,time)

    plot_population_mean(hObject,time,1)
    plot_population_mean(hObject,time,0)

end

function plot_population_mean(hObject,time,neuron_type)

    handles = guidata(hObject) ;
    h = handles.data_holder ;
    
    if neuron_type
        Title = 'Go neurons' ;
    else
        Title = 'NoGo neurons' ;
    end
        
    sig_signs = {'n.s.','*','**','***'} ;
    pp = [inf,0.05,0.005,0.0005] ;
    colors = handles.prob_colors ;
        
    mice = handles.phys_mice ;
    
    traces = cell(2,3) ;
    resps = [] ;
    
    temp = {} ;
    
    for mouse = mice
        switch time
            case 'baseline'
                sessions = 1 ;
            case 'bias'
                sessions = 3:5 ;
        end

        neurons = 1 : h{mouse}(1).n_neurons ;
        for sess = sessions
            hh = h{mouse}(sess) ;
            for neur = neurons
                n = hh.Neurons(neur) ;
                if n.responsive && n.pref_sound==neuron_type
                    neurons(neurons == neur) = [] ;
                    temp{end+1} = cell(2,3) ;
                    resps(end+1,:,:) = [n.mean_resps(:,1:2) mean(n.mean_resps(:,1:2),2)] ;
                    for type = 1:2
                        for odor = 1 : 3
                            traces{type,odor} = [traces{type,odor} ; mean(n.get_dFF(odor,type),1)] ;
                        end
                    end
                end
            end
        end

    end 
    
    figure()
    
    y_lim = 0.2 ;
    for odor = 1:3
        x = hh.Time ;
        y1 = mean(traces{1,odor},1) ;
        dy1 = std(traces{1,odor},[],1)/sqrt(size(traces{1,odor},1)) ;
        y2 = mean(traces{2,odor},1) ;
        dy2 = std(traces{2,odor},[],1)/sqrt(size(traces{2,odor},1)) ;
        y1 = y1(x>-1 & x<2) ;
        dy1 = dy1(x>-1 & x<2) ;
        y2 = y2(x>-1 & x<2) ;
        dy2 = dy2(x>-1 & x<2) ;
        x = x(x>-1 & x<2) ;
        
%         y_lim = max([y_lim 1.2*max(y1) 1.2*max(y2)]) ;
        
        subplot(3,2,1)
        if odor ~=2
            plot(x, y1 , 'color' , colors{odor} , 'linewidth' , 2)
            patch([x flip(x)] , [y1+dy1 flip(y1-dy1)]  , colors{odor} , 'facealpha' , 0.2)
            axis([-1 2 -0.04 0.11])
            hold on
        end
        
        set(gca,'XTick',[],'YTick',[],'box', 'off','xcolor','none','ycolor','none')
        plot(gca,[-0.6 -0.6],[0.07 0.09],'k','linewidth',2)
        plot(gca,[-0.6 0.4],[0.07 0.07],'k','linewidth',2)
        ht = text(gca,-0.85,0.07,'0.02 dF/F','fontsize',10) ;
        text(gca,-0.75,0.065,'1 sec','fontsize',10)
        set(ht , 'rotation' , 90)
        rectangle(gca,'Position',[-0.5 -0.015 0.5 0.005],'FaceColor','g')
        rectangle(gca,'Position',[0 -0.015 0.15 0.005],'FaceColor',(type==1)*ones(1,3))
        
        subplot(3,2,[3 5])
        hold on
        scatter(odor - 0.3 + 0.6*rand(1,size(resps,1)), min(resps(:,odor,1), y_lim), 6, colors{odor},'filled')
        plot(odor + [-0.3 0.3],mean(resps(:,odor,1))*[1 1], 'k', 'LineWidth', 2)
        set(gca,'XTick',[])
        axis([0.5 3.5 0 y_lim])     
        
        subplot(3,2,2)
%         title(Title,'FontSize',16)
        if odor ~=2
            plot(x, y2 , 'color' , colors{odor} , 'linewidth' , 2)
            patch([x flip(x)] , [y2+dy2 flip(y2-dy2)]  , colors{odor} , 'facealpha' , 0.2)
            axis([-1 2 -0.04 0.11])
            hold on
        end
        
        set(gca,'XTick',[],'YTick',[],'box', 'off','xcolor','none','ycolor','none')
        rectangle(gca,'Position',[-0.5 -0.015 0.5 0.005],'FaceColor','g')
        rectangle(gca,'Position',[0 -0.015 0.15 0.005],'FaceColor',(type==1)*ones(1,3))
        
        subplot(3,2,[4 6])
        hold on
        scatter(odor - 0.3 + 0.6*rand(1,size(resps,1)), min(resps(:,odor,2), y_lim), 6, colors{odor},'filled')
        plot(odor + [-0.3 0.3],mean(resps(:,odor,2))*[1 1], 'k', 'LineWidth', 2)
        set(gca,'XTick',[])
        axis([0.5 3.5 0 y_lim])
        
    end
    
    mcc = 3 ;
    subplot(3,2,[3 5])
    [~,p] = ttest(resps(:,1,1),resps(:,2,1)) ;
    text(1.3,0.088,sig_signs{sum(mcc*p<=pp)},'fontsize',12)
    plot([0.9 1.9],0.077*[1 1],'k')
    
    [~,p] = ttest(resps(:,2,1),resps(:,3,1)) ;
    text(2.3,0.088,sig_signs{sum(mcc*p<=pp)},'fontsize',12)
    plot([2.1 3.1],0.077*[1 1],'k')
    
    [~,p] = ttest(resps(:,1,1),resps(:,3,1)) ;
    text(1.8,0.11,sig_signs{sum(mcc*p<=pp)},'fontsize',12)
    plot([1 3],0.1*[1 1],'k')
    
    subplot(3,2,[4 6])
    [~,p] = ttest(resps(:,1,2),resps(:,2,2)) ;
    text(1.3,0.088,sig_signs{sum(mcc*p<=pp)},'fontsize',12)
    plot([0.9 1.9],0.077*[1 1],'k')
    
    [~,p] = ttest(resps(:,2,2),resps(:,3,2)) ;
    text(2.3,0.088,sig_signs{sum(mcc*p<=pp)},'fontsize',12)
    plot([2.1 3.1],0.077*[1 1],'k')
    
    [~,p] = ttest(resps(:,1,2),resps(:,3,2)) ;

    text(1.8,0.11,sig_signs{sum(mcc*p<=pp)},'fontsize',12)
    plot([1 3],0.1*[1 1],'k') 
    
end

function plot_dDI_to_dprime(Ax,hObject)

    handles = guidata(hObject) ;
    hh = handles.data_holder ;
    
    lims = [-1 1 -0.1 0.1] ;
    
    mice = handles.phys_mice ; 
    
    d_d_prime = [] ;
    dDI = [] ;
    handles.h_scat = {} ;
    
    for mouse = 1 : length(mice)
        for sess = 1 : (numel(hh{mice(mouse)}))
            h = hh{mice(mouse)}(sess) ;
            if ~isnan(h.mean_dDI(1))
                d_d_prime(end+1) = (h.d_prime_exp(1) - h.d_prime_exp(3))/(h.d_prime_exp(1) + h.d_prime_exp(3)) ;
                dDI(end+1) = -h.mean_dDI ;
                handles.h_scat{mouse}{sess} = scatter(Ax,d_d_prime(end),dDI(end),40,...
                    'MarkerEdgeColor','k','MarkerFaceColor','k') ;
            end
        end
    end
    
    
    [R,p] = corrcoef(d_d_prime,dDI) ;
    
    plot(Ax,lims(1:2),[0 0],'k')
    plot(Ax,[0 0],lims(3:4),'k')
    
    par = polyfit(d_d_prime,dDI,1) ;
    
    plot(Ax,[-1 1],polyval(par,[-1.45 1.45]) , 'r','LineWidth',2)
    
%     xlabel(Ax,"d' CI",'FontSize',16)
%     ylabel(Ax,"<\DeltaDI>",'FontSize',16)

    set(Ax,'XTick',[lims(1) 0 lims(2)],'YTick',[lims(3) 0 lims(4)])
    
    guidata(hObject,handles)

end

% ========================== Helper functions ============================

function plot_dFF(hObject)

    handles = guidata(hObject) ;
    
    delete(get(handles.neuron_panel,'children'))
    
    h_fig = handles.neuron_panel ;

    if strcmp(handles.neuron_menu.String{handles.neuron_menu.Value},'No neurons')
        return
    end
    
    neur = str2double(handles.neuron_menu.String{handles.neuron_menu.Value}(8:end)) ;    
    sess = handles.session_menu.Value ;
    mouse = handles.mouse_menu.Value ;
    h = handles.data_holder{mouse}(sess) ;
    if isnan(neur)
        return
    else
        n = h.Neurons(neur) ;
    end
    
    odors = handles.odor_menu.Value - 1 ;
    if odors == 0
        odors = 1 : 3 ;
    elseif odors == 4
        odors = [1 3] ;        
    end

    titles = {'Hit' , 'CR'} ;
    colors = handles.prob_colors ;
       
    lims = [] ;
    n_frames = h.frames_before_stim + h.frames_for_amp ;
        
    for type = 1 : 2
        subplot(1,2,type,'Parent',h_fig)
        hold on
        for odor = odors
            dFF = n.get_dFF(odor,type) ;
            if size(dFF,1) > 0
                x = h.Time ;
                y = smoothdata(mean(dFF,1),'movmean',3) ;

                dy = std(dFF,[],1)/sqrt(size(dFF,1)) ;
                plot(x, y , 'color' , colors{odor} , 'linewidth' , 2)
                patch([x flip(x)] , [y+dy flip(y-dy)]  , colors{odor} , 'facealpha' , 0.2)
                if type < 3
                    lims(:,end+1) = [max(y(1:n_frames)) ; min(y(1:n_frames))] ;
                end
            end
        end
        title(titles{type})
        if type == 1
            xlabel('Time (seconds)')
            ylabel('\DeltaF/F_0')
        end

    end

    for type = 1 : 2
        subplot(1,2,type,'Parent',h_fig)
        axis([-1 2 min(lims(2,:))-0.1 max(lims(1,:))+0.1]) ;
        if sum(n.responsiveness(:,type)>0)
            plot(-0.4 , max(lims(1,:)) , '*k')
        end
    end
        
    
end

function scroll_func(hObject,eventdata)

    handles = guidata(hObject) ;
    
    fig_pos = get(hObject,'Position') ;
    mouse_pos = get(hObject,'CurrentPoint') ;
    
    anlss_type_menu_pos(1) = fig_pos(3) * handles.all_mice_anlss_menu.Position(1) ;
    anlss_type_menu_pos(2) = fig_pos(3) * (handles.all_mice_anlss_menu.Position(1) + handles.all_mice_anlss_menu.Position(3)) ;
    anlss_type_menu_pos(3) = fig_pos(4) * handles.all_mice_anlss_menu.Position(2) ;
    anlss_type_menu_pos(4) = fig_pos(4) * (handles.all_mice_anlss_menu.Position(2) + handles.all_mice_anlss_menu.Position(4)) ;
    
    anlss_type_menu_pos2(1) = fig_pos(3) * handles.one_mouse_anlss_menu.Position(1) ;
    anlss_type_menu_pos2(2) = fig_pos(3) * (handles.one_mouse_anlss_menu.Position(1) + handles.one_mouse_anlss_menu.Position(3)) ;
    anlss_type_menu_pos2(3) = fig_pos(4) * handles.one_mouse_anlss_menu.Position(2) ;
    anlss_type_menu_pos2(4) = fig_pos(4) * (handles.one_mouse_anlss_menu.Position(2) + handles.one_mouse_anlss_menu.Position(4)) ;
    
    session_menu_pos(1) = fig_pos(3) * handles.session_menu.Position(1) ;
    session_menu_pos(2) = fig_pos(3) * (handles.session_menu.Position(1) + handles.session_menu.Position(3)) ;
    session_menu_pos(3) = fig_pos(4) * handles.session_menu.Position(2) ;
    session_menu_pos(4) = fig_pos(4) * (handles.session_menu.Position(2) + handles.session_menu.Position(4)) ;
    
    neuron_menu_pos(1) = fig_pos(3) * handles.neuron_menu.Position(1) ;
    neuron_menu_pos(2) = fig_pos(3) * (handles.neuron_menu.Position(1) + handles.neuron_menu.Position(3)) ;
    neuron_menu_pos(3) = fig_pos(4) * handles.neuron_menu.Position(2) ;
    neuron_menu_pos(4) = fig_pos(4) * (handles.neuron_menu.Position(2) + handles.neuron_menu.Position(4)) ;
    
    neuron_filter_menu_pos(1) = fig_pos(3) * handles.neuron_filter_menu.Position(1) ;
    neuron_filter_menu_pos(2) = fig_pos(3) * (handles.neuron_filter_menu.Position(1) + handles.neuron_filter_menu.Position(3)) ;
    neuron_filter_menu_pos(3) = fig_pos(4) * handles.neuron_filter_menu.Position(2) ;
    neuron_filter_menu_pos(4) = fig_pos(4) * (handles.neuron_filter_menu.Position(2) + handles.neuron_filter_menu.Position(4)) ;
    
    if (mouse_pos(1) > anlss_type_menu_pos(1)) && (mouse_pos(1) < anlss_type_menu_pos(2)) &&...
       (mouse_pos(2) > anlss_type_menu_pos(3)) && (mouse_pos(2) < anlss_type_menu_pos(4))
   
        val = handles.all_mice_anlss_menu.Value ;
        new_val = val + eventdata.VerticalScrollCount ;
        max_val = numel(handles.all_mice_anlss_menu.String) ;
   
        handles.all_mice_anlss_menu.Value = min([max([new_val 1]) max_val]) ;
        
        all_mice_anlss_menu_Callback(hObject,0,handles)
        
    elseif (mouse_pos(1) > anlss_type_menu_pos2(1)) && (mouse_pos(1) < anlss_type_menu_pos2(2)) &&...
       (mouse_pos(2) > anlss_type_menu_pos2(3)) && (mouse_pos(2) < anlss_type_menu_pos2(4))
   
        val = handles.one_mouse_anlss_menu.Value ;
        new_val = val + eventdata.VerticalScrollCount ;
        max_val = numel(handles.one_mouse_anlss_menu.String) ;
   
        handles.one_mouse_anlss_menu.Value = min([max([new_val 1]) max_val]) ;
        
        one_mouse_anlss_menu_Callback(hObject,0,handles)
        
    elseif (mouse_pos(1) > session_menu_pos(1)) && (mouse_pos(1) < session_menu_pos(2)) &&...
           (mouse_pos(2) > session_menu_pos(3)) && (mouse_pos(2) < session_menu_pos(4))
   
        val = handles.session_menu.Value ;
        new_val = val + eventdata.VerticalScrollCount ;
        max_val = numel(handles.session_menu.String) ;
   
        handles.session_menu.Value = min([max([new_val 1]) max_val]) ;
        
        session_menu_Callback(hObject,0,handles)
        
    elseif (mouse_pos(1) > neuron_menu_pos(1)) && (mouse_pos(1) < neuron_menu_pos(2)) &&...
           (mouse_pos(2) > neuron_menu_pos(3)) && (mouse_pos(2) < neuron_menu_pos(4))
   
        val = handles.neuron_menu.Value ;
        new_val = val + eventdata.VerticalScrollCount ;
        max_val = numel(handles.neuron_menu.String) ;
   
        handles.neuron_menu.Value = min([max([new_val 1]) max_val]) ;
        
        neuron_menu_Callback(hObject,0,handles)
        
    elseif (mouse_pos(1) > neuron_filter_menu_pos(1)) && (mouse_pos(1) < neuron_filter_menu_pos(2)) &&...
           (mouse_pos(2) > neuron_filter_menu_pos(3)) && (mouse_pos(2) < neuron_filter_menu_pos(4))
   
        val = handles.neuron_filter_menu.Value ;
        new_val = val + eventdata.VerticalScrollCount ;
        max_val = numel(handles.neuron_filter_menu.String) ;
   
        handles.neuron_filter_menu.Value = min([max([new_val 1]) max_val]) ;
        
        neuron_filter_menu_Callback(hObject,0,handles)
               
    end
    
end

function update_session_menu(hObject)

    handles = guidata(hObject) ;
    
    mouse = handles.mouse_menu.Value ;
    h = handles.data_holder{mouse} ;
    
    handles.session_menu.String = {'Baseline'} ;
    for i = 2 : numel(h)
        handles.session_menu.String{i} = ['Bias ' num2str(i-1)] ;
    end
    
end

% ========================== Callbacks ===================================

function load_data_button_Callback(hObject, eventdata, handles)
  

    [data_file,data_path] = uigetfile('*.mat' , 'Choose data file' , handles.main_folder) ;
    if data_file
        
        load([data_path data_file])
        update_anlss_menu(hObject)
        update_neuron_filter_menu(hObject)
        handles.odor_menu.String = {'All odors','Odor 1','Odor 2','Odor 3','Odors 1 & 3'} ;

        handles.data_holder = d ;

        handles.session_menu.Enable = 'on' ;
        handles.neuron_menu.Enable = 'on' ;
        handles.neuron_filter_menu.Enable = 'on' ;
        handles.mouse_menu.Enable = 'on' ;
        handles.all_mice_anlss_menu.Enable = 'on' ;
        handles.one_mouse_anlss_menu.Enable = 'on' ;
        handles.odor_menu.Enable = 'on' ;

        handles.mouse_menu.String = {} ;
        for i = 1 : numel(handles.data_holder)
            handles.mouse_menu.String{end+1} = handles.data_holder{i}(1).Parameters.mouse ;

            p = handles.data_holder{i}(end).Parameters ;
            odor_probs(1) = p.Go_odor1_prob/(p.Go_odor1_prob + p.NoGo_odor1_prob) ;
            odor_probs(2) = p.Go_odor2_prob/(p.Go_odor2_prob + p.NoGo_odor2_prob) ;
            odor_probs(3) = p.Go_odor3_prob/(p.Go_odor3_prob + p.NoGo_odor3_prob) ;
            [~,odor_order] = sort(odor_probs,'descend') ;
            odor_order(4) = 4 ;
            for j = 1 : numel(handles.data_holder{i})
                handles.data_holder{i}(j).odor_order = odor_order ;
                handles.data_holder{i}(j).organize_data()
            end

        end
        handles.mouse_menu.Value = 1 ;
                
        
        guidata(hObject,handles)

        mouse_menu_Callback(hObject,0,0)
        
    end
    
end


function mouse_menu_Callback(hObject, eventdata, handles)

    handles = guidata(hObject) ;
    
    handles.neuron_menu.Value = 1 ;
    update_session_menu(hObject)
    handles.session_menu.Value = 1 ;
    session_menu_Callback(hObject,0,0)
        
    guidata(hObject,handles)
    
end

function neuron_menu_Callback(hObject, eventdata, handles)

    handles = guidata(hObject) ;

    plot_dFF(hObject)
    
end

function odor_menu_Callback(hObject, eventdata, handles)
    neuron_menu_Callback(hObject,0,0)
end

function session_menu_Callback(hObject, eventdata, handles)

    handles = guidata(hObject) ;
    sess = handles.session_menu.Value ;
    mouse = handles.mouse_menu.Value ;
    if sess <= numel(handles.session_menu.String)
        h = handles.data_holder{mouse}(sess) ;

        handles.first_trial.String = num2str(h.ft) ;
        handles.first_trial.Enable = 'on' ;
        handles.last_trial.String = num2str(h.lt) ;
        handles.last_trial.Enable = 'on' ;
        handles.max_trials.String = ['/' num2str(length(h.trial_types))] ;

        neuron_filter_menu_Callback(hObject)

    end

end

% ========================== Creation functions ==========================

function all_mice_anlss_menu_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
end

function mouse_menu_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
end

function neuron_menu_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
end

function neuron_filter_menu_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
end

function odor_menu_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
end

function one_mouse_anlss_menu_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

function session_menu_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
end
