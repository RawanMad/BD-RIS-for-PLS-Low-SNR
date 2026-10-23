        clear; clc;
        rng(3);
        tic;
        
        %% ===================== Parameters =====================
        nMonte = 100;
        Ma = 4;
        Mb = 1;
        Me = 1;
        P = eye(Ma);                      % Kept for compatibility with channel function
        K = 0;                           % Linear Rician K-factor
        NI = 44;                         % Fixed number of RIS elements
        NI_max = NI;
        PT = 15;                         % Transmit power in dBm
        pt = 10^((PT-30)/10);            % Transmit power in watts
        LS = 1;
        
        if LS == 0
            N_B = 1;
        else
            N_B_dBm = -80;
            N_B = 10^((N_B_dBm-30)/10);
        end
        
        kappaE = 1;
        N_E = kappaE*N_B;
        noiseRatio = N_B/N_E;
        
        locB = [150;0];
        locE_list = 115:5:145;             % Eve x-location sweep (m)
        locE = [locE_list(1);0];
        locA = [0;0];
        locI = [145;5];
        
        
        % Select one transmit antenna BEFORE RIS optimization:
        % ant_sel = 1: maximum Bob CS 
        % ant_sel = 2: minimum Eve CS 
        % ant_sel = 3: maximum Bob/Eve SNR ratio, reference Theta = I
        % ant_sel = 4: maximum Bob-Eve SNR difference, reference Theta = I
        % ant_sel = 5: maximum secrecy ratio (1+gammaB)/(1+gammaE), Theta = I
        ant_sel = 1;
        rho =0.9;
        
        
        %% ===================== Algorithms =====================
        % Storage indices 1,2,3 represent the original algorithm IDs 1,2,4.
        idxMRT = 1;
        idxNull = 2;
        idxBobEve = 3;
        algIDs = [1 2 4];
        algNames = {'MRT','MRT + Nulling','Bob-Eve (preselected antenna)'};
        nAlg = numel(algNames);
        
        %% ===================== Nulling Parameters =====================
        num_iter = 100;
        epsilon = 1e-6;
        varepsilon = 1e-6;
        
        %% ===================== RIS Types =====================
        % 0: FC-RIS, 1: SC-RIS, 2/4: GC-RIS group size.
        NG_list = [0 1 2 4];
        RIS_names = {'FC','SC','GC-2','GC-4'};
        nRIS = numel(NG_list);
        nN = numel(locE_list);
        
        for iRIS = 1:nRIS
            group = NG_list(iRIS);
            if group > 1
                assert(mod(NI,group)==0, ...
                    'NI must be divisible by the GC-RIS group size.');
            end
        end
        
        %% ===================== Allocate Results =====================
        SNR_B = nan(nMonte,nN,nAlg,nRIS);
        SNR_E = nan(nMonte,nN,nAlg,nRIS);
        OBJ = nan(nMonte,nN,nAlg,nRIS);
        SR = nan(nMonte,nN,nAlg,nRIS);
        ANT_selected = nan(nMonte,nN,nAlg,nRIS);
        CS_B = nan(nMonte,nN,Ma);
        CS_E = nan(nMonte,nN,Ma);
        ANT_CS = nan(nMonte,nN);  
        SNR_B_ref = nan(nMonte,nN,Ma);
        SNR_E_ref = nan(nMonte,nN,Ma);
        RATIO_ref = nan(nMonte,nN,Ma);
        DIFF_ref = nan(nMonte,nN,Ma);
        
        %% ===================== Monte Carlo Simulation =====================
        for iMonte = 1:nMonte
            if mod(iMonte,100) == 0
                fprintf('Monte Carlo trial: %d / %d\n',iMonte,nMonte);
            end
        
            % Generate channels for each Eve location, shared across all algorithms.
            for iN = 1:nN
                locE = [locE_list(iN);0];
                [hBI_max,hIA_max,hEI_max,GBI,GIA,GEI] = generate_rician_channels( ...
                    Mb,Ma,Me,P,NI_max,K,rho,LS,locB,locA,locI,locE);
                hBI_N = hBI_max(:,1:NI);       % 1 x NI
                HIA_full = hIA_max(1:NI,:);    % NI x Ma
                hEI_N = hEI_max(:,1:NI);      % Me x NI
        
                %% ----- Preselect antenna using the chosen criterion -----
           
                alice_ants = sum(abs(HIA_full).^2,1);
                csBoundB = alice_ants;
                csBoundE = alice_ants;
                % csBoundB = norm(hBI_N,2)^2 * alice_ants;
                % 
                % csBoundE = norm(hEI_N,'fro')^2 * alice_ants;
                ratio_ref=csBoundB/csBoundE;
                switch ant_sel
                    case 1
                        [~,iAnt] = max(csBoundB);
                    case 2
                        [~,iAnt] = min(csBoundE);
                    case 3
                        [~,iAnt] = max(ratio_ref);
                    case 4
                        [~,iAnt] = max(diff_ref);
                    case 5
                        [~,iAnt] = max(secrecy_ratio_ref);
                    otherwise
                        error('ant_sel must be an integer from 1 to 5.');
                end
                ANT_CS(iMonte,iN) = iAnt;
                hIA_N = HIA_full(:,iAnt);
        
                
                    hBI_norm = hBI_N/norm(hBI_N);
                           hEI_norm = hEI_N/norm(hEI_N,'fro');
                      
                    hIA_norm = hIA_N/norm(hIA_N);
                
        
                %% ----- Loop over RIS architectures -----
                for iRIS = 1:nRIS
                    NG = NG_list(iRIS);
                    if NG == 0
                        NG_eff = NI;             % FC-RIS
                    else
                        NG_eff = NG;             % SC/GC-RIS
                    end
                    scale = 1;
        
                    %% =========================================================
                    % Algorithm 1: MRT for the PRESELECTED antenna
                    %% =========================================================
                    
                        [Theta_MRT,~] = func_MRT_GC( ...
                            hIA_norm,hBI_norm,NG_eff);
                    
        
                    gB_MRT = hBI_N*Theta_MRT*hIA_N;
                    gE_MRT = hEI_N*Theta_MRT*hIA_N;
                    gammaB_MRT = pt/N_B*abs(gB_MRT)^2/scale;
                    gammaE_MRT = pt/N_E*noiseRatio*sum(abs(gE_MRT).^2)/scale;
                    SR_MRT = max(log2(1+gammaB_MRT)-log2(1+gammaE_MRT),0);
        
                    ANT_selected(iMonte,iN,idxMRT,iRIS) = iAnt;
                    SNR_B(iMonte,iN,idxMRT,iRIS) = gammaB_MRT;
                    SNR_E(iMonte,iN,idxMRT,iRIS) = gammaE_MRT;
                    OBJ(iMonte,iN,idxMRT,iRIS) = gammaB_MRT-gammaE_MRT;
                    SR(iMonte,iN,idxMRT,iRIS) = SR_MRT;
        
                    %% =========================================================
                    % Algorithm 2: MRT-initialized nulling, SAME selected antenna
                    %% =========================================================
                    
                        init = 0;
                        [Theta_null,~] = func_Nulling_GC( ...
                            hIA_norm,hEI_norm,NG_eff,Theta_MRT, ...
                            num_iter,epsilon,varepsilon,init);
                    
        
                    gB_null = hBI_N*Theta_null*hIA_N;
                    gE_null = hEI_N*Theta_null*hIA_N;
                    gammaB_null = pt/N_B*abs(gB_null)^2/scale;
                    gammaE_null = pt/N_E*noiseRatio*sum(abs(gE_null).^2)/scale;
                    SR_null = max(log2(1+gammaB_null)-log2(1+gammaE_null),0);
        
                    ANT_selected(iMonte,iN,idxNull,iRIS) = iAnt;
                    SNR_B(iMonte,iN,idxNull,iRIS) = gammaB_null;
                    SNR_E(iMonte,iN,idxNull,iRIS) = gammaE_null;
                    OBJ(iMonte,iN,idxNull,iRIS) = gammaB_null-gammaE_null;
                    SR(iMonte,iN,idxNull,iRIS) = SR_null;
        
                    %% =========================================================
                    % Algorithm 4: Bob-Eve RIS design for the SAME preselected
                    % antenna used by Algorithms 1 and 2.
                    % No loop over transmit antennas, no secrecy-based TAS.
                    %% =========================================================
                    Q = hBI_N'*hBI_N - noiseRatio*(hEI_N'*hEI_N);
                    [V,D] = eig(Q);
                    [~,idxEig] = max(real(diag(D)));
                    u_max = V(:,idxEig);
                    u_max = u_max/norm(u_max);
        
                   
                        % Design using ONLY the preselected Alice antenna.
                        B = u_max*hIA_norm';
                        if NG_eff == 1
                            Theta_BobEve = diag(exp(1i*angle(diag(B))));
                        else
                            Theta_BobEve = zeros(NI,NI);
                            nGroups = NI/NG_eff;
                            for g = 1:nGroups
                                idx_g = (g-1)*NG_eff+(1:NG_eff);
                                B_g = B(idx_g,idx_g);
                                if norm(B_g,'fro') == 0
                                    Theta_g = eye(NG_eff);
                                else
                                    Theta_g = alg_symuni(B_g);
                                end
                                Theta_BobEve(idx_g,idx_g) = Theta_g;
                            end
                        end
                    
        
                    % Evaluate the exact same selected antenna used by Alg. 1/2.
                    gB_BobEve = hBI_N*Theta_BobEve*hIA_N;
                    gE_BobEve = hEI_N*Theta_BobEve*hIA_N;
                    gammaB_BobEve = pt/N_B*abs(gB_BobEve)^2/scale;
                    gammaE_BobEve = pt/N_E*noiseRatio*sum(abs(gE_BobEve).^2)/scale;
                    SR_BobEve = max(log2(1+gammaB_BobEve) ...
                        - log2(1+gammaE_BobEve),0);
        
                    ANT_selected(iMonte,iN,idxBobEve,iRIS) = iAnt;
                    SNR_B(iMonte,iN,idxBobEve,iRIS) = gammaB_BobEve;
                    SNR_E(iMonte,iN,idxBobEve,iRIS) = gammaE_BobEve;
                    OBJ(iMonte,iN,idxBobEve,iRIS) = ...
                        gammaB_BobEve-gammaE_BobEve;
                    SR(iMonte,iN,idxBobEve,iRIS) = SR_BobEve;
                end
            end
        end
    
        
        %% ===================== Taylor Approximations =====================
        SR_first = max((SNR_B-SNR_E)/log(2),0);
        SR_second = max(((SNR_B-SNR_E) ...
            - 0.5*(SNR_B.^2-SNR_E.^2))/log(2),0);
        
        %% ===================== Average Results =====================
        PR_B_avg = reshape(mean(SNR_B,1,'omitnan'),nN,nAlg,nRIS);
        PR_E_avg = reshape(mean(SNR_E,1,'omitnan'),nN,nAlg,nRIS);
        OBJ_avg = reshape(mean(OBJ,1,'omitnan'),nN,nAlg,nRIS);
        SR_avg = reshape(mean(SR,1,'omitnan'),nN,nAlg,nRIS);
        SR_first_avg = reshape(mean(SR_first,1,'omitnan'),nN,nAlg,nRIS);
        SR_second_avg = reshape(mean(SR_second,1,'omitnan'),nN,nAlg,nRIS);
        
        fprintf('Elapsed time: %.2f seconds\n',toc);
        
        %% ===================== Common Figure Settings =====================
        LineW = 1.25;
        MarkS = 7;
        MarkInd = 1:2:nN;
        C = colororder;
        RIS_color = C(1:nRIS,:);
        RIS_line = {'-','--','-.',':'};
        markers = {'x','d','*'};
        plotAlgs = 1:nAlg;
        
        %% ===================== Figure 1: Bob and Eve SNR =====================
        figure('Color','w','Units','inches','Position',[0 0 10 5]);
        t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
        ax1 = nexttile(t,1); hold(ax1,'on');
        for iAlg = plotAlgs
            for iRIS = 1:nRIS
                semilogy(ax1,locE_list,PR_B_avg(:,iAlg,iRIS), ...
                    [RIS_line{iRIS} markers{iAlg}], ...
                    'Color',RIS_color(iRIS,:),'LineWidth',LineW, ...
                    'MarkerSize',MarkS,'MarkerIndices',MarkInd);
            end
        end
        addLocalLegend(ax1,RIS_color,RIS_line,RIS_names,markers, ...
            algNames,algIDs,LineW,MarkS,plotAlgs);
        grid(ax1,'on'); box(ax1,'on');
        xlabel(ax1,{'Eve x-location (m)','(a) Bob'}, ...
            'Interpreter','latex','FontSize',12);
        ylabel(ax1,'Average received SNR','Interpreter','latex','FontSize',12);
        formatAxes(ax1,locE_list);
        
        ax2 = nexttile(t,2); hold(ax2,'on');
        for iAlg = plotAlgs
            for iRIS = 1:nRIS
                semilogy(ax2,locE_list,PR_E_avg(:,iAlg,iRIS), ...
                    [RIS_line{iRIS} markers{iAlg}], ...
                    'Color',RIS_color(iRIS,:),'LineWidth',LineW, ...
                    'MarkerSize',MarkS,'MarkerIndices',MarkInd);
            end
        end
        addLocalLegend(ax2,RIS_color,RIS_line,RIS_names,markers, ...
            algNames,algIDs,LineW,MarkS,plotAlgs);
        grid(ax2,'on'); box(ax2,'on');
        xlabel(ax2,{'Eve x-location (m)','(b) Eve'}, ...
            'Interpreter','latex','FontSize',12);
        ylabel(ax2,'Average received SNR','Interpreter','latex','FontSize',12);
        formatAxes(ax2,locE_list);
        
        %% ===================== Figure 2: Difference and Secrecy Rate =====================
        figure('Color','w','Units','inches','Position',[0 0 10 5]);
        t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
        ax1 = nexttile(t,1); hold(ax1,'on');
        for iAlg = plotAlgs
            for iRIS = 1:nRIS
                plot(ax1,locE_list,OBJ_avg(:,iAlg,iRIS), ...
                    [RIS_line{iRIS} markers{iAlg}], ...
                    'Color',RIS_color(iRIS,:),'LineWidth',LineW, ...
                    'MarkerSize',MarkS,'MarkerIndices',MarkInd);
            end
        end
        addLocalLegend(ax1,RIS_color,RIS_line,RIS_names,markers, ...
            algNames,algIDs,LineW,MarkS,plotAlgs);
        grid(ax1,'on'); box(ax1,'on');
        xlabel(ax1,{'Eve x-location (m)','(a) Difference'}, ...
            'Interpreter','latex','FontSize',12);
        ylabel(ax1,'Average SNR difference, $\gamma_B-\gamma_E$', ...
            'Interpreter','latex','FontSize',12);
        formatAxes(ax1,locE_list);
        
        ax2 = nexttile(t,2); hold(ax2,'on');
        for iAlg = plotAlgs
            for iRIS = 1:nRIS
                plot(ax2,locE_list,SR_avg(:,iAlg,iRIS), ...
                    [RIS_line{iRIS} markers{iAlg}], ...
                    'Color',RIS_color(iRIS,:),'LineWidth',LineW, ...
                    'MarkerSize',MarkS,'MarkerIndices',MarkInd);
            end
        end
        addLocalLegend(ax2,RIS_color,RIS_line,RIS_names,markers, ...
            algNames,algIDs,LineW,MarkS,plotAlgs);
        grid(ax2,'on'); box(ax2,'on');
        xlabel(ax2,{'Eve x-location (m)','(b) Secrecy rate'}, ...
            'Interpreter','latex','FontSize',12);
        ylabel(ax2,'Secrecy rate (bits/s/Hz)','Interpreter','latex','FontSize',12);
        formatAxes(ax2,locE_list);
        
        %% ===================== Figure 3: Taylor Accuracy =====================
        figure('Color','w','Units','inches','Position',[0 0 15 5]);
        t3 = tiledlayout(1,nAlg,'TileSpacing','compact','Padding','compact');
        for iAlg = plotAlgs
            ax = nexttile(t3,iAlg); hold(ax,'on');
            for iRIS = 1:nRIS
                plot(ax,locE_list,SR_avg(:,iAlg,iRIS),'-', ...
                    'Color',RIS_color(iRIS,:),'LineWidth',LineW);
                plot(ax,locE_list,SR_first_avg(:,iAlg,iRIS),'--', ...
                    'Color',RIS_color(iRIS,:),'LineWidth',LineW);
                plot(ax,locE_list,SR_second_avg(:,iAlg,iRIS),':', ...
                    'Color',RIS_color(iRIS,:),'LineWidth',1.7);
            end
            grid(ax,'on'); box(ax,'on');
            xlabel(ax,'Eve x-location (m)', ...
                'Interpreter','latex','FontSize',11);
            title(ax,sprintf('Alg. %d: %s',algIDs(iAlg),algNames{iAlg}), ...
                'Interpreter','none');
            formatAxes(ax,locE_list);
            if iAlg == idxMRT
                legendText = cell(1,3*nRIS);
                for iRIS = 1:nRIS
                    j = 3*(iRIS-1);
                    legendText{j+1} = [RIS_names{iRIS} ' exact'];
                    legendText{j+2} = [RIS_names{iRIS} ' first order'];
                    legendText{j+3} = [RIS_names{iRIS} ' first + second'];
                end
                legend(ax,legendText,'Location','northwest', ...
                    'Interpreter','none','FontSize',7,'NumColumns',2);
                ylabel(ax,'Average secrecy rate (bits/s/Hz)', ...
                    'Interpreter','latex','FontSize',11);
            end
        end
        
        %% ===================== Local Functions =====================
        function addLocalLegend(ax,RIS_color,RIS_line,RIS_names, ...
            markers,algNames,algIDs,LineW,MarkS,plotAlgs)
            h = gobjects(1,numel(RIS_names)+numel(plotAlgs));
            legendText = cell(1,numel(h));
            for iRIS = 1:numel(RIS_names)
                h(iRIS) = plot(ax,nan,nan,RIS_line{iRIS}, ...
                    'Color',RIS_color(iRIS,:),'LineWidth',LineW);
                legendText{iRIS} = [RIS_names{iRIS} '-RIS'];
            end
            for iPlot = 1:numel(plotAlgs)
                idx = plotAlgs(iPlot);
                k = numel(RIS_names)+iPlot;
                h(k) = plot(ax,nan,nan,['k' markers{idx}], ...
                    'LineStyle','none','MarkerSize',MarkS);
                legendText{k} = ['Alg. ' num2str(algIDs(idx)) ': ' algNames{idx}];
            end
            legend(ax,h,legendText,'Location','northwest', ...
                'Interpreter','none','FontSize',8,'NumColumns',2);
        end
        
        function formatAxes(ax,locE_list)
            ax.XLim = [min(locE_list) max(locE_list)];
            ax.XTick = locE_list;
            ax.FontName = 'Times New Roman';
            ax.FontSize = 12;
            ax.LineWidth = 0.75;
        end

