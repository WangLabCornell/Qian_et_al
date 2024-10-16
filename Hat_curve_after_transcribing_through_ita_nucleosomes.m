function new_C0 = Hat_curve_after_transcribing_through_ita_nucleosomes (ita, dn, dH, C0_i)
Li = 12667;
lempty = 3260.2545/64;
lnuc = 530.70309/64;
mu = (dH-(ita*lnuc))/lempty;
dx_naked = (mu+ita)*197;
Lf = Li-dx_naked;
[~,p_i] = getNucQuality(C0_i);
dp = -ita;
C0_i(:,1) = C0_i(:,1)+dn;
new_C0 = Predict_Hatcurve_NonTC (C0_i, Li, Lf, p_i, dp);


%% Nested functions

    function [nucQuality,n_nuc] = getNucQuality(r0)
        %obtains goodness of a nucleosome array based on Tung's method described in
        %190119_Tung_Selection criteria for single_double chromatin array data on
        %AOT_MT.pptx assumes nuc array is from 197-64ex plasmid
        %modified based on 190205_Seong ha_Hat curve dimensions in different buffer conditions.pptx
        %load('x_f_WLC.mat');
        
        r = reshape(r0,4,2);
        %load('nVLatF.mat');
        %slope = interp1(nVLatF(:,1),nVLatF(:,2),F); % N Vs Length relationship obtained on 190520
        slope = -0.0414;
        %intercept = interp1(nVLatF(:,1),nVLatF(:,3),F);
        intercept = 3.2306;
        [Height, ~, ~, buckling_positive, ~,~] = HCpara(r);
        %dL = 0.200; % length measurement uncertainty on MT
        n_nuc =  (Height - intercept)/slope;
        n_nucRange = [n_nuc-5 n_nuc+5];
        widthRange = 22.7936+ 0.5973*n_nucRange;
        
        if buckling_positive > widthRange(1) && buckling_positive < widthRange(2) %abs(n_nuc  - n_nuc_width) < dn_nuc_width
            nucQuality = 'good';
        else
            nucQuality = 'bad';
        end
        
    end

    function [Height, hatcenter, buckling_negative, buckling_positive, slope_negative, slope_positive] = HCpara(r0)
        %obtains hat curve parameters from the result of a 5 piece fit
        % buckling positive is buckling point relative to the hat center
        
        
        r = reshape(r0,4,2);
        M_tl = [2*r(2,1) 1 0;  r(1,1)^2  r(1,1) 1; r(2,1)^2  r(2,1) 1];
        C_tl = [(r(3,2)-r(2,2))/(r(3,1)-r(2,1)) ;  r(1,2);  r(2,2)];
        para_tl= M_tl\C_tl;
        
        M_tr = [2*r(3,1) 1 0;  r(4,1)^2  r(4,1) 1;  r(3,1)^2  r(3,1) 1];
        C_tr = [(r(3,2)-r(2,2))/(r(3,1)-r(2,1)) ;  r(4,2);  r(3,2)];
        para_tr = M_tr\C_tr;
        alpha_n = 2 * para_tl(1) * r(1,1) + para_tl(2);
        alpha_p = 2 * para_tr(1) * r(4,1) + para_tr(2);
        slope_negative = alpha_n;
        slope_positive = alpha_p;
        alpha_c = (r(3,2)-r(2,2))/(r(3,1)-r(2,1));
        
        buckling_negative = (r(2,2) - r(1,2) + alpha_n * r(1,1) - alpha_c * r(2,1))/(alpha_n-alpha_c);
        buckling_positive = (r(4,2) - r(3,2) - alpha_p * r(4,1) + alpha_c * r(3,1))/(alpha_c-alpha_p);
        
        %% top left curve
        M_tl = [2*r(2,1) 1 0;
            r(1,1)^2  r(1,1) 1;
            r(2,1)^2  r(2,1) 1];
        C_tl = [(r(3,2)-r(2,2))/(r(3,1)-r(2,1)) ;
            r(1,2);
            r(2,2)];
        para_tl= M_tl\C_tl;
        Height = para_tl(3) - para_tl(2)^2/(4*para_tl(1));
        hatcenter =  -para_tl(2)/2/para_tl(1);
        buckling_positive = buckling_positive - hatcenter;
        buckling_negative = buckling_negative - hatcenter;
    end

    function new_C0 = Predict_Hatcurve_NonTC (C0_i, Li, Lf, p_i, dp)
        L_standard = 12667;
        C0_scaled = C0_i/Li*L_standard;
        dp_scaled = dp/Lf*L_standard+p_i*L_standard*(1/Lf-1/Li);
        new_C0_scaled = delta_hat_curve (C0_scaled, dp_scaled);
        new_C0 = new_C0_scaled*Lf/L_standard;
        
    end

    function new_C0 = delta_hat_curve (C0, dp)
        
        [~, ~, wL, wR, kL, kR] = HCpara(C0);
        
        kM = (C0(2,2)-C0(3,2))/(C0(2,1)-C0(3,1));
        Zero_turn_height = f_5piece(C0,0);
        
        k_height = -42.6492/1000;
        k_y2 = -38.4353/1000;
        k_kL = -(61.76892-17.56061)/65/1000;
        k_kM = -(-4.39694+4.17612)/65/1000;
        k_kR = -(-57.32513+17.63203)/65/1000;
        k_wL = 0.03915;
        k_wR = 0.54257;
        k_x2 = (1.23/4.97*25)/(13.22/6.4*70);
        k_x3 = (6.39/5.28*60)/(11.38/6.96*70);
        
        new_kL = kL+k_kL*dp;
        new_kR = kR+k_kR*dp;
        new_kM = kM+k_kM*dp;
        new_wL = wL+k_wL*dp;
        new_wR = wR+k_wR*dp;
        new_x2 = C0(2,1)+k_x2*dp;
        new_x3 = C0(3,1)+k_x3*dp;
        new_height = Zero_turn_height+k_height*dp;
        
        xmu = -wL;
        QL = new_kL*xmu;
        QM = new_kM*xmu;
        ALPHAL = 0.25*(new_kL^2-new_kM^2)/(QM-QL);
        %CL = QL+new_kL^2/(4*ALPHAL);
        %x1_tilde = new_kL/(2*ALPHAL)+xmu;
        x2_tilde = new_kM/(2*ALPHAL)+xmu;
        wJL = new_x2-x2_tilde;
        wJR = new_wR-new_wL+wJL;
        
        new_C0 = construct_yshifted_hatcurve_xpara (new_x2, new_x3, new_kL, new_kM, new_kR, wJL, wJR);
        new_C0(:,2) = new_C0(:,2)+new_height;
        
    end

    function r = fit5piece(turn, zBead)
        
        height = max(zBead);
        indexmax = find(zBead == max(zBead)); indexmax = indexmax(1);
        peakTurn = mean(turn(indexmax));
        F_positive = @(para, x) (height + para(1) * (x-peakTurn)) .* (x <= para(3)) + ((height + para(1) * (para(3) - peakTurn)) + para(2) * (x - para(3))) .* (x > para(3));
        para0 = [-0.004, -0.020 , 40+peakTurn];
        [para_positive,~,~,~,~, ~, ~] = lsqcurvefit(F_positive,para0,turn(turn >= peakTurn), zBead(turn >= peakTurn));
        turn_buckling_positive = para_positive(3) ;
        
        F_negative = @(para, x) (height + para(1) * (x - peakTurn)) .* (x >= para(3)) + ((height + para(1) * (para(3)-peakTurn)) + para(2) * (x - para(3))) .* (x < para(3));
        para0 = [0.001, 0.020 , -10+peakTurn];
        [para_negative,~,~,~,~, ~, ~] = lsqcurvefit(F_negative,para0, turn(turn < peakTurn), zBead(turn < peakTurn));
        turn_buckling_negative = para_negative(3);
        
        r0 = [ turn_buckling_negative F_negative(para_negative,turn_buckling_negative);
            peakTurn F_negative(para_negative,peakTurn);
            turn_buckling_positive-1 F_positive(para_positive,turn_buckling_positive-1);
            turn_buckling_positive+1 F_positive(para_positive,turn_buckling_positive+1)];
        
        
        r = lsqcurvefit(@f_5piece,r0,turn,zBead);
        
    end
    function z = f_5piece(r0,x)
        % r is the Catersian coordinates of 4 points that define the 5 piece
        % function
        r = reshape(r0,4,2);
        % top piece
        y3 = r(2,2) + (r(3,2)-r(2,2))/(r(3,1)-r(2,1)).*(x-r(2,1));
        %% the two curving parts
        %% top left curve
        M_tl = [2*r(2,1) 1 0;
            r(1,1)^2  r(1,1) 1;
            r(2,1)^2  r(2,1) 1];
        C_tl = [(r(3,2)-r(2,2))/(r(3,1)-r(2,1)) ;
            r(1,2);
            r(2,2)];
        para_tl= M_tl\C_tl;
        y2 = para_tl(1) * x.^2 + para_tl(2) * x + para_tl(3);
        
        %% top right curve
        M_tr = [2*r(3,1) 1 0;
            r(4,1)^2  r(4,1) 1;
            r(3,1)^2  r(3,1) 1];
        C_tr = [(r(3,2)-r(2,2))/(r(3,1)-r(2,1)) ;
            r(4,2);
            r(3,2)];
        para_tr = M_tr\C_tr;
        y4 = para_tr(1) * x.^2 + para_tr(2) * x + para_tr(3);
        
        %% calculating alpha_n and alpha_p
        % lower left piece
        alpha_n = 2 * para_tl(1) * r(1,1) + para_tl(2);
        y1 = r(1,2) + alpha_n * (x-r(1,1));
        
        % lower right piece
        alpha_p = 2 * para_tr(1) * r(4,1) + para_tr(2);
        y5 = r(4,2) + alpha_p * (x-r(4,1));
        
        %% whole curve
        z = y1 .* (x<r(1,1)) + y2 .* (x>=r(1,1) & x<r(2,1)) + y3 .* (x>=r(2,1) & x<r(3,1)) + y4 .* (x>=r(3,1) & x<r(4,1)) + y5.* (x>=r(4,1));
        
    end

    function C_output = construct_yshifted_hatcurve_xpara (new_x2, new_x3, new_kL, new_kM, new_kR, wJL, wJR)
        
        YL = new_kM*wJL;
        bL = YL-new_kL*wJL;
        qL = bL+new_kL*new_x2;
        qM1 = new_kM*new_x2;
        new_x1 = new_x2 + 2*(qM1-qL)/(new_kL-new_kM);
        alpha_L = 0.5*(new_kL-new_kM)/(new_x1-new_x2);
        
        YR = new_kM*wJR;
        bR = YR-new_kR*wJR;
        qR = bR+new_kR*new_x3;
        qM2 = new_kM*new_x3;
        new_x4 = new_x3 + 2*(qM2-qR)/(new_kR-new_kM);
        alpha_R = 0.5*(new_kR-new_kM)/(new_x4-new_x3);
        
        turns = [new_x1, new_x2, new_x3, new_x4, 0];
        
        SI = heaviside(new_x1-turns).*(new_kL*turns+bL);
        SII = heaviside(turns-new_x1).*heaviside(new_x2-turns).*( alpha_L*(turns-new_x2).^2+new_kM*(turns-new_x2)+qM1 );
        SIII = heaviside(turns-new_x2).*heaviside(new_x3-turns).*( new_kM*turns );
        SIV = heaviside(turns-new_x3).*heaviside(new_x4-turns).*( alpha_R*(turns-new_x3).^2+new_kM*(turns-new_x3)+qM2 );
        SV = heaviside(turns-new_x4).*(new_kR*turns+bR);
        S_total = SI+SII+SIII+SIV+SV;
        
        mock_z0 = S_total(5);
        mock_y1 = S_total(1);
        mock_y2 = S_total(2);
        mock_y3 = S_total(3);
        mock_y4 = S_total(4);
        new_y1 = mock_y1-mock_z0;
        new_y2 = mock_y2-mock_z0;
        new_y3 = mock_y3-mock_z0;
        new_y4 = mock_y4-mock_z0;
        
        C_output = [new_x1, new_y1; new_x2, new_y2; new_x3, new_y3; new_x4, new_y4];
        
    end

end