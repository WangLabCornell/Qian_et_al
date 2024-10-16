function z = Nucleosome_array_hat_curve(r0,x) 
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
