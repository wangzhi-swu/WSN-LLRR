function  [SigmaX, svp] =  WSNM(Opt, SigmaY, C, p )

    p1    =   p;
% 	Temp = SigmaY;
    Temp = Opt.siagS1;
	s = SigmaY;
    s1 = zeros(size(s));
    for i=1:3
%         W_Vec = C.*ones(size(SigmaY));
        W_Vec    =   (C*sqrt(Opt.d*Opt.n))./( Temp + eps );
        [s1, svp] =   solve_Lp_w(s, W_Vec, p1);
        Temp     =   s1;
    end
    SigmaX = s1;
return;
