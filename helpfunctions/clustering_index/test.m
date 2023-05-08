function a = test()
A = [1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 4 1 3 2];
B = [1 2 1 1 1 1 1 2 2 2 2 3 1 1 3 3 3 1 1 2 2];
NMI_1 = nmi_1(A, B);
NMI   = nmi(A, B);
R = rand_index(A,B);
aR = rand_index(A,B,'adjusted');

[AR, RI ,MI , HI ]=RandIndex(A,B);
c=0;
end










