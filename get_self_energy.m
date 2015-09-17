function [SL,SR] = get_self_energy(Energy,M,HL0,HL1,HT)

NE = length(Energy);
SL = zeros(2*M,2*M,NE);
SR = zeros(2*M,2*M,NE);
j = 0;
for j=1:NE
    E = Energy(j);
    [is_div_L,is_div_R,sL,sR] = cal_self_energy(E,HL0,HL1,HT); % self-energy, retarded        
    SL(:,:,j) = sL;
    SR(:,:,j) = sR;
    clear sL sR     
    j = j + 1;   
end