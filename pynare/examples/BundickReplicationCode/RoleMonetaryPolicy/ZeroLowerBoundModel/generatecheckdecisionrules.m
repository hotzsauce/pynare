% generatecheckdecisionrules.m
%

clear;

parameters;
format long g;
disp (' ')

generaterule = 1;

if ((exist('decisionrulen.dat') == 2) && (exist('decisionrulei.dat') == 2) && (exist('decisionrulepie.dat') == 2) && (exist('decisionruleu.dat') == 2) &&(exist('decisionruleexpvf1sigma.dat') == 2)) 
    load decisionrulen.dat
    load decisionrulei.dat    
    load decisionrulepie.dat    
    load decisionruleu.dat   
    load decisionruleexpvf1sigma.dat    
    
    if ((length(decisionrulen) == na*nk*nvola*nlagy) && (length(decisionrulei) == na*nk*nvola*nlagy) && (length(decisionrulepie) == na*nk*nvola*nlagy) && (length(decisionruleu) == na*nk*nvola*nlagy) && (length(decisionruleexpvf1sigma) == na*nk*nvola*nlagy))
       generaterule = 0;
       disp('Using Previously Saved Decision Rules')
       disp (' ')
    end
end


if generaterule == 1
    disp('Generating New Guesses Decision Rules')
    disp(' ')

    rule = zeros(na,nk,nvola,nlagy,5);
    
    nrulefile = fopen('decisionrulen.dat','w');
    irulefile = fopen('decisionrulei.dat','w');
    pierulefile = fopen('decisionrulepie.dat','w');
    urulefile = fopen('decisionruleu.dat','w');
    expvf1sigmarulefile = fopen('decisionruleexpvf1sigma.dat','w');

% POLICY AND TRANSITION FUNCTIONS
%                              n               pie             inv             u               expvf1sigma
% Constant                    0.325994        1.005000        0.268235        1.000000        1.000000
% k(-1)                      -0.002525        0.000135       -0.011732       -0.086183       -0.035866
% y(-1)                       0.155076        0.040096        0.581139        0.717684       -0.018570
% ea                         -0.133696        0.049025       -0.754723        0.221756       13.665664

    for i1 = 1:na
        for i2 = 1:nk
            for i3 = 1:nvola
	            for i4 = 1:nlagy                    
                    
				rule(i1,i2,i3,i4,1) = nss - 0.002525*(gridk(i2)-kss) - 0.133696*(grida(i1)-ass) + 0.155076*(gridlagy(i4)-yss);
				rule(i1,i2,i3,i4,2) = iss - 0.011732*(gridk(i2)-kss) - 0.754723*(grida(i1)-ass) + 0.581139*(gridlagy(i4)-yss);
				rule(i1,i2,i3,i4,3) = piess + 0.000135*(gridk(i2)-kss) + 0.049025*(grida(i1)-ass) + 0.040096*(gridlagy(i4)-yss);
                rule(i1,i2,i3,i4,4) = uss - 0.086183*(gridk(i2)-kss) + 0.221756*(grida(i1)-ass) + 0.717684*(gridlagy(i4)-yss);
                rule(i1,i2,i3,i4,5) = max(0.2,expvf1sigmass - 0.034866*(gridk(i2)-kss) + 13.665664*(grida(i1)-ass)) - 0.018570*(gridlagy(i4)-yss);
	
				fprintf(nrulefile,'%.16f\n',rule(i1,i2,i3,i4,1));
				fprintf(irulefile,'%.16f\n',rule(i1,i2,i3,i4,2));
				fprintf(pierulefile,'%.16f\n',rule(i1,i2,i3,i4,3));
                fprintf(urulefile,'%.16f\n',rule(i1,i2,i3,i4,4));
				fprintf(expvf1sigmarulefile,'%.16f\n',rule(i1,i2,i3,i4,5));

				end
            end
        end
    end

    fclose(nrulefile);
    fclose(irulefile);
    fclose(pierulefile);
    fclose(urulefile);
    fclose(expvf1sigmarulefile);
    
end

%---------------------------------------------------------------------
