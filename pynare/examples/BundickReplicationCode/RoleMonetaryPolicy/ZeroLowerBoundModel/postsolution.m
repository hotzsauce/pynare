
parameters


rule = zeros(na,nk,nvola,nlagy,5);

load decisionrulen.dat
load decisionrulei.dat
load decisionrulepie.dat
load decisionruleu.dat
load decisionruleexpvf1sigma.dat

vectorcounter = 1;

for i1 = 1:na
    for i2 = 1:nk
        for i3 = 1:nvola
            for i4 = 1:nlagy

                rule(i1,i2,i3,i4,1) = decisionrulen(vectorcounter,1);
                rule(i1,i2,i3,i4,2) = decisionrulei(vectorcounter,1);
                rule(i1,i2,i3,i4,3) = decisionrulepie(vectorcounter,1);
                rule(i1,i2,i3,i4,4) = decisionruleu(vectorcounter,1);
                rule(i1,i2,i3,i4,5) = decisionruleexpvf1sigma(vectorcounter,1);

                vectorcounter = vectorcounter + 1;

            end
        end
    end
end
