function [distance] = irfdistance(parametersestimation,irfdata,weights)

phikestimate            = 10*exp(parametersestimation(1));
rhovolaestimate 		= exp(parametersestimation(2));
volvolaestimate			= exp(parametersestimation(3))/100;
volassestimate			= exp(parametersestimation(4))/100;
volzssestimate			= exp(parametersestimation(5))/100;
rhoaestimate            = exp(parametersestimation(6));
rhozestimate            = exp(parametersestimation(7));

bbeffectivedemandmatchirf;

stdydata    = 1.10678;
stdcdata    = 0.710136;
stdinvdata  = 3.79229;
stdndata    = 1.41389;

unconditionalmomentsdata = [stdydata;stdcdata;stdinvdata;stdndata];

unconditionalmomentsmodel = [stdymodel;stdcmodel;stdinvmodel;stdnmodel];

unconditionalmomentsweightingmatrix = diag(unconditionalmomentsdata.^(-2));

weightonunconditionalmoments = 400;

if checkresult ~= 0
    irfmodel = [logvxo_evola;y_evola;c_evola;inv_evola;n_evola;pl_evola;r_evola];
    distance = (irfdata - irfmodel)'*weights*(irfdata - irfmodel) + ...
        weightonunconditionalmoments*(unconditionalmomentsdata - unconditionalmomentsmodel)'*unconditionalmomentsweightingmatrix*(unconditionalmomentsdata - unconditionalmomentsmodel);
else
    distance = 1e8;
end

disp([distance phikestimate rhovolaestimate volvolaestimate volassestimate volzssestimate rhoaestimate rhozestimate])


disp([unconditionalmomentsmodel]')
