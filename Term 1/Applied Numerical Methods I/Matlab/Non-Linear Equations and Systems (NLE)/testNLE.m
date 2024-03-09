clear, clc
% testing Non-linear methods for Equation Systems
x0 = [2;-1];
tol = 1e-20;
maxiter = 40;
[sol_Newton,iter_Newton,ACOC_Newton,incre1_Newton,incre2_Newton] = Newton_Sist('testFunctionNLSE',vpa(x0),tol, maxiter);
[sol_Trapecios,iter_Trapecios,ACOC_Trapecios,incre1_Trapecios,incre2_Trapecios] = Trapecios_Sist('testFunctionNLSE',vpa(x0),tol, maxiter);
[sol_PuntoMedio,iter_PuntoMedio,ACOC_PuntoMedio,incre1_PuntoMedio,incre2_PuntoMedio] = PuntoMedio_Sist('testFunctionNLSE',vpa(x0),tol, maxiter);
[sol_Simpson,iter_Simpson,ACOC_Simpson,incre1_Simpson,incre2_Simpson] = Simpson_Sist('testFunctionNLSE',vpa(x0),tol, maxiter);
[sol_RN,iter_RN,ACOC_RN,incre1_RN,incre2_RN] = RN_Sist('testFunctionNLSE',vpa(x0),tol, maxiter);
[sol_Traub,iter_Traub,ACOC_Traub,incre1_Traub,incre2_Traub] = Traub_Sist('testFunctionNLSE',vpa(x0),tol, maxiter);
[sol_GR,iter_GR,ACOC_GR,incre1_GR,incre2_GR] = GR_Sist('testFunctionNLSE',vpa(x0),tol, maxiter);
[sol_NA,iter_NA,ACOC_NA,incre1_NA,incre2_NA] = NA_Sist('testFunctionNLSE',vpa(x0),tol, maxiter);

metodos = {'Newton', 'Trapecios', 'PuntoMedio', 'Simpson', 'RN', 'Traub', 'GR', 'NA'};

resultados = [iter_Newton, incre1_Newton, incre2_Newton, ACOC_Newton(end);
              iter_Trapecios, incre1_Trapecios, incre2_Trapecios, ACOC_Trapecios(end);
              iter_PuntoMedio, incre1_PuntoMedio, incre2_PuntoMedio, ACOC_PuntoMedio(end);
              iter_Simpson, incre1_Simpson, incre2_Simpson, ACOC_Simpson(end);
              iter_RN, incre1_RN, incre2_RN, ACOC_RN(end);
              iter_Traub, incre1_Traub, incre2_Traub, ACOC_Traub(end);
              iter_GR, incre1_GR, incre2_GR, ACOC_GR(end);
              iter_NA, incre1_NA, incre2_NA, ACOC_NA(end)];

tabla_resultados = table(metodos', resultados(:,1), resultados(:,2), ...
    resultados(:,3), resultados(:,4), 'VariableNames', ...
    {'Metodo', 'Iteraciones', 'Incre1', 'Incre2', 'ACOC'});

disp(tabla_resultados);

