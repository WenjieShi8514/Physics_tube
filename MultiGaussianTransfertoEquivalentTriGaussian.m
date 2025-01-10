function [EquivalentTriGaussian_parameter,area_sumsum] = MultiGaussianTransfertoEquivalentTriGaussian(x,y,rho0_best,rho_inf_array,E_array,w_array,break_point)

%%% Convert multi-Gaussian fitting results to equivalent three-Gaussian distribution
%%% Steps:
%%% 1. Calculate the areas of the fitted multi-Gaussian distribution in the basal domain 
%%%    anterior/posterior-lateral domain regions.
%%% 2. Fix the area of the basal domain region and fit the equivalent mono-Gaussian 
%%%    distribution for the basal domain region to obtain its characteristic parameters.
%%% 3. Recalculate the areas of the anterior/posterior-lateral domain regions and fine-tune
%%%    the parameters of the three-Gaussian distribution to ensure area conservation across 
%%%    the basal domain and the anterior/posterior-lateral domain regions.


%%% step 1
basal_area=integral(@y_MultiGaussian,break_point(1),break_point(2));
anterior_lateral_area=integral(@y_MultiGaussian,min(x),break_point(1));
posterior_lateral_area=integral(@y_MultiGaussian,break_point(2),max(x));
area_sum1=[basal_area anterior_lateral_area posterior_lateral_area];
    function yy = y_MultiGaussian(x)
        yy = rho0_best;
        for i=1:length(rho_inf_array)
            yy=yy+rho_inf_array(i)*exp(-1/2*((x-E_array(i))/w_array(i)).^2);
        end
    end


%%% step 2
ft=fittype(['(basal_area-integral(@(x,break_point1,break_point2) rho0+rho_inf_anterior*exp(-1/2*((x-E_anterior)/w_anterior).^2)+rho_inf_posterior*exp(-1/2*((x-E_posterior)/w_posterior).^2),break_point1,break_point2))/' ...
    'integral(@(x,break_point1,break_point2) exp(-1/2*((x-E)/w).^2),break_point1,break_point2)*exp(-1/2*((x-E)/w).^2)+rho_inf_anterior*exp(-1/2*((x-E_anterior)/w_anterior).^2)+rho_inf_posterior*exp(-1/2*((x-E_posterior)/w_posterior).^2)'],'independent','x','dependent','y');
opts=fitoptions('Method','NonlinearLeastSquares');
opts.Algorithm='Levenberg-Marquardt';
opts.Display='Off';
opts.MaxFunEvals=600;
opts.MaxIter=400;
opts.Lower=[-0.1 E_array(end-1) E_array(end) basal_area break_point(1) break_point(2)...
    rho0_best rho_inf_array(end-1) rho_inf_array(end) 0.001 w_array(end-1) w_array(end)];
opts.Upper=[0.1 E_array(end-1) E_array(end) basal_area break_point(1) break_point(2)...
    rho0_best rho_inf_array(end-1) rho_inf_array(end) 1 w_array(end-1) w_array(end)];
[xData,yData]=prepareCurveData(x,y);
% Fit model to data.
wucha_best=1000000000;
for i_fit=1:10  % 循环10次，寻找最优拟合解
    wucha=0;
    [fitresult,~]=fit(xData,yData,ft,opts);
    % Fit parameter result to variables
    E=fitresult.E;
    w=fitresult.w;
    rho_inf=(basal_area-integral(@(x) rho0_best+rho_inf_array(end-1)*exp(-1/2*((x-E_array(end-1))/w_array(end-1)).^2)+rho_inf_array(end)*exp(-1/2*((x-E_array(end))/w_array(end)).^2),break_point(1),break_point(2)))/...
        integral(@(x) exp(-1/2*((x-E)/w).^2),break_point(1),break_point(2));
    for ii=1:length(x)
        wucha=wucha+(rho0_best+rho_inf*exp(-1/2*((x(ii)-E)/w).^2)+rho_inf_array(end-1)*exp(-1/2*((x(ii)-E_array(end-1))/w_array(end-1)).^2)+rho_inf_array(end)*exp(-1/2*((x(ii)-E_array(end))/w_array(end)).^2)-y(ii))^2;
    end

    if wucha<wucha_best
        wucha_best=wucha;
        rho_inf_best=rho_inf;
        w_best=w;
        E_best=E;
    end
end
EquivalentTriGaussian_parameter=[rho0_best rho_inf_best rho_inf_array(end-1) rho_inf_array(end) w_best w_array(end-1) w_array(end) E_best E_array(end-1) E_array(end)];
basal_area1=integral(@y_TriGaussian_temp,break_point(1),break_point(2));
anterior_lateral_area1=integral(@y_TriGaussian_temp,min(x),break_point(1));
posterior_lateral_area1=integral(@y_TriGaussian_temp,break_point(2),max(x));
area_sum2=[basal_area1 anterior_lateral_area1 posterior_lateral_area1];
    function yy = y_TriGaussian_temp(x)
        yy = rho0_best;
        for i=1:3
            yy=yy+EquivalentTriGaussian_parameter(i+1)*exp(-1/2*((x-EquivalentTriGaussian_parameter(i+7))/EquivalentTriGaussian_parameter(i+4)).^2);
        end
    end


%%% step 3
ft=fittype(['rho0+rho_inf*exp(-1/2*((x-E)/w).^2)+' ...
    '(anterior_lateral_area-integral(@(x,break_point1) rho0+rho_inf*exp(-1/2*((x-E)/w).^2),-0.5,break_point1))/' ...
    'integral(@(x,break_point1) exp(-1/2*((x-E_anterior)/w_anterior).^2),-0.5,break_point1)*exp(-1/2*((x-E_anterior)/w_anterior).^2)+' ...
    '(posterior_lateral_area-integral(@(x,break_point2) rho0+rho_inf*exp(-1/2*((x-E)/w).^2),break_point2,0.5))/' ...
    'integral(@(x,break_point2) exp(-1/2*((x-E_posterior)/w_posterior).^2),break_point2,0.5)*exp(-1/2*((x-E_posterior)/w_posterior).^2)'],'independent','x','dependent','y');
opts=fitoptions('Method','NonlinearLeastSquares');
opts.Algorithm='Levenberg-Marquardt';
opts.Display='Off';
opts.MaxFunEvals=600;
opts.MaxIter=400;
opts.Lower=[E_best E_array(end-1) E_array(end) anterior_lateral_area break_point(1) break_point(2) posterior_lateral_area...
    rho0_best rho_inf_best w_best w_array(end-1)*0.9 w_array(end)*0.9];
opts.Upper=[E_best E_array(end-1) E_array(end) anterior_lateral_area break_point(1) break_point(2) posterior_lateral_area...
    rho0_best rho_inf_best w_best max(w_array(end-1)*1.1,w_array(end-1)+0.1) max(w_array(end)*1.1,w_array(end)+0.1)];
[xData,yData]=prepareCurveData(x,y);
% Fit model to data.
wucha_best=1000000000;
for i_fit=1:10  % 循环10次，寻找最优拟合解
    wucha=0;
    [fitresult,~]=fit(xData,yData,ft,opts);
    % Fit parameter result to variables
    E_anterior=fitresult.E_anterior;
    w_anterior=fitresult.w_anterior;
    E_posterior=fitresult.E_posterior;
    w_posterior=fitresult.w_posterior;
    rho_inf_anterior=(anterior_lateral_area-integral(@(x) rho0_best+rho_inf_best*exp(-1/2*((x-E_best)/w_best).^2),-0.5,break_point(1)))/...
        integral(@(x) exp(-1/2*((x-E_anterior)/w_anterior).^2),-0.5,break_point(1));
    rho_inf_posterior=(posterior_lateral_area-integral(@(x) rho0_best+rho_inf_best*exp(-1/2*((x-E_best)/w_best).^2),break_point(2),0.5))/...
        integral(@(x) exp(-1/2*((x-E_posterior)/w_posterior).^2),break_point(2),0.5);
    for ii=1:length(x)
        wucha=wucha+(rho0_best+rho_inf_best*exp(-1/2*((x(ii)-E_best)/w_best).^2)+rho_inf_anterior*exp(-1/2*((x(ii)-E_array(end-1))/w_anterior).^2)+rho_inf_posterior*exp(-1/2*((x(ii)-E_array(end))/w_posterior).^2)-y(ii))^2;
    end

    if wucha<wucha_best
        wucha_best=wucha;
        rho_inf_anterior_best=rho_inf_anterior;
        rho_inf_posterior_best=rho_inf_posterior;
        w_anterior_best=w_anterior;
        w_posterior_best=w_posterior;
    end
end
EquivalentTriGaussian_parameter=[rho0_best rho_inf_best rho_inf_anterior_best rho_inf_posterior_best w_best w_anterior_best w_posterior_best E_best E_array(end-1) E_array(end)];
basal_area2=integral(@y_TriGaussian,break_point(1),break_point(2));
anterior_lateral_area2=integral(@y_TriGaussian,min(x),break_point(1));
posterior_lateral_area2=integral(@y_TriGaussian,break_point(2),max(x));
area_sum3=[basal_area2 anterior_lateral_area2 posterior_lateral_area2];
    function yy = y_TriGaussian(x)
        yy = rho0_best;
        for i=1:3
            yy=yy+EquivalentTriGaussian_parameter(i+1)*exp(-1/2*((x-EquivalentTriGaussian_parameter(i+7))/EquivalentTriGaussian_parameter(i+4)).^2);
        end
    end
rho_inf_best=rho_inf_best*basal_area1/basal_area2;
EquivalentTriGaussian_parameter=[rho0_best rho_inf_best w_best E_best rho_inf_anterior_best w_anterior_best E_array(end-1) rho_inf_posterior_best w_posterior_best E_array(end)];

basal_area3=integral(@y_TriGaussian2,break_point(1),break_point(2));
anterior_lateral_area3=integral(@y_TriGaussian2,min(x),break_point(1));
posterior_lateral_area3=integral(@y_TriGaussian2,break_point(2),max(x));
area_sum4=[basal_area3 anterior_lateral_area3 posterior_lateral_area3];
    function yy = y_TriGaussian2(x)
        yy = rho0_best;
        for i=1:3
            yy=yy+EquivalentTriGaussian_parameter(3*i-1)*exp(-1/2*((x-EquivalentTriGaussian_parameter(3*i+1))/EquivalentTriGaussian_parameter(3*i)).^2);
        end
    end
area_sumsum=[sum(area_sum1) sum(area_sum2) sum(area_sum3) sum(area_sum4)];
end
