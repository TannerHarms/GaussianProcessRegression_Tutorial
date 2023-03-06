clear all; close all; clc; 
set(groot, 'defaultAxesTickLabelInterpreter','latex')
% set up the parameters we'll use here
W          = linspace(0,1,100); 
choice_flag= 0; %set to 1 if you want an adaptive scheme based on acquisition function
der_flag   = 0; %set to 1 if you want to use derivative info
kernel_flag= 'RBF'; %'1/2' or 'RBF'  - gives matern 3/2, 1/2 or radial basis function
f          = @(x) .1*sin(5*pi*x)+abs(5*(x-.5)); %function we want to model 
df         = @(x) .1*5*pi*cos(5*pi*x)+5*sign(5*(x-.5));  %derivative
S1         = f(W);
%
clc; close all;
rng(1)
ii         = randperm(length(W));
%kernel functions: 
if strcmp(kernel_flag,'3/2')
%matern 3/2 kernel: 
    kernel_fun = @(x1,x2,s,l) s.^2.*(1 + sqrt(3)./l.*abs(x1-x2)).*exp(-sqrt(3)/l*abs(x1-x2)); %kernel function
    kernel_der = @(x1,x2,s,l) s.^2.*(-sqrt(3)./(l.^2).*(x1-x2)).*exp(-sqrt(3)/l*abs(x1-x2)); %derivative of kernel function
    kernel_hes = @(x1,x2,s,l) s.^2.*(3./(l.^2) - 3*sqrt(3)./(l.^3).*abs(x1-x2)).*exp(-sqrt(3)/l*abs(x1-x2)); %second derivative of kernel function
elseif strcmp(kernel_flag,'1/2')
%matern 1/2 kernel: 
    kernel_fun = @(x1,x2,s,l)  s.^2.*(exp(-abs(x1-x2)./l)); %kernel function
    kernel_der = @(x1,x2,s,l) -s.^2.*(exp(-abs(x1-x2)./l))./l.*sign(x1-x2); %derivative of kernel function
    kernel_hes = @(x1,x2,s,l) -s.^2.*(exp(-abs(x1-x2)./l))./l.^2; %second derivative of kernel function
elseif strcmp(kernel_flag,'RBF')
    %radial basis function:
    kernel_fun = @(x1,x2,s,l) s.^2.*exp(-(x1-x2).^2./(2*l.^2)); %kernel function
    kernel_der = @(x1,x2,s,l) (x2-x1)/l.^2.*kernel_fun(x1,x2,s,l); %derivative of kernel function
    kernel_hes = @(x1,x2,s,l) kernel_fun(x1,x2,s,l).*(1./l.^2 - (x2-x1).^2./l.^4); %second derivative of kernel function
end
%hyper parameters 
s = .1; l = .5; i = 10^-2; %in practice, a robust code optimizes the choice os s and l and i 
kappa = 1; %only for the acquisition function
x_data = 0; y_data = 0; d_data = 0; %initialize
figure(1)
colors = [0.8941    0.1020    0.1098;0.2157    0.4941    0.7216;0.3020    0.6863    0.2902;0.5961    0.3059    0.6392];%cbrewer('qual','Set1',4);
%this gets the original set of points, nothing adaptive occurs here
for ind = 1:4
    %get data for the GPR:
    w0                     = W(ii(ind));
    eval                   = f(w0); %evaluate at sample point
    deriv                  = df(w0); %get derivative of sample point
    x_data(ind,1)          = w0; 
    y_data(ind,1)          = (eval); 
    d_data(ind,1)          = deriv; 
    m_data                 = [y_data;d_data];
    plot(W,S1,'LineWidth',2.5,'color','k'); hold on; %plot the actual data
    [mean_m, cov_m, var_m] = plot_results(x_data,m_data,W,s,l,i,kernel_fun,kernel_der,kernel_hes,der_flag); %compute the GPR and plot
    xlabel('$x$','Interpreter','latex'); ylabel('$\sigma_{1}$','Interpreter','latex');
    set(gca,'LineWidth',2,'FontSize',20)
    hold off;
    pause(1)
end
%previous points were to initalize the method
%this gets the next points
for ind = 5:14
    %get data for the GPR:
    if choice_flag == 1
        %if you want adaptive scheme
        [x_next,f_acquis]      = acquisition_function(W,mean_m,sqrt(abs(var_m)),ind*kappa);
    else
        %if you want random points
        x_next = ii(ind);
    end
    if min(abs(x_next - x_data))<10^-2
        break
    else
    end
    w0                     = x_next;
    w0                     = W(ii(ind));
    eval                   = f(w0);     
    x_data(ind,1)          = w0; 
    y_data(ind,1)          = (eval); 
    d_data(ind,1)          = df(w0); 
    m_data                 = [y_data;d_data];
    figure(1)
    plot(W,S1,'LineWidth',2.5,'color','k'); hold on; 
    [mean_m, cov_m, var_m] = plot_results(x_data,m_data,W,s,l,i,kernel_fun,kernel_der,kernel_hes,der_flag);
    xlabel('$x$','Interpreter','latex'); ylabel('$f$','Interpreter','latex');
    set(gca,'LineWidth',2,'FontSize',20)
    hold off;
%     figure(2) %pots the acquisition function
%     [x_next,f_acquis]      = acquisition_function(W,mean_m,sqrt(abs(var_m)),ind*kappa);
%     plot(W,(f_acquis),'LineWidth',2.5,'color','k'); hold on; 
%     yticks([10^-2,10^-1,1,10,100,1000,10000]); %ylim([10^-1,max((f_acquis))])
%     xticks([-4:4]);%xlim([-4,4]);
%     %legend({'$\sigma_{1}$',['Guess ',num2str(j)]},'Interpreter','latex')
%     xlabel('$x$','Interpreter','latex'); ylabel('$f$','Interpreter','latex');
%     set(gca,'LineWidth',2,'FontSize',20)
%     set(gcf,'Position',[600,100,600,400]); 
%     hold off;
     pause(.1)
end

%% functions to use
%interpolant is commented out so not necessary
function [interp,xout] = interpolant_1(xdata,ydata,ddata,Npoints)
    x = (xdata-min(xdata))/range(xdata); 
    ddata = ddata*range(xdata);
    %get the coefficients
    d = ydata(1); c = ddata(1); a = -2*ydata(end) + 2*ydata(1) + ddata(1) + ddata(end); 
    b = ydata(end) - ydata(1) - ddata(1) - a; 
    %get the inerpolant
    X = linspace(0,1,Npoints);
    interp = polyval([a,b,c,d],X);
    xout = linspace(min(xdata),max(xdata),Npoints); 
end
%covariance matrix for the derivtives, includes derivative info
function cov_mat = covariance_matrix_der(x1,x2,s,l,kernel_fun,kernel_der,kernel_hes)
%covariance matrix with derivative info
    cov_mat1 = zeros(length(x1),length(x2));
    cov_mat2 = cov_mat1; 
    cov_mat3 = cov_mat1;
    for ind1 = 1:length(x1)
        for ind2 = 1:length(x2)
            cov_mat1(ind1,ind2) = kernel_fun(x1(ind1),x2(ind2),s,l);
            cov_mat2(ind1,ind2) = kernel_der(x1(ind1),x2(ind2),s,l);
            cov_mat3(ind1,ind2) = kernel_hes(x1(ind1),x2(ind2),s,l);
        end
    end
    cov_mat = [cov_mat1,-cov_mat2;cov_mat2,cov_mat3]; 
end
%covariance matrix without derivative info
function cov_mat = covariance_matrix(x1,x2,s,l,kernel_fun)
%covariance matrix without derivative info
    cov_mat1 = zeros(length(x1),length(x2));
    for ind1 = 1:length(x1)
        for ind2 = 1:length(x2)
            cov_mat1(ind1,ind2) = kernel_fun(x1(ind1),x2(ind2),s,l);
        end
    end
    cov_mat = cov_mat1; 
end
%GPR with derivative info
function [mean_at_values, cov_at_values, var_at_values] = GPR_der(xdata,ydata,x,s,l,i,kernel_fun,kernel_der,kernel_hes)
%set up the kernels, using derivatives
    K  = covariance_matrix_der(xdata,xdata,s,l,kernel_fun,kernel_der,kernel_hes);
    %Ki = inv(K + 10^-5*diag([zeros(length(xdata),1);ones(length(xdata),1)]));
    %the second term models noise and helps regularize the inversion
    Kl = chol(K + i*diag([zeros(length(xdata),1);ones(length(xdata),1)]));
    Kli= inv(Kl);
    Ks = covariance_matrix_der(x,xdata,s,l,kernel_fun,kernel_der,kernel_hes);
    Kss= covariance_matrix_der(x,x,s,l,kernel_fun,kernel_der,kernel_hes);
%get the mean
    %mean_at_values = Ks*Ki*ydata(:); 
    mean_at_values = Kli'*ydata(:);
    mean_at_values = Kli*mean_at_values;
    mean_at_values = Ks*mean_at_values; 
%get the covariance
    KsKli = Ks*Kli;
    %cov_at_values  = Kss - Ks*Ki*(Ks)'; cov_at_values = (cov_at_values + (cov_at_values)')/2 + 0*eps*eye(size(cov_at_values));
    cov_at_values  = Kss - KsKli*KsKli';  cov_at_values = (cov_at_values + (cov_at_values)')/2 + 0*eps*eye(size(cov_at_values));
    var_at_values  = diag(cov_at_values); 
    mean_at_values = mean_at_values(1:length(x));
    cov_at_values  = cov_at_values(1:length(x)); 
    var_at_values  = var_at_values(1:length(x)); 
end
%GPR without drivatives
function [mean_at_values, cov_at_values, var_at_values] = GPR(xdata,ydata,x,s,l,i,kernel_fun)
%set up the kernels, without using derivatives
    K  = covariance_matrix(xdata,xdata,s,l,kernel_fun);
    %the second term models noise and helps regularize the inversion
    Kl = chol(K + i*diag(zeros(length(xdata),1)));
    Kli= inv(Kl);
    Ks = covariance_matrix(x,xdata,s,l,kernel_fun);
    Kss= covariance_matrix(x,x,s,l,kernel_fun);
%get the mean
    %mean_at_values = Ks*Ki*ydata(:); 
    mean_at_values = Kli'*ydata(:);
    mean_at_values = Kli*mean_at_values;
    mean_at_values = Ks*mean_at_values; 
%get the covariance
    KsKli = Ks*Kli;
    %cov_at_values  = Kss - Ks*Ki*(Ks)'; cov_at_values = (cov_at_values + (cov_at_values)')/2 + 0*eps*eye(size(cov_at_values));
    cov_at_values  = Kss - KsKli*KsKli';  cov_at_values = (cov_at_values + (cov_at_values)')/2;
    var_at_values  = diag(cov_at_values); 
    mean_at_values = mean_at_values(1:length(x));
    cov_at_values  = cov_at_values(1:length(x)); 
    var_at_values  = var_at_values(1:length(x)); 
end

function [mean_m, cov_m, var_m] = plot_results(x_data,y_data,x,s,l,i,kernel_fun,kernel_der,kernel_hes,flag)
%flag determines if we use derivative info or not
    if flag == 1 %use derivative info
        disp('using derivative info')
        [mean_m, cov_m, var_m] = GPR_der(x_data,y_data,x,s,l,i,kernel_fun,kernel_der,kernel_hes);
    else %dont use derivative info
        disp('not using derivative info')
        [mean_m, cov_m, var_m] = GPR(x_data,y_data(1:length(x_data),:),x,s,l,i,kernel_fun);
    end
    std_m = sqrt(abs(real(var_m))); 
    x2 = [x,flip(x)];
    for ind1 = 1:5
        ind = 5-ind1+1;
        top = mean_m + ind*std_m; 
        bot = mean_m - ind*std_m; 
        inb = ([top;flip(bot)]);
        fill(x2(:),inb(:),[0.5569    0.5686    0.5608],'FaceAlpha',ind1*.1,'EdgeColor','none'); hold on;
    end
    hold on;
    colors = [0.8941    0.1020    0.1098;0.2157    0.4941    0.7216;0.3020    0.6863    0.2902;0.5961    0.3059    0.6392];%cbrewer('qual','Set1',4);%cbrewer('qual','Set1',4); 
    plot(x,(mean_m),'LineWidth',2,'Color',colors(1,:))
    hold on; 
    plot(x_data,(y_data(1:length(x_data))),'.','LineWidth',2,'MarkerSize',20,'Color',colors(2,:))
    plot(x_data(end),(y_data(length(x_data))),'+','LineWidth',2,'MarkerSize',8,'Color',colors(3,:))
    xlabel('$x$','Interpreter','latex'); 
    set(gca,'LineWidth',2,'FontSize',18)
end

function [x_next,f_acquis] = acquisition_function(x,mean,std,kappa)
    f_star    = max(mean); 
    f_acquis  = ((mean).*(std).^(1/6)) + std; %the acquisition function
    %raising kappa means we value the uncertainty more
    [~,index] = max(f_acquis); %the maximum of this
    x_next    = x(index); %the best point to go to 
end