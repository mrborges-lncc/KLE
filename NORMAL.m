function NORMAL(Y,mu,std,tipo)
%LOGNORMAL    Create plot of datasets and fits
%   LOGNORMAL(Y)
%   Creates a plot, similar to the plot in the main distribution fitting
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with dfittool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  1
%   Number of fits:  1

% This function was automatically generated on 02-Sep-2011 14:34:34
 
% Data from dataset "Y data":
%    Y = Y

% Force all inputs to be column vectors
ymax = 0.4;

Y = Y(:);

% Set up figure to receive datasets and fits
f_ = 4;%clf;
figure(f_);

set(f_,'Units','Pixels','Position',[655 332 674 469.975]);
legh_ = []; legt_ = {};   % handles and text for legend
ax_ = newplot;
set(ax_,'Box','on','FontSize',20,'FontName','Times New Roman');
hold on;

% --- Plot data originally in dataset "Y data"
t_ = ~isnan(Y);
Data_ = Y(t_);
[F_,X_] = ecdf(Data_,'Function','cdf'...
              );  % compute empirical cdf
Bin_.rule = 3;
Bin_.nbins = 40;
[C_,E_] = dfswitchyard('dfhistbins',Data_,[],[],Bin_,F_,X_);
[N_,C_] = ecdfhist(F_,X_,'edges',E_); % empirical pdf from cdf
h_ = bar(C_,N_,'hist');
set(h_,'FaceColor','none','EdgeColor',[0.6 0.6 0.6],...
       'LineStyle','-', 'LineWidth',1);
xlabel(tipo,'FontSize',22,'Interpreter','latex');
ylabel('Density','FontSize',22,'Interpreter','tex')
legh_(end+1) = h_;
legt_{end+1} = 'computed';

% Nudge axis limits beyond data limits
xlim_ = get(ax_,'XLim');
if all(isfinite(xlim_))
   xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
   set(ax_,'XLim',xlim_)
end
%xlim_=[(mu-5*std) (mu+5*std)]
xlim_=[-5 5];
x_ = linspace(xlim_(1),xlim_(2),300);

% --- Create fit "lognormal"

% Fit this distribution to get parameter values
t_ = ~isnan(Y);
Data_ = Y(t_);
% To use parameter estimates from the original fit:
%     p_ = [ 0.3970417308125, 0.2607364479685];
pargs_ = cell(1,2);
[pargs_{:}] = normfit(Data_, 0.05);
p_ = [pargs_{:}];
y_ = normpdf(x_,p_(1), p_(2));
h_ = plot(x_,y_,'Color',[0 0 0],...
          'LineStyle','-', 'LineWidth',2,...
          'Marker','none', 'MarkerSize',6);
legh_(end+1) = h_;
%legt_{end+1} = ['$\mathbf{N}(' num2str(mu,'%1.2e') ',\,' num2str(std*std,'%1.2e') ')$'];
%legt_{end+1} = ['N(' num2str(mu,'%1.2e') ',' num2str(std*std,'%1.2e') ')'];
legt_{end+1} = ['N(' num2str(mu,'%1.2f') ', ' num2str(std*std,'%1.2f') ')'];

hold off;
% Create axes
% Uncomment the following line to preserve the Y-limits of the axes
% ylim([min(y_) max(y_)*1.6]);
ylim([0 ymax]);
xlim(xlim_);

% set(ax_,'DataAspectRatio',[1 (max(y_)*1.6-min(y_))/(max(x_)-min(x_)) 1]);
set(ax_,'DataAspectRatio',[1 (ymax-0)/(10) 1]);
leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast','Interpreter','latex'}; 
h_ = legend(ax_,legh_,legt_,leginfo_{:});  % create legend
set(h_,'FontSize',14,'Interpreter','tex','EdgeColor',[1 1 1]);