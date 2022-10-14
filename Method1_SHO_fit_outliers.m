%  function SHOfit_new_guess
% close all 
clear all

global n_wbin n_dc_steps n_loops n_row n_col w_vec v_dc
global save_path


load_path='***\';
load_name='****';
showFigures=1;
show_step=50;

save_path=['****\']; 

%mkdir([load_path 'Analysis\Figures_read']);
%mkdir([load_path 'Analysis\Figures_write']);
%save_path=[load_path 'Analysis\'];

n_fig=1;

max_fun_evals=200;
max_iterations=200;
tolerance=1e-9;


fitting_function='SHO_fit_01';
options = optimset('MaxFunEvals',max_fun_evals,'MaxIter',max_iterations,'TolFun',tolerance,'display','off');


import_param(load_path, load_name)
data=import_data_read(load_path, load_name);
%data=import_data_write(load_path, load_name);

size(data)
data=reshape(data,n_wbin,n_dc_steps./n_loops,n_loops,n_col,n_row);

data=permute(data,[2 1 3 4 5]);  
size(data)
w_vec_p=w_vec/1000; % to display kHz

row = 12;
col = 2;
loop = 3;
dc_step = 7;
%dc_step = 7;
%dc_step = 9;



       
            response_vec=data(dc_step,:,loop,col,row);
            response_vec = detrend(response_vec); %linear background removal
            Q_guess = 120;
            [A_max,wr_ind] = max(abs(response_vec));
            A_guess = A_max/Q_guess;
            wr_guess = w_vec(wr_ind);
            phi_guess = angle(response_vec(wr_ind))+pi/2;
            
            w_lb=w_vec(1)*0.9;
            w_ub=w_vec(end)*1.1;
            
            lower_b=[0.1 w_lb 1 -pi+0.05];
            upper_b=[inf w_ub 500 +pi];
            
            ig=[A_guess, wr_guess, Q_guess, phi_guess];
                   
            H_data=[real(response_vec); imag(response_vec)];
           
            [fp, resnorm, residual,exitflag,output]=lsqcurvefit(eval(['@' fitting_function]),ig,w_vec,H_data,lower_b,upper_b,options);
            
           
            
                    [row col loop dc_step]
                    fit_data_h=feval(fitting_function,fp,w_vec);
                    fit_data=fit_data_h(1,:)+i*fit_data_h(2,:);
                    
                    guess_data_h=feval(fitting_function,ig,w_vec);
                    guess_data=guess_data_h(1,:)+i*guess_data_h(2,:);
                    
                    fit_amp=abs(fit_data);
                    fit_phase=angle(fit_data);
                    guess_amp=abs(guess_data);
                    guess_phase=angle(guess_data);

                    %**********************************************FIGURES************************************************************
                    flim=[min(w_vec_p)*0.99 max(w_vec_p)*1.01];
                    
                    fh=figure(48);
                    clf
                    fh.Color='White';
                   annotation('textbox', [0 0.9 1 0.1], ...
                        'String', ['R ' num2str(row) ' C ' num2str(col) ' L ' num2str(loop) ' V ' num2str(v_dc(dc_step))  ' Q = ' num2str(fp(3))], ...
                        'EdgeColor', 'none', ...
                        'HorizontalAlignment', 'center', 'FontSize',8)
                                        
                   subplot(2,1,1)
                    plot(w_vec_p,fit_amp,'-b.')
                    hold on
                    plot(w_vec_p,abs(response_vec),'rx')
                    %plot(w_vec_p,guess_amp,'-.g')
                    xlabel('\omega [kHz]')
                    ylabel('R [a.u.]')
                    xlim(flim);
                                        
                    subplot(2,1,2)
                    plot(w_vec_p,fit_phase,'-b.')
                    hold on
                    plot(w_vec_p,angle(response_vec),'rx')
                    %plot(w_vec_p,guess_phase,'-.g')
                    xlabel('\omega [kHz]')
                    ylabel('\phi [rad]')
                    ylim([-2*pi 2*pi]);
                    xlim(flim);
                                        
%                     subplot(2,2,2)
%                     plot(w_vec_p,real(fit_data),'-b.')
%                     hold on
%                     plot(w_vec_p,real(response_vec),'rx')
%                     %plot(w_vec_p,real(guess_data),'-.g')
%                     xlabel('\omega [kHz]')
%                     ylabel('Real [a.u.]')
%                     xlim(flim);
%                     
%                     
%                     subplot(2,2,4)
%                     plot(w_vec_p,imag(fit_data),'-b.')
%                     hold on
%                     plot(w_vec_p,imag(response_vec),'rx')
%                     %plot(w_vec_p,imag(guess_data),'-.g')
%                     xlabel('\omega [kHz]')
%                     ylabel('Imaginary [a.u.]')
%                     xlim(flim);
                    
                    print(fh,[ save_path 'Row = ' num2str(row) ' Col = ' num2str(col)] ,'-dpng')
                
%                 figure(2)
%                 yyaxis left
%                 plot(w_vec_p,fit_amp,'-k.')
%                 hold on
%                 plot(w_vec_p,abs(response_vec),'rx')
%                 xlabel('\omega [kHz]')
%                 ylabel('R [a.u.]')
%                 yyaxis right
%                  plot(w_vec_p,fit_phase,'-k.')
%                     hold on
%                     plot(w_vec_p,angle(response_vec),'bx')
%                     %plot(w_vec_p,guess_phase,'-.g')
%                     ylabel('\phi [rad]')
%                     ylim([-2*pi 2*pi]);
%                                     xlim(flim);

                %if(mod(n_fig,show_step)==0) 
                 %   [row col loop dc_step]
                %end
            %end
            
            fit_mat(col, row, loop, dc_step,:)=fp;
            ef_mat(col, row, loop, dc_step,:)=exitflag;
            %n_fig=n_fig+1;
       
   %SHOcoeff = squeeze(fit_mat(4,3,2,14,:));

%fid=fopen([save_path load_name '_SHOcoeff_read.dat'],'w');
%fwrite(fid,fit_mat,'real*4');
%fclose(fid);

%fid=fopen([save_path load_name '_Exit_flags_read.dat'],'w');
%fwrite(fid,ef_mat,'real*4');
%fclose(fid);





function H = sho_IQ_fit_01(a,w)

    w = -w;
    A = a(1);
    wo = a(2);
    Q = a(3);
    phi = a(4);

    Omega = w.^2 - wo.^2;
    Xi = (w.*wo)./Q;

    denom = (Omega.^2 + Xi.^2);
    H_I = A.*(wo.^2).*(Omega*cos(phi)-Xi*sin(phi))./denom;
    H_Q = A.*(wo.^2).*(Xi*cos(phi)+Omega*sin(phi))./denom;

    H = [H_I; H_Q];

end


function H_mat = SHO_fit_01(a,w)
    A = a(1);
    wo = a(2);
    Q = a(3);
    phi = a(4);


    H = A*exp(1i*phi).*wo.^2./( w.^2 + 1i*w.*wo/Q - wo.^2 );

    H_real = real(H);
    H_imag = imag(H);

    H_mat = [ H_real ; H_imag ];

end

function [data]= import_data_read(load_path, load_name)
    global i_read r_read

    file_id = fopen([load_path load_name '_F0_read_imag.dat'], 'r');
    i_read = fread(file_id,'real*4');
    fclose(file_id);

    file_id = fopen([load_path load_name '_F0_read_real.dat'], 'r');
    r_read = fread(file_id,'real*4');
    fclose(file_id);

    data=[r_read+i*i_read];

end

function [data]= import_data_write(load_path, load_name)
    global i_read r_read

    file_id = fopen([load_path load_name '_F0_write_imag.dat'], 'r');
    i_read = fread(file_id,'real*4');
    fclose(file_id);

    file_id = fopen([load_path load_name '_F0_write_real.dat'], 'r');
    r_read = fread(file_id,'real*4');
    fclose(file_id);

    data=[r_read+i*i_read];

end
   
function import_param(file_path, file_name)
    global n_wbin n_dc_steps n_row n_col w_vec v_dc n_loops

    param_cell=load([file_path file_name]);

    w_vec=param_cell.bin_w;
    size_bin_ind=size(w_vec);
    n_wbin=size_bin_ind(2)

    v_dc=param_cell.dc_amp_vec_full;
    n_dc_steps=length(v_dc)

    n_loops=param_cell.SS_parm_vec(3)

    % size_AI2mat=size(param_cell.AI2_read_mat3);
    n_row=param_cell.position_vec(3)
    n_col=param_cell.position_vec(4)

end