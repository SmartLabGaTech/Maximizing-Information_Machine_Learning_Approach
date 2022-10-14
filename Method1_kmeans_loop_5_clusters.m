% function mean_of_the_fullyswitched_read

    close all
    clear all

    global n_wbin n_dc_steps n_row n_col w_vec v_dc n_loops Filename 
    global save_path 
    
 % Where (which folder) is the file is on your computer
   load_path='****\';
    
    % Name of main .mat file
    load_name='****';
    
%       % Where (which folder) is the file is on your computer
%    load_path='C:\Users\User\Documents\PTP Samples\Batch 2\Batch 2 S 1-2\190618\BEPS_0005_d\';
%     
%     % Name of main .mat file
%     load_name='BEPS_0005';
    
    mkdir([load_path 'BEPS_IMAGES\kmeans_loops']);
    save_path=[load_path 'BEPS_IMAGES\kmeans_loops'];


    %This function imports important parameters such as number of voltage
    %steps, number of rows and cols,....
    import_param(load_path,load_name);

    % This function imports SHO coefficients phase, amplitude, resonance
    % frequency and quality factor for all pixels at all voltage steps
    SHOcoef = import_SHOcoeff(load_path,[load_name '_SHOcoeff_read.dat']);
    SHOcoef = permute(SHOcoef,[2,1,3,4,5]);
    SHOcoefwrite = import_SHOcoeff(load_path,[load_name '_SHOcoeff_write.dat']);
    SHOcoefwrite = permute(SHOcoefwrite,[2,1,3,4,5]);


    % SHO coefficients phase, amplitude, resonance
    % frequency and quality factor are assigned to arrays
    % the dimensions are e.g. phase_mat(row, column, cycle, dc voltage step)
    offset = -1.55; %phase offset should be selected so that the resulting phase is between 0 and pi
    phase_mat = SHOcoef(:,:,:,:,4) + offset;
    amp_mat = SHOcoef(:,:,:,:,1);
    wres_mat = SHOcoef(:,:,:,:,2); 
    q_mat = SHOcoef(:,:,:,:,3);
    mixed_mat = amp_mat.*cos(phase_mat);
    
    phase_mat_w = SHOcoefwrite(:,:,:,:,4)+offset ;
    amp_mat_w = SHOcoefwrite(:,:,:,:,1);
    wres_mat_w = SHOcoefwrite(:,:,:,:,2); 
    q_mat_w = SHOcoefwrite(:,:,:,:,3);
    mixed_mat_w = amp_mat_w.*cos(phase_mat_w);
    
    %set low fit point as nan
    phase_mat(q_mat<5|q_mat>1000)=NaN;
    amp_mat(q_mat<5|q_mat>1000)=NaN;
    mixed_mat(q_mat<5|q_mat>1000)=NaN;
    wres_mat(q_mat<5|q_mat>1000)=NaN;
    q_mat(q_mat<5|q_mat>1000)=NaN;
    
    phase_mat_w(q_mat_w<5|q_mat_w>1000)=NaN;
    amp_mat_w(q_mat_w<5|q_mat_w>1000)=NaN;
    mixed_mat_w(q_mat_w<5|q_mat_w>1000)=NaN;
    wres_mat_w(q_mat_w<5|q_mat_w>1000)=NaN;
    q_mat_w(q_mat_w<5|q_mat_w>1000)=NaN;

 
    % The dimensions of SHO coefficients can be seen here:
    size(amp_mat)% Specify tolerance and voltage step and cycle of data used for masking
    tolerance = 0.1
    cycle = 2:3
    voltage_step = 40
    
    
    %v_dc Loop for one cycle 
    steps_loop=length(v_dc)/n_loops
    start_idx = (cycle-1)*steps_loop+1
    v_dc_loop = v_dc(start_idx:start_idx+steps_loop-1);


    % Create the feature array for mixed read
    feature_array = zeros(length(v_dc_loop),(n_row*n_col));
    for i = 1:n_row
        for j = 1:n_col
        feature_array(:,((i-1)*n_row+j)) = nanmean(mixed_mat(i,j,cycle,:));        
        end
    end
coeff1 = pca(feature_array); %carry out principal component analysis
test1 = coeff1(:,1:5); 

    % Create the feature array for mixed write

    feature_array = zeros(length(v_dc_loop),(n_row*n_col));
    for i = 1:n_row
        for j = 1:n_col
        feature_array(:,((i-1)*n_row+j)) = nanmean(mixed_mat_w(i,j,cycle,:));        
        end
    end
coeff2 = pca(feature_array); %carry out principal component analysis
test2 = coeff2(:,1:5); 



   

    rng(1); % For reproducibility
coeffall = horzcat(test1,test2);
[idx,C] = kmeans(coeffall,5);
idxrot = rot90(idx);
    idxre = reshape(idxrot,[n_row,n_col]);
    mask = rot90(rot90(fliplr(rot90(idxre))));   
    
    coeffplot = vertcat(test1,test2);
    coeffplotx = coeffplot(:,1);
    coeffploty = coeffplot(:,2);


    figure(9)
    clf
    plot(coeffplotx(idx==1,1),coeffploty(idx==1,1),'r.','MarkerSize',12)
hold on
    plot(coeffplotx(idx==2,1),coeffploty(idx==2,1),'y.','MarkerSize',12)
    plot(coeffplotx(idx==3,1),coeffploty(idx==3,1),'g.','MarkerSize',12)
    plot(coeffplotx(idx==4,1),coeffploty(idx==4,1),'b.','MarkerSize',12)
    plot(coeffplotx(idx==5,1),coeffploty(idx==5,1),'k.','MarkerSize',12)
hold off

cycle = 2;

 figure(10) 
    subplot(1,1,1)
    imagesc(mask)
    axis square
    title('Mask')
    axis off
    mymap = [1 1 0
    1 0 0
    0 1 0
    0 0.5 0.75
    0 0 1]
    colormap(mymap) 
     
    % Mask is 2 dimensional so here phase data is sliced up

    % Mask is 2 dimensional so here mixed piezoresponse data is sliced up
for l = 1:n_loops
    for s = 1:n_dc_steps/n_loops
        mixed_slice = squeeze(mixed_mat(:,:,l,s));
        mixed_masked_1_mean(l,s) = mean(mixed_slice(mask == 1),1,'omitnan');
        mixed_masked_1_std(l,s) = std(mixed_slice(mask == 1),0,1,'omitnan');
        mixed_masked_2_mean(l,s) = mean(mixed_slice(mask == 2),1,'omitnan');
        mixed_masked_2_std(l,s) = std(mixed_slice(mask == 2),0,1,'omitnan');
        mixed_masked_3_mean(l,s) = mean(mixed_slice(mask == 3),1,'omitnan');
        mixed_masked_3_std(l,s) = std(mixed_slice(mask == 3),0,1,'omitnan');
        mixed_masked_4_mean(l,s) = mean(mixed_slice(mask == 4),1,'omitnan');
        mixed_masked_4_std(l,s) = std(mixed_slice(mask == 4),0,1,'omitnan');
        mixed_masked_5_mean(l,s) = mean(mixed_slice(mask == 5),1,'omitnan');
        mixed_masked_5_std(l,s) = std(mixed_slice(mask == 5),0,1,'omitnan');
    end
end
 % Mask is 2 dimensional so here mixed piezoresponse data is sliced up
for l = 1:n_loops
    for s = 1:n_dc_steps/n_loops
        mixed_slice_w = squeeze(mixed_mat_w(:,:,l,s));
        mixed_masked_1_mean_w(l,s) = mean(mixed_slice_w(mask == 1),1,'omitnan');
        mixed_masked_1_std_w(l,s) = std(mixed_slice_w(mask == 1),0,1,'omitnan');
        mixed_masked_2_mean_w(l,s) = mean(mixed_slice_w(mask == 2),1,'omitnan');
        mixed_masked_2_std_w(l,s) = std(mixed_slice_w(mask == 2),0,1,'omitnan');
        mixed_masked_3_mean_w(l,s) = mean(mixed_slice_w(mask == 3),1,'omitnan');
        mixed_masked_3_std_w(l,s) = std(mixed_slice_w(mask == 3),0,1,'omitnan');
        mixed_masked_4_mean_w(l,s) = mean(mixed_slice_w(mask == 4),1,'omitnan');
        mixed_masked_4_std_w(l,s) = std(mixed_slice_w(mask == 4),0,1,'omitnan');
        mixed_masked_5_mean_w(l,s) = mean(mixed_slice_w(mask == 5),1,'omitnan');
        mixed_masked_5_std_w(l,s) = std(mixed_slice_w(mask == 5),0,1,'omitnan');
    end
end

    

    figure(11)
    clf
    plot(v_dc(1:steps_loop), mixed_masked_1_mean(cycle, :),'k.-','LineWidth',1,'MarkerSize',20)
    hold on
    plot(v_dc(1:steps_loop), mixed_masked_1_mean_w(cycle, :),'Marker','s','MarkerFaceColor',[0.5,0.5,0.5],'color',[0.5,0.5,0.5],'LineWidth',1.5,'MarkerSize',6)
    xlim([-10, 10])
% legend('Read',' Write', 'Location','best')
%     legend boxoff
%     xlabel('DC Voltage (V)')
%     ylabel('Average PR (a.u.)')
    ylim ([-10,10])
    
    set(gca, 'FontSize', 14)
    set(gcf,'unit','normalized','position',[0.2,0.2,0.5,0.5])


   
    
    
     
    
    %Type 2
    
    figure(21)
    clf
    plot(v_dc(1:steps_loop), mixed_masked_2_mean(cycle, :),'k.-','LineWidth',1,'MarkerSize',20)
    hold on
    plot(v_dc(1:steps_loop), mixed_masked_2_mean_w(cycle, :),'Marker','s','MarkerFaceColor',[0.5,0.5,0.5],'color',[0.5,0.5,0.5],'LineWidth',1.5,'MarkerSize',6)
    xlim([-10, 10])
% legend('Read',' Write', 'Location','best')
%     legend boxoff
%     xlabel('DC Voltage (V)')
%     ylabel('Average PR (a.u.)')
    ylim ([-10.5,10])
    set(gca, 'FontSize', 14)
    set(gcf,'unit','normalized','position',[0.2,0.2,0.5,0.5])
    

   
    
   %Type 3
    
    figure(31)
    clf
    plot(v_dc(1:steps_loop), mixed_masked_3_mean(cycle, :),'k.-','LineWidth',1,'MarkerSize',20)
    hold on
    plot(v_dc(1:steps_loop), mixed_masked_3_mean_w(cycle, :),'Marker','s','MarkerFaceColor',[0.5,0.5,0.5],'color',[0.5,0.5,0.5],'LineWidth',1.5,'MarkerSize',6)
    xlim([-10, 10])
% legend('Read',' Write', 'Location','best')
%     legend boxoff
%     xlabel('DC Voltage (V)')
%     ylabel('Average PR (a.u.)')
    ylim ([-10,10])
    set(gca, 'FontSize', 14)
    set(gcf,'unit','normalized','position',[0.2,0.2,0.5,0.5])
    

   
  
     
    
      %Type 4
    
    figure(41)
    clf
    plot(v_dc(1:steps_loop), mixed_masked_4_mean(cycle, :),'k.-','LineWidth',1,'MarkerSize',20)
    hold on
    plot(v_dc(1:steps_loop), mixed_masked_4_mean_w(cycle, :),'Marker','s','MarkerFaceColor',[0.5,0.5,0.5],'color',[0.5,0.5,0.5],'LineWidth',1.5,'MarkerSize',6)
    xlim([-10, 10])
%     legend('Read',' Write', 'Location','best')
%     legend boxoff
%     xlabel('DC Voltage (V)')
%     ylabel('Average PR (a.u.)')
    ylim ([-10,10])
    set(gca, 'FontSize', 14)
    set(gcf,'unit','normalized','position',[0.2,0.2,0.5,0.5])
 

   
    %Type 5
    
    figure(51)
    clf
    plot(v_dc(1:steps_loop), mixed_masked_5_mean(cycle, :),'k.-','LineWidth',1,'MarkerSize',20)
    hold on
    plot(v_dc(1:steps_loop), mixed_masked_5_mean_w(cycle, :),'Marker','s','MarkerFaceColor',[0.5,0.5,0.5],'color',[0.5,0.5,0.5],'LineWidth',1.5,'MarkerSize',6)
    xlim([-10, 10])
%     legend('Read',' Write', 'Location','best')
%     legend boxoff

   
%     xlabel('DC Voltage (V)')
%     ylabel('Average PR (a.u.)')
    ylim ([-10,10])
    set(gca, 'FontSize', 14)
    set(gcf,'unit','normalized','position',[0.2,0.2,0.5,0.5])
    
     
     

    
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

function [mat] = import_SHOcoeff(file_path, file_name)
    global n_wbin n_dc_steps n_row n_col w_vec v_dc n_loops
    [file_path file_name]
    fid=fopen([file_path file_name]);
    mat=fread(fid,'real*4');
    fclose(fid);

    size(mat)

    mat=reshape(mat, n_row, n_col, n_loops, length(v_dc)./n_loops, 4);

end