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
    phase_mat = SHOcoef(:,:,:,:,4)+ offset;
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
    phase_mat(q_mat<5|q_mat>500)=NaN;
    amp_mat(q_mat<5|q_mat>500)=NaN;
    mixed_mat(q_mat<5|q_mat>500)=NaN;
    wres_mat(q_mat<5|q_mat>500)=NaN;
    q_mat(q_mat<5|q_mat>500)=NaN;
    
    phase_mat_w(q_mat_w<5|q_mat_w>500)=NaN;
    amp_mat_w(q_mat_w<5|q_mat_w>500)=NaN;
    mixed_mat_w(q_mat_w<5|q_mat_w>500)=NaN;
    wres_mat_w(q_mat_w<5|q_mat_w>500)=NaN;
    q_mat_w(q_mat_w<5|q_mat_w>500)=NaN;

 
    % The dimensions of SHO coefficients can be seen here:
    size(amp_mat)% Specify tolerance and voltage step and cycle of data used for masking
    tolerance = 0.1
    cycle = 2:3;    
    
    %v_dc Loop for one cycle 
    steps_loop=length(v_dc)/n_loops
    start_idx = (cycle-1)*steps_loop+1
    v_dc_loop = v_dc(start_idx:start_idx+steps_loop-1);

    feature_array1 = zeros(length(v_dc_loop),(n_row*n_col));
      feature_array2 = zeros(length(v_dc_loop),(n_row*n_col));
    
   
%     % Create a mask matrix of value 0, then set every pixel that is >
%     % pi-tolerance to 1
    mask = zeros(n_row, n_col);
    feature_array = zeros(length(v_dc_loop),(n_row*n_col));
    for i = 1:n_row
        for j = 1:n_col
        feature_array1(:,((i-1)*n_row+j)) = nanmean(mixed_mat(i,j,cycle,:),3);        
        end
    end
coeff1 = pca(feature_array1);
test1 = coeff1(:,1:5);

    feature_array = zeros(length(v_dc_loop),(n_row*n_col));
    for i = 1:n_row
        for j = 1:n_col
        feature_array2(:,((i-1)*n_row+j)) = nanmean(mixed_mat_w(i,j,cycle,:),3);        
        end
    end
coeff2 = pca(feature_array2);
test2 = coeff2(:,1:5);


   

    rng(1); % For reproducibility
coeffall = horzcat(test1,test2);
[idx,C] = kmeans(coeffall,4);
idxrot = rot90(idx);
    idxre = reshape(idxrot,[n_row,n_col]);
    mask = rot90(rot90(fliplr(rot90(idxre))));   
    
    coeffplot = vertcat(test1,test2);
    coeffplotx = coeffplot(:,1);
    coeffploty = coeffplot(:,2);
    
    %PCA 
    [coeff1,~,latent1,~,explained1] = pca(feature_array1);
    [coeff2,~,latent2,~,explained2] = pca(feature_array2);
    figure(6)
    plot((explained1(1:10))/100,'Marker','.','Color','k','LineWidth',1.5,'MarkerSize',30)
    hold on;
    plot((explained2(1:10))/100,'Marker','s','Color',[0.5,0.5,0.5],'LineWidth',1.5,'MarkerSize',8)
    legend('Off-F PR','ON-F PR', 'Location','northeast')
    legend boxoff
    xlabel('Principal Component Number')
    ylabel('Nomrmalized Eigenvalue')
    figure(7)
    [~,score1,~,] = pca(feature_array1');
    for i=1:6
        subplot(2,3,i)
        imagesc(reshape(abs(score1(:,i)/10),30,30));axis square; axis off; caxis ([0,1]),colorbar% 
    end
    figure(8)
    [~,score2,~,] = pca(feature_array2');
    for i=1:6
        subplot(2,3,i)
        imagesc(reshape(abs(score2(:,i)/15),30,30));axis square; axis off; caxis ([0,1]),colorbar% 
    end
    %coeff
    figure(9)
    clf
    plot(coeffplotx(idx==1,1),coeffploty(idx==1,1),'y.','MarkerSize',12)
hold on
    plot(coeffplotx(idx==2,1),coeffploty(idx==2,1),'r.','MarkerSize',12)
    plot(coeffplotx(idx==3,1),coeffploty(idx==3,1),'g.','MarkerSize',12)
    plot(coeffplotx(idx==4,1),coeffploty(idx==4,1),'b.','MarkerSize',12)
hold off

    %K-means map
 figure(10) 
    subplot(1,2,1)
    imagesc(mask)
    axis square
    axis off
    title('Mask')
    subplot(1,2,2)
    imagesc(idxre)
    axis square
    %caxis([-0.5 pi+0.5])
    title('amplitude data')
     mymap = [1 1 0
    1 1 1
    1 1 1
    1 1 1]
    colormap(mymap) 
 figure(11) 
    subplot(1,2,1)
    imagesc(mask)
    axis square
    axis off
    title('Mask')
    subplot(1,2,2)
    imagesc(idxre)
    axis square
    %caxis([-0.5 pi+0.5])
    title('amplitude data')
     mymap = [1 1 1
    1 0 0
    1 1 1
    1 1 1]
    colormap(mymap) 
  figure(12) 
    subplot(1,2,1)
    imagesc(mask)
    axis square
    axis off
    title('Mask')
    subplot(1,2,2)
    imagesc(idxre)
    axis square
    axis off
    %caxis([-0.5 pi+0.5])
    title('amplitude data')
     mymap = [1 1 1
    1 1 1
    0 1 0
    1 1 1]
    colormap(mymap) 
  figure(13) 
    subplot(1,2,1)
    imagesc(mask)
    axis square
    axis off
    title('Mask')
    subplot(1,2,2)
    imagesc(idxre)
    axis square
    %caxis([-0.5 pi+0.5])
    title('amplitude data')
     mymap = [1 1 1
    1 1 1
    1 1 1
    0 0 1]
    colormap(mymap) 
    
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