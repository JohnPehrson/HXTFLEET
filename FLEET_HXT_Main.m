clear all;close all;clc;
folderpath = "D:\HXT_FLEET";
scalefolder = "SCALES";
datapaths = ["apr1";"apr1_second";"apr2";"apr4";"apr5";"apr6";"apr7";"apr9";"apr10";"apr11"];
prerun_foldernames = ["nitro3torr";"nitro3torr";"nitro3torr";"8.6TorrAir"; "air_8_6torr";"air_x_torr";"air";"3torrair";"air_2Torr";"apr11_air_2Torr"];

%% Reference image processing and image scales

% folder_folderpath = fullfile(folderpath,scalefolder);
% [filteredfiles] = Folders_in_Folder(folder_folderpath,"","");
% 
% for i = 1:length(filteredfiles)
%     folder_folderpath_days = fullfile(folderpath,scalefolder,filteredfiles(i));
%     [filteredfiles_days] = Folders_in_Folder(folder_folderpath_days,"","");
% 
%     for j = 1:length(filteredfiles_days)
%         %If I've already got a time-average
%         folder_folderpath_ims = fullfile(folder_folderpath_days,filteredfiles_days(j));
%         have = ["mat"];
%         [filteredfiles_ims] = Folders_in_Folder(folder_folderpath_ims,have,"");
% 
%         av_filepath = fullfile(folder_folderpath_ims,"AVERAGE.mat");
%         if ~isempty(filteredfiles_ims)
%             load(av_filepath);
%         else %If I haven't already got a time-average
%             folder_folderpath_ims = fullfile(folder_folderpath_days,filteredfiles_days(j));
%             canthave = ["rec","zip","mat","ave"];
%             [filteredfiles_ims] = Folders_in_Folder(folder_folderpath_ims,"",canthave);
%             
%             max_ims = min(length(filteredfiles_ims),400);
%             scale_ims = zeros(2048,2048,max_ims);
%         for k = 1:max_ims
%                     fileNameCat = fullfile(folder_folderpath_ims,filteredfiles_ims(k));
%                     scale_ims(:,:,k) = double(imread(fileNameCat));
%         end
%             av_scale = mean(scale_ims,3);
%             clear scale_ims
%             save(av_filepath,'av_scale')
%         end
% 
%             figure;
%             image(av_scale)
%             colorbar;
%             colormap(turbo(max(av_scale(:))));
%             grid on;
%             axis equal;
%             xlim([0,2048]);
%             ylim([600,1600]);
%             title(strcat(filteredfiles(i), ":", filteredfiles_days(j)))
% 
% 
%             disp('Pick an initial and displaced point')
%             [x,y] = ginput(2);
%             dx = x(2)-x(1);
%             dy = y(2)-y(1);
%             total_dist = sqrt(dx.^2+dy.^2);
% 
%             prompt = "Number of mm between the two points";
%             mm = input(prompt);
%             imagescale = total_dist/mm; %pix/mm
% 
%     end
% 
% end

images_scales = [21.171;21.58;21.2185;21.17;21.271]; %pix/mm 
g3g4_ave = mean([images_scales(3),images_scales(4)]);

run_scales = [images_scales(1);images_scales(1);images_scales(1);...
              images_scales(2);g3g4_ave;g3g4_ave;g3g4_ave;g3g4_ave;g3g4_ave;images_scales(5)]; %pix/mm

imaging_settings = [1500,150;...
                    1500,150;...
                    1500,150;...
                    1500,150;...
                    2000,200;...
                    3000,200;...
                    1500,160;...
                    1500,160;...
                    1500,160;...
                    1500,160;]; %[Delay, gate] with rows of runs. The delay time is given in ns to the center of the gate also given in ns.

% Load Prerun Images
prerun_row_averages = zeros(length(datapaths),2048); %indexes are run number, column

for i = 1:length(datapaths)

    % Prerun
        %get folder name
        run_filepath = fullfile(folderpath,datapaths(i));
        nothave = ["run","readme"];
        [filteredfiles_days] = Folders_in_Folder(run_filepath,"",nothave);
        prerun_foldernames(i) = filteredfiles_days;
        prerun_filefolder = fullfile(run_filepath,filteredfiles_days);
        
        av_prerun_filepath = fullfile(prerun_filefolder,"AVERAGE.mat");
        have = ["mat"];
        [filteredfiles_ims] = Folders_in_Folder(prerun_filefolder,have,"");
     if ~isempty(filteredfiles_ims)
        load(av_prerun_filepath);
     else %If I haven't already got a time-average

        %get files in folder
        nothave = ["rec","mat","u"];
        [prerun_files] = Folders_in_Folder(prerun_filefolder,"",nothave);

        max_ims = min(length(prerun_files),2000);

        prerun_ims = zeros(200,2048,max_ims);
        prerun_row_averages = zeros(length(datapaths),2048); %indexes are run number, column
        for k = 1:max_ims
            fileNameCat = fullfile(prerun_filefolder,prerun_files(k));
            prerun_ims(:,:,k) = double(imread(fileNameCat));
        end
        prerun_ave = mean(prerun_ims,3);
        row_sum_prerun = mean(prerun_ave,1);
        clear prerun_ims

        save(av_prerun_filepath,'prerun_ave','row_sum_prerun');
     end

        %get rid of background 
        prerun_row_averages(i,:) = row_sum_prerun-min(row_sum_prerun);

%             figure;
%             image(prerun_ave)
%             colorbar;
%             colormap(turbo(max(prerun_ave(:))));
%             grid on;
%             axis equal;
%             title(["Pre-run average from ",prerun_filefolder])
end

%% Display Prerun images
% for i = 1:length(datapaths)
%     filepath = fullfile(folderpath,datapaths(i),prerun_foldernames(i),"AVERAGE.mat");
%     load(filepath);
% 
%         figure;
%         subplot(2,1,1);
%         image(prerun_ave)
%         colorbar;
%         colormap(turbo(max(prerun_ave(:))));
%         grid on;
%         axis equal;
%         title(["Pre-run average from ",datapaths(i)])
%         subplot(2,1,2);
%         plot(1:length(prerun_row_averages(i,:)),prerun_row_averages(i,:));
% 
% end

% Load Run Images
        run_ims_count = 11;
        run_ims = zeros(200,2048,run_ims_count);
        run_row_averages = zeros(length(datapaths),2048,run_ims_count);
for i = 1:length(datapaths)
        run_filepath = fullfile(folderpath,datapaths(i),"run");
        have = ["tif"];
        nothave = ["rec"];
        [filteredfiles_run] = Folders_in_Folder(run_filepath,have,nothave);
        process_runs = filteredfiles_run(end-(run_ims_count-1):end);
    
%         figure;
        for k = 1:length(process_runs)
            fileNameCat = fullfile(run_filepath,process_runs(k));
            run_ims(:,:,k) = double(imread(fileNameCat));
            run_row_averages(i,:,k) = mean(run_ims(:,:,k),1);
% 
%             subplot(2,1,1);
%             image(run_ims(:,:,k))
%             colorbar;
%             colormap(turbo(max(max(run_ims(:,:,k)))));
%             grid on;
%             axis equal;
%             title(['Run number', num2str(i), process_runs(k)])
%             subplot(2,1,2);
%             plot(1:length(run_row_averages(i,:,k)),run_row_averages(i,:,k));
        end

end

%% Identify displacements in units of pixels
% close all;
% f1 = figure;
% f2 = figure;
% 
% disp_ims = zeros(length(datapaths),size(run_row_averages,3));
% disp_unc_ims = zeros(length(datapaths),size(run_row_averages,3));
% 
%     %cross correlation
%     for i = 1:length(datapaths) %runs
%         for j = 1:size(run_row_averages,3) %images from each run
%         prerun = prerun_row_averages(i,:);
%         run = run_row_averages(i,:,j);
%         col_length = length(prerun);
%         corr_step = (-1*(col_length-1)):(col_length-1);
% 
%         if min(run(500:1000)) <65000
% 
%             %normalization
%             prerun = prerun./max(prerun);
%             run = run-min(run(500:1000));
%             run = run./max(run);
%             corr = xcorr(prerun,run);
% 
%             %plotting
%             figure(f1);
%             subplot(2,1,1);
%             plot(1:col_length,prerun);
%             title(strcat("Prerun average from run ",num2str(i)));
%             subplot(2,1,2);
%             plot(1:col_length,run);
%             title(strcat("Run ",num2str(i)," image number ",num2str(j)));
% 
%             figure(f2);
%             subplot(2,1,1);
%             plot(corr_step,corr);
%             xlim([-200,200]);
%             title(strcat("Correlation from Run ",num2str(i)," image number ",num2str(j)));
%             xline(0)
%             subplot(2,1,2);
%             plot(corr_step,corr);
%             title(strcat("Correlation from Run ",num2str(i)," image number ",num2str(j)));
%             xline(0)
% 
%         %Pre-fitting cropping and backgrounding
%             pix_width = 2048;
%             t0line = 575;
%             region_background = (-pix_width+t0line):t0line-1;
%             Lia = ismember(corr_step,region_background);
% 
%             crop_corr = corr(Lia);
%             crop_step = corr_step(Lia);
% 
%         %Fitting
%             lims = [2*max(corr),1,200,50,100,max(corr);...
%                     max(corr),0.5,-50,25,10,median(crop_corr);...
%                     max(corr)/2,0,-200,5,2,median(crop_corr)/2]; %[h,n,x0,sigma,R,bkg];
% 
%         [fitvariables,outfit,fit_unc] = MaxFitting(crop_corr,crop_step,lims);
%             figure(f2);
%             subplot(2,1,1);
%             hold on;
%             plot(crop_step,outfit,'r','Linewidth',2);
%             xlim([-200,50]);
%             hold off;
% 
%             %residual near the fitting location
%             region_near_max = round(fitvariables(3));
%             region_fit = (region_near_max-50):(region_near_max+50);
%             Lia = ismember(crop_step,region_fit);
%             fit_region = crop_step(Lia);
%             fit_corr = crop_corr(Lia);
%             fit_outfit = outfit(Lia);
%             resid = sum(abs(fit_corr-fit_outfit));
% 
%             resid_thresh = 50;
%             if resid>resid_thresh
%                 x = input("Can you identify a peak point accurately?");
%                 if strcmp(x,"y")
%                     disp(strcat("Residual was ", num2str(resid)," so manually pick estimate of point height"))
%                     %manually pick
%                     figure(f2);
%                     subplot(2,1,1);
%                     roi = drawpoint('Color','r');
%                     displacement = -roi.Position(1);
% 
%                     roi = drawline('Color','k');
%                     displacement_unc = abs(roi.Position(1)-roi.Position(2))/2;
%                 else
%                     displacement = NaN;
%                     displacement_unc = NaN;
%                 end
%   
%             else %use fitted data, send data out
%                 disp(strcat("Residual was ", num2str(resid)," so use fitted displacement"))
%                 displacement = -fitvariables(3);
%                 displacement_unc = 2.*fit_unc;
% 
%             end
%         disp_ims(i,j) = displacement;
%         disp_unc_ims(i,j) = displacement_unc;
% 
%         else
%         disp_ims(i,j) = NaN;
%         disp_unc_ims(i,j) = NaN;
%         end
% 
%         end
%     end
% 
%     save_disp_name = "Displacements.mat";
%     save(save_disp_name,'disp_ims','disp_unc_ims');


%% Velocity Calculation
%run_scales
%imaging_settings
% disp_ims disp_unc_ims

%load in relevant displacements
save_disp_name = "Displacements.mat";
load(save_disp_name);

%conversions to standard units
run_scales_m = run_scales.*1000; %pix/m 
delays = imaging_settings./10^9; %ns

%velocity calculation
velo = zeros(size(disp_ims));
for i = 1:length(datapaths)
s = run_scales_m(i);
t = delays(i,1);
velo(i,:) = disp_ims(i,:)./(s.*t);

end

%velocity uncertainty calculation
velo_unc = zeros(size(disp_ims));
for i = 1:length(datapaths)
    for j= 1:size(disp_ims,2)
        s = run_scales_m(i);
        del_s = run_scales_m(i).*0.01;
        t = delays(i,1)-delays(i,2)./2;
        del_t = delays(i,2)./4;
        x = disp_ims(i,j);
        del_x = disp_unc_ims(i,j);

        velo_unc(i,j) = sqrt(((del_x)/(t*s)).^2+((x*del_t)/(s*(t^2))).^2+((x*del_s)/(t*(s^2))).^2);
    end
end

%velo
%velo_unc
timing = linspace(1,size(disp_ims,2),size(disp_ims,2));

%% Plotting
marker_list = ["o";"*";"square";"diamond";"^";"v";"<";">";"pentagram";"hexagram"];
datapaths_legend = datapaths;
datapaths_legend(2) = "apr1 run 2";
datapaths_legend = strcat("Run on ",datapaths_legend);

figure;
for i= 1:length(datapaths)
    e = errorbar(timing,velo(i,:),velo_unc(i,:),'o');
    e.Marker = marker_list(i);
    xlabel('Processed Frames [$\mu$s]','Interpreter','Latex')
    ylabel('Freestream Velocity [m/s]')
    hold on;
end

legend(datapaths_legend)
xlim([min(timing)-0.5,max(timing)+0.5]);
grid on;
set(gca,'FontSize', 18);
set(gca,'fontname','times')  % Set it to times

%% Saving

velo_save = round(velo, 1);
velo_unc_save = round(velo_unc, 1);

csv_filename = 'FLEET_HXT_Velocities.txt';
writematrix(velo,csv_filename,'Delimiter','tab')
csv_filename = 'FLEET_HXT_VelocityUncertainties.txt';
writematrix(velo_unc_save,csv_filename,'Delimiter','tab')


















