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

images_scales = [21.171;21.58;21.2185;21.17;21.271]; %pix/mm for 
g3g4_ave = mean([images_scales(3),images_scales(4)]);

run_scales = zeros(length(datapaths),1);
run_scales = [images_scales(1);images_scales(1);images_scales(1);...
              images_scales(2);g3g4_ave;g3g4_ave;g3g4_ave;g3g4_ave;g3g4_ave;images_scales(5)];

%% Load Prerun Images
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
% 
%             figure;
%             image(prerun_ave)
%             colorbar;
%             colormap(turbo(max(prerun_ave(:))));
%             grid on;
%             axis equal;
%             title(["Pre-run average from ",prerun_filefolder])
end

%% Fitting Prerun Images
for i = 1:length(datapaths)
    filepath = fullfile(folderpath,datapaths(i),prerun_foldernames(i),"AVERAGE.mat");
    load(filepath);

        figure;
        subplot(2,1,1);
        image(prerun_ave)
        colorbar;
        colormap(turbo(max(prerun_ave(:))));
        grid on;
        axis equal;
        title(["Pre-run average from ",datapaths(i)])
        subplot(2,1,2);
        plot(1:length(prerun_row_averages(i,:)),prerun_row_averages(i,:));

end

%% Load Run Images
        run_ims_count = 7;
        run_ims = zeros(200,2048,run_ims_count);
        run_row_averages = zeros(length(datapaths),2048,run_ims_count);
for i = 1:length(datapaths)
        run_filepath = fullfile(folderpath,datapaths(i),"run");
        have = ["tif"];
        nothave = ["rec"];
        [filteredfiles_run] = Folders_in_Folder(run_filepath,have,nothave);
        process_runs = filteredfiles_run(end-6:end);
    
%         figure;
        for k = 1:length(process_runs)
            fileNameCat = fullfile(run_filepath,process_runs(k));
            run_ims(:,:,k) = double(imread(fileNameCat));
            run_row_averages(i,:,k) = mean(run_ims(:,:,k),1);

%             subplot(2,1,1);
%             image(run_ims(:,:,k))
%             colorbar;
%             colormap(turbo(max(max(run_ims(:,:,k)))));
%             grid on;
%             axis equal;
%             title(['Image number ', num2str(k), process_runs(k)])
%             subplot(2,1,2);
%             plot(1:length(run_col_averages(:,:,k)),run_col_averages(:,:,k));
        end

end

%% Identify displacements

    %cross correlation
    for i = 1:length(datapaths) %runs
        for j = 1:size(run_row_averages,3) %images from each run
        prerun = prerun_row_averages(i,:);
        run = run_row_averages(i,:,j);
        corr = xcorr(prerun,run);

            figure;
            subplot(3,1,1);


        end
    end

    %Seperate fitting...

%% Velocity Calculation