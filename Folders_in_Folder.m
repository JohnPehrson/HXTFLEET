function [filteredfiles] = Folders_in_Folder(folderpath,matching,unmatching)
%This function reports a list of folders in a specified folder as a string
%vector. The user can also input a matching string for the folder that
%folder names must match.
%--------------- Inputs ---------------%
% folderpath: Complete filepath to the folder
% matching: Vector of strings that I want each filtered file to use
% unmatching: Vector of strings I don't want to appear in any filtered file
%--------------- Output ---------------%
% filteredfiles: All files or folders in the original folder that passed
    % filtering
    
%% Function

    %check to see if the folder exists
    if exist(folderpath, 'file') == 7 %data exists
    
            %getting all contents in the folder
            a=dir(folderpath);
            if size(a,1)<3
                disp('ERROR: Specified folder has no contents');
            else
                a = a(3:end,:);
                foldernames = strings(size(a,1),1);
                for i = 1:length(foldernames)
                    foldernames(i) = a(i).name;
                end
            end
    
            %Compare the folder contents against the desired matching or
            %non-matching strings
            tf_match = zeros(length(foldernames),length(matching));
            for i = 1:length(matching)
                tf_match(:,i) = contains(foldernames,matching(i));
            end
            tf_unmatch = zeros(length(foldernames),length(unmatching));
            if unmatching(1)~=""
                for i = 1:length(unmatching)
                    tf_unmatch(:,i) = not(contains(foldernames,unmatching(i)));
                end
            else
                tf_unmatch = ones(length(foldernames),length(unmatching));
            end
            all_use = (sum(tf_match,2)+sum(tf_unmatch,2))==(length(matching)+length(unmatching));
            filteredfiles = foldernames(all_use);
    else
    disp("Specified folder does not exist");
    end

end