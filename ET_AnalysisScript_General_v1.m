%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Eye Tracking Analysis Script for ET-Data acquired with Vive Pro Eye
%
% Leon Kroczek, 2024
%
% CC-BY 4.0
%
% ReadMe
% 
% 1. Import
%
% This script requires a custom import function.
% Becasue data file may vary between experiemnts it is best to create the
% import function for each study

% This can be done "Right-clicking" on a raw file and selecting "Import
% Data" in MATLAB

% An import window opens: Check if columns are labeleld correctly and all
% columns and rows are selected.
% If so press on the arrow at "Import Selection" -> Generate Function
% Save the Fucntion as "importET.m" in the folder of this script 

% 2. Parameters

% The current procedure use a two-step method to identify fixations

% First, fixations are estimated based on a velocity criterion
% Everything slower than a pre-specified threshold is counted as fixation
% Parameter = threshold_vel
% This procedure uses the Savitzky-Golay Filter function savGol
%        W. H. Press and S. A. Teukolsky,
%        Savitzky-Golay Smoothing Filters,
%        Computers in Physics, 4 (1990), pp. 669-672.

% This procedure often produces a high number of fixations, because
% fixations can be split in several parts following when data is noisy 
% If one is only interested in dwell times (sum of all time on a spefici AOI) this should produce valid
% results. However, if on is interested in mean fixation duration or
% fication count a second step is recommended:
% Based on teh initial fixations, the algorithm checks whether subsequent
% fixations ahve the same target, are close in time, and close in space
% If two fixations have the same target and are below a threshold they are
% summarized as one

% 4. Marker
% The script requires that for all segments of interest both the
% onset and the offset are defined in the correct order.
% The script assumes that every marker appears once
% Other options can be easily inlcuded:
%   - Marker repetitions
%   - Marker onset but end with time criterion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Subjects and Folders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subjects to Process
subjects        = {'M24','M8','F41'};

% Path where raw data is located
path_eye        = '';

% How are ET files named after the individual Subject code?
% "XXXXXXX_F7_eye.txt" -> '_eye.txt'
file_suffix     = '_eye.txt';

% Get list of all relevant files in data folder
files_eye = dir(path_eye);

% Save Options: Path where CSV shoudl be saved
out_dir = '';

% Subject Logfile suffix: is appended after Subject name in CSV file: ["Subject + suffix.csv"]
output_suffix = '_ET';

%%% Settings for Fixation identification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters to identify fixations/saccades (intiial step)
threshold_vel   = 75; %in degree/second - everything above is counted as saccade

% Time criterion in ms
% How much time can be between fixations that they are counted as one 
time_c = 500;

% Space criterion in cm
% What distance between fixations is allowed so that they count as one
space_c = 30;

% Can the first 100 lines of datafile (~ 100*10ms = 1 sec) be ommited?
% set variable to 'yes'
% The first lines after initiating the ET logfile often contain
% inconsistencies

clean_initial = 'yes';

% Markers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Begin and End Markers for segments of interest
% Note this works only when there are unique begin and end markers for each segement
% and each marker appears only one time in the whole experiment

% Set Begin and End Markers
% Important: order must be the same

segment_begin_marker    = {'FirstViewOnComitee',            'ViewOnComitteeStartStress',        'BeginVorstellungsrunde',   'PhilipStellenSiesichvor_EBUR128_over'};
segment_end_marker      = {'PhilipVorstellungPeter_EBUR128','PhilipStellenSiesichvor_EBUR128',  'EndVorstellungsrunde',     'PhilipvielenDankdaswars_EBUR128'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Part 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for sub = 1:length(subjects) % Subject loop

    % Read Eye Tracking File (looks for File which matches Subject Code)
    eye_idx = find(cellfun(@(x) ~isempty(x), strfind({files_eye.name},[subjects{sub} file_suffix])));

    % Read data using custom import function
    data    = importET(fullfile(path_eye,files_eye(eye_idx).name));

    % Delete the first 100 data rows 
    % These lines often contain irregularities
    if strcmp(clean_initial,'yes')
        data(1:100,:) = [];
    end

    % recalc Unreal Timestamps
    % Reason: Sometimes Timestamps reset to zero
    % This line recalculates the timestamps to be continuously increasing  
    data.Timestamp(find(diff(data.Timestamp)<0)+1:end) = data.Timestamp(find(diff(data.Timestamp)<0)+1:end)+ data.Timestamp(find(diff(data.Timestamp)<0));

    % Diagnostics only
    % % Get time per frame
    % n_frames        = max(data.Frame)-data.Frame(1);
    % duration_all    = max(data.Timestamp)-min(data.Timestamp);
    % 
    % % Average Frame Duration of Unreal Presentatiobn
    % avg_sampleT     = duration_all/n_frames;


    % Calcuale time difference between samples (in Msec)
    data.deltaT = [0; diff(data.Timestamp)];

    % Export Counter
    l = 1;

    for segment = 1:length(segment_begin_marker)

        % Find Line indices corresponding to markers
        event_LineBegin    = find(strcmp(data.Marker,segment_begin_marker{segment}));
        event_LineEnd      = find(strcmp(data.Marker,segment_end_marker{segment}));

        % Extract data between Begin and End Markers
        seg = data(event_LineBegin:event_LineEnd,:);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Marker
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Indices at whixh marker is sent
        marker_change = find(~strcmp(seg.Marker,''));

        % Repeat Marker in every row until the next Marker is sent
        for m = 1:length(marker_change)
            if m < length(marker_change)
                seg.Marker(marker_change(m):marker_change(m+1)-1) = seg.Marker(marker_change(m));
            else
                seg.Marker(marker_change(m):end) = seg.Marker(marker_change(m));
            end
        end 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Fixations and saccades
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         % Calc angles between subsequent samples
         Vec1 = [(seg.LeftGazeX+seg.RightGazeX)/2 ...
                 (seg.LeftGazeY+seg.RightGazeY)/2 ...
                 (seg.LeftGazeZ+seg.RightGazeZ)/2];

         for j = 2:height(seg)
            seg.Angle_diff(j) = atan2d(norm(cross(Vec1(j,:),Vec1(j-1,:))),dot(Vec1(j,:),Vec1(j-1,:))); % Angle in degree
         end

         % Remove angles with invalid data
         bad_data = find(seg.LeftValidity~=31 | seg.RightValidity~=31);
         seg.Angle_diff(bad_data) = nan; 

        % Caluclate velocities and filter with savitzky golay filter
        seg.v       = seg.Angle_diff./(seg.deltaT/1000);
        seg.v2      = sgolayfilt(seg.v,2,5); %Requires Signal Processing Toolbox; Alternatively see savGol fucntion:  savGol(seg.v,2,2,2);

        % Labels as saccades, blinks and fixations
        % Write numeric Variable event_code = 1 for fixations and = 0 for
        % everything else
        
        %seg.event(1)                                                    = {'sac'};
        seg.event(find(seg.v2>threshold_vel))                           = {'sac'};
        seg.event(find(seg.v2<threshold_vel))                           = {'fix'};
        seg.event_code(find(seg.v2<threshold_vel))                      = 1;
        seg.event(find(seg.LeftGazeX==0 | seg.RightGazeX==0))           = {'blk'};
        seg.event_code(find(seg.LeftGazeX==0 | seg.RightGazeX==0))      = 0;


        % Delete ultra-short fixation/saccades
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Find changes between sacades and fixations
        tmp_diff        = [0;diff(seg.event_code)];
        change_idx      = [1;find(tmp_diff ~= 0)];

        tmp_event = [];
        % Calculate length of single fixations/saccades
        for c = 1:length(change_idx)-1
            tmp_event(c,1) = change_idx(c);
            tmp_event(c,2) = change_idx(c+1)-1;
        end

        tmp_event(:,3) = tmp_event(:,2) - tmp_event(:,1)+1;

        % Delete fixation/saccades <= two samples

        to_delete = [];
        for c = 1:size(tmp_event,1)

            if tmp_event(c,3) <= 4
                to_delete = [to_delete,tmp_event(c,1):tmp_event(c,2)];
            end 

        end

        seg(to_delete,:) = [];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get Fixations: Initial Computation based on Velocity
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Find changes between sacades and fixations (new)
        tmp_diff        = [0;diff(seg.event_code)];
        change_idx      = find(tmp_diff ~= 0);

        % Enumerate Fixations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fix_num = 1;

        % no fixations found
        if ~isempty(find(strcmp(seg.event,'fix'))) % at least one fixation has to be found

            if ~isempty(change_idx) % if more than 1 fixation was found

                for k = 1:length(change_idx)
                    if tmp_diff(change_idx(k)) == -1 && k == 1 % first change: fix -> sac
                        seg.fix_num(1:change_idx(k)-1) = fix_num;
                        fix_num = fix_num +1;
                    elseif tmp_diff(change_idx(k)) == -1 && k ~= 1 % intermediate change fix -> sac
                         seg.fix_num(change_idx(k-1):change_idx(k)-1) = fix_num;
                         fix_num = fix_num +1;
                    elseif tmp_diff(change_idx(k)) == 1 && k == length(change_idx) % last change: sac -> fix
                         seg.fix_num(change_idx(k):end) = fix_num;
                         fix_num = fix_num +1;
                    end
                end

            else % if only one big fixations was found
                seg.fix_num = repmat(fix_num,height(seg),1);
            end
        end


        % Calc Duration,Latency, Impact Component and Impact Coordiantes (mean) of each fixation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Total number of fixations
        FIX_NUM = setdiff(unique(seg.fix_num),[0]);

        % Loop through all fixations
        for k = 1:length(FIX_NUM)

            
            tmp_on  = find(seg.fix_num == FIX_NUM(k),1,'first');    % Get index of start of fixatiom
            tmp_off = find(seg.fix_num == FIX_NUM(k),1,'last');     % Get index of end of fixation

            % Calc Duration and Latency
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            seg.fix_dur(find(seg.fix_num == FIX_NUM(k))) =  sum(seg.deltaT(tmp_on:tmp_off));
            seg.fix_lat(find(seg.fix_num == FIX_NUM(k))) =  seg.Timestamp(tmp_on)  - seg.Timestamp(1);

            % Get Name of Impact Component
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Find all names of imaopct components that are associated with
            % one fixation
           [s,~,j] = unique(setdiff(seg.Hit(tmp_on:tmp_off),''));
            %s{mode(j)};

            % Use the most frequent one as real impact component 
            if ~isempty(s)
                seg.fix_impact(find(seg.fix_num == FIX_NUM(k))) = {s{mode(j)}};
            else
                seg.fix_impact(find(seg.fix_num == FIX_NUM(k))) = {''};
            end

            % Coordinates of Impact Point
            % Caution: This uses the geometric mean of all coordinates
            % within a fixation based on the Impact Point/Hit Point
            % In 3D the coordinates might therefore be shifted away from
            % the participant position 
            % Alternatvely, one might use the point closest to the
            % participant (not implemented here)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            seg.fix_x(find(seg.fix_num == FIX_NUM(k))) = mean(seg.HitX(find(seg.fix_num == FIX_NUM(k))),'omitnan');
            seg.fix_y(find(seg.fix_num == FIX_NUM(k))) = mean(seg.HitY(find(seg.fix_num == FIX_NUM(k))),'omitnan');
            seg.fix_z(find(seg.fix_num == FIX_NUM(k))) = mean(seg.HitZ(find(seg.fix_num == FIX_NUM(k))),'omitnan');

            % Coordinates of VP Position
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            seg.cam_x(find(seg.fix_num == FIX_NUM(k))) = mean(seg.CameraX(find(seg.fix_num == FIX_NUM(k))),'omitnan');
            seg.cam_y(find(seg.fix_num == FIX_NUM(k))) = mean(seg.CameraY(find(seg.fix_num == FIX_NUM(k))),'omitnan');
            seg.cam_z(find(seg.fix_num == FIX_NUM(k))) = mean(seg.CameraZ(find(seg.fix_num == FIX_NUM(k))),'omitnan');

            % angle deviation from target
            %seg1.data.fix_angle(find(seg1.data.fix_num == FIX_NUM(k))) =  mean(seg1.data.angle_target(find(seg1.data.fix_num == FIX_NUM(k))),'omitnan');

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get Fixations2: Second Computation based on Time and Dispersion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check if consecutive fixations are similar in name, time or space
        % Reprocess fixations and calculate new values for updated fixations
        % Combines fixations that are similar
        % parameters set under settings
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % More than one fixation?
        if length(FIX_NUM) > 1

            % helper structure where info about fixation pairs is stored
            flag = []; % Stores pairs not single fixations

            % Loop through consecutive fixations
            for f = 1:length(FIX_NUM)-1

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Name of Impact Component 
                % Do they have the same target (name)?
                % Check: name next == name_before?
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                name_before = string(unique(seg.fix_impact(seg.fix_num == FIX_NUM(f))));
                name_next   = string(unique(seg.fix_impact(seg.fix_num == FIX_NUM(f+1))));
                
                if  name_before == name_next 
                    flag(f,1) = 1;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Time
                % How much time lies between end of first and start of
                % second fixation in msec
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % End first = Latency + Duration
                end1        = unique(seg.fix_lat(seg.fix_num == FIX_NUM(f))) - unique(seg.fix_dur(seg.fix_num == FIX_NUM(f)));
                % Start Second = Latency
                start2      = unique(seg.fix_lat(seg.fix_num == FIX_NUM(f+1)));

                flag(f,2)   = start2 - end1;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Space
                % What is the spatial distance between second and first fixation (in cm)?
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                coord_fixation_before   = [unique(seg.fix_x(seg.fix_num == FIX_NUM(f))) unique(seg.fix_y(seg.fix_num == FIX_NUM(f))) unique(seg.fix_z(seg.fix_num == FIX_NUM(f)))];
                coord_fixation_next     = [unique(seg.fix_x(seg.fix_num == FIX_NUM(f+1))) unique(seg.fix_y(seg.fix_num == FIX_NUM(f+1))) unique(seg.fix_z(seg.fix_num == FIX_NUM(f+1)))];
                
                
                flag(f,3) = norm(coord_fixation_next- coord_fixation_before);

            end 

            % Define fixation that are obviously connected (= same)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Criteria - see above (currently checks all three)
                
            same(1,1)           = 1; % Stores fixtions not fixation pairs!
            new_fixation_index  = 1;
            
            for f = 1:size(flag,1)
                
                if flag(f,1) == 1 && flag(f,2) < time_c && flag(f,3) < space_c
                    same(f+1,1) = new_fixation_index;
                else
                    new_fixation_index = new_fixation_index+1;
                    same(f+1,1) = new_fixation_index;
                end
            end
            
           

            % Calculate new parameters for updated fixations
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for j = 1:max(same)

                % Recalc fixation number, latency (= earliest begin),
                % duration ( = sum of all fixations that are summarized)
                seg.fix_num2(find(ismember(seg.fix_num,FIX_NUM(same == j)))) = j;
                seg.fix_lat2(find(ismember(seg.fix_num,FIX_NUM(same == j)))) = min(seg.fix_lat(find(ismember(seg.fix_num,FIX_NUM(same == j)))));
                seg.fix_dur2(find(ismember(seg.fix_num,FIX_NUM(same == j)))) = sum(unique(seg.fix_dur(find(ismember(seg.fix_num,FIX_NUM(same == j))))));

                % 
                seg.fix_x2(find(ismember(seg.fix_num,FIX_NUM(same == j)))) = mean(unique(seg.fix_x(find(ismember(seg.fix_num,FIX_NUM(same == j))))));
                seg.fix_y2(find(ismember(seg.fix_num,FIX_NUM(same == j)))) = mean(unique(seg.fix_y(find(ismember(seg.fix_num,FIX_NUM(same == j))))));
                seg.fix_z2(find(ismember(seg.fix_num,FIX_NUM(same == j)))) = mean(unique(seg.fix_z(find(ismember(seg.fix_num,FIX_NUM(same == j))))));

                seg.cam_x2(find(ismember(seg.fix_num,FIX_NUM(same == j)))) = mean(unique(seg.cam_x(find(ismember(seg.fix_num,FIX_NUM(same == j))))));
                seg.cam_y2(find(ismember(seg.fix_num,FIX_NUM(same == j)))) = mean(unique(seg.cam_y(find(ismember(seg.fix_num,FIX_NUM(same == j))))));
                seg.cam_z2(find(ismember(seg.fix_num,FIX_NUM(same == j)))) = mean(unique(seg.cam_z(find(ismember(seg.fix_num,FIX_NUM(same == j))))));
                
            end

        % If there is only one fixation in the segment just copy the data
        else

            seg.fix_num2 = seg.fix_num;
            seg.fix_lat2 = seg.fix_lat;
            seg.fix_dur2 = seg.fix_num;

            seg.fix_x2 = seg.fix_x;
            seg.fix_y2 = seg.fix_y;
            seg.fix_z2 = seg.fix_z;

            seg.cam_x2 = seg.cam_x;
            seg.cam_y2 = seg.cam_y;
            seg.cam_z2 = seg.cam_z;

        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Write Data in Export structure
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if ~isempty(find(strcmp(seg.event,'fix'))) 

            for s = 1:max(seg.fix_num2)

                exp_new{l,1} = subjects{sub};      % Subject
                exp_new{l,2} = segment;            % Trial Num   
                exp_new{l,3} = seg.Marker(find(seg.fix_num2 == s,1,'first'));   % Writes marker code to correspondign fixation (only relevant when markers change within a segment)
                exp_new{l,4} = s;                   % Fixation Num

                exp_new{l,5} = seg.fix_x(find(seg.fix_num2 == s,1,'first'));       % Fixation Coord X
                exp_new{l,6} = seg.fix_y(find(seg.fix_num2 == s,1,'first'));       % Fixation Coord Y
                exp_new{l,7} = seg.fix_z(find(seg.fix_num2 == s,1,'first'));       % Fixation Coord Z

                exp_new{l,8} = seg.cam_x(find(seg.fix_num2 == s,1,'first'));       % Camera Coord X
                exp_new{l,9} = seg.cam_y(find(seg.fix_num2 == s,1,'first'));       % Camera Coord Y
                exp_new{l,10} = seg.cam_z(find(seg.fix_num2 == s,1,'first'));      % Camera Coord Z

                exp_new{l,11} = seg.fix_dur2(find(seg.fix_num2 == s,1,'first'));      % Fixation Duration
                exp_new{l,12} = seg.fix_lat2(find(seg.fix_num2 == s,1,'first'));      % Fixation Latency

                exp_new{l,13} = cell2mat(seg.fix_impact(find(seg.fix_num2 == s,1,'first')));      % Fixation Impact Name

                exp_new{l,14} = cell2mat(seg.event(find(seg.fix_num2 == s,1,'first'))); % Type
                l = l+1;


            end

        else % If no fixation was found write empty columns
             
            exp_new{l,1} = subjects{sub};      % Subject
            exp_new{l,2} = segment;                           % Trial Num   
            exp_new{l,3} = seg.Marker(1);               % Marker Code 
            exp_new{l,4} =[];               % Sound

            exp_new{l,5} = 0;                           % Fixation Num
            exp_new{l,6} = [];
            exp_new{l,7} = [];       % Fixation Coord Y
            exp_new{l,8} = [];
            exp_new{l,9} = [];
            exp_new{l,10} = [];
            exp_new{l,11} = [];

            exp_new{l,12} = [];
            exp_new{l,13} = [];

            exp_new{l,14} = [];

            l = l+1;

        end  

        clear seg same flag

    end

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save Subject Specific Logfile
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    name = [subjects{sub} output_suffix];
    fid = fopen(fullfile(out_dir,[name '.csv']),'w');                                         % Open the file

    % Column Names
    NAME = {'Subject','Trial','Marker','Fix_num','Fix_Coord_X','Fix_Coord_Y','Fix_Coord_Z','Camera_Pos_X','Camera_Pos_Y','Camera_Pos_Z','Fix_dur','Fix_lat','Fix_impact_name','Type'};
    fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',NAME{1:end});
    
    % Write Data: Line by Line
    for line = 1:size(exp_new,1)
         fprintf(fid,'%s,%g,%s,%g,%g,%g,%g,%g,%g,%g,%g,%g,%s,%s\n',exp_new{line,1:end});                                  % Write line per line
    end
    fclose(fid);    

    
    clear data exp_new   
    
end % Subject End





        


        
        
    
    
   
    
    











