%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Aug 2015.
%
% Track individuals for John Hedley's visualization tool
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function f_generate_track_files(META, REEF, CORAL, ALGAL, RECORD, colony_list, environ_list)

disturb_list = [0:META.nb_time_steps ; iseven(0:META.nb_time_steps)]' ;
disturb_list(2:META.nb_time_steps+1,3) = RECORD.bleaching_events(1,1:META.nb_time_steps)' ;
disturb_list(2:META.nb_time_steps+1,4) = RECORD.hurricane_events(1,1:META.nb_time_steps)' ;

[data, result]= readtext('track_info_template','[,\t]');
data(10,1) = {META.grid_x_count * META.grid_y_count};
data(13,1) = {META. cell_area_cm2} ;
data(16,1) = {META.nb_time_steps};
data(19,1) = {META.nb_coral_types} ;
data(22,1) = {META.max_colonies} ;
data(25,1) = {META.nb_algal_types} ;
data(28,1) = {100*REEF.nongrazable_substratum};
data(31,1:META.nb_coral_types)= {100*REEF.initial_coral_cover};

dlmcell('track_info.txt',data,' ')
csvwrite('track_colony.csv', colony_list(2:end,:)) ; % excludes first row which was a primer
csvwrite('track_environ.csv', environ_list(2:end,:)) ; % excludes first row which was a primer
csvwrite('track_disturb.csv', disturb_list);

zip('track_outputs.zip',{'track_colony.csv','track_environ.csv','track_disturb.csv','track_info.txt'})

delete('track_colony.csv','track_environ.csv','track_disturb.csv','track_info.txt')