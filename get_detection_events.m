% 'detection_case' leads to the desired detection criterion
% 'whe_single' decides the intention of picking out independent detected
% events. For example, from the detection criterion, five consecutive
% events are picked out, if 'whe_single == 1', then only one event would be
% selected from those five consecutive events. 
function [pick_events] = get_detection_events(all_events_array,k_scale,detection_case,whe_single)
    
    switch detection_case
        case 1
            pick_events = all_events_array.^2 > (k_scale * var(all_events_array));
        case 2
            pick_events = all_events_array    > (k_scale * mean(all_events_array));
        case 3
            pick_events = all_events_array    < (k_scale * mean(all_events_array));
        case 4
            pick_events = all_events_array    > (k_scale * max(all_events_array));
        otherwise
            disp('This detection criterion has not been coded.')
    end
    

    if  whe_single == 1
    repeat_left  = zeros(length(pick_events)-1,1);
    repeat_right = zeros(length(pick_events)-1,1);
    % left
    k_array1 = 1;
    if pick_events(1)==1
        repeat_left(k_array1) = 1;
        k_array1 = k_array1 + 1;
    end
    for k_array =1: length(pick_events)-1
        if pick_events(k_array)==0 && pick_events(k_array+1)==1
            repeat_left(k_array1) = k_array;
            k_array1 = k_array1 + 1;
        end
    end
    % right
    k_array1 = 1;
    for k_array =2: length(pick_events)
        if pick_events(k_array-1)==1 && pick_events(k_array)==0
            repeat_right(k_array1) = k_array;
            k_array1 = k_array1 + 1;
        end
    end
    if pick_events(end)==1
        repeat_right(k_array1) = length(pick_events);
        k_array1 = k_array1 + 1;
    end

    repeat_left(repeat_left==0)   = [];
    repeat_right(repeat_right==0) = [];

    repeat_index = round((repeat_left + repeat_right) ./ 2);
    pick_events = 0.* pick_events;
    pick_events(repeat_index) = 1;
    end

end