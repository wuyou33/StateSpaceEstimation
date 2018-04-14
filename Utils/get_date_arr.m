function [ val ] = get_date_arr( startTime, timeFormat,  date)
    val = datevec(startTime, timeFormat);
    
    fields_date = fieldnames(date);
    
    for i = 1 : size(fields_date, 1)
        val(4-i) = date.( char( fields_date(i) ) );
    end
end
