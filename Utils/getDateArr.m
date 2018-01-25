function [ val ] = getDateArr( startTime, timeFormat,  date)
    val = datevec(startTime, timeFormat);
    
    fieldsDate = fieldnames(date);
    
    for i = 1 : size(fieldsDate, 1)
        val(4-i) = date.( char( fieldsDate(i) ) );
    end
end
