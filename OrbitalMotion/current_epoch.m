function [ epoch_in_julian_age ] = current_epoch( julian_date, time_stamp )
    % current_epoch. Convert current Julian Date and provided time stamp to the ephemeresis time in Julian age.
    %   Provided time stamp converted to the UTC+3 time zone and then converted to the fractional day part.
    %   Provided Julian date with converted time stamp converted to the Julian age from 01 Jan 2000 12:00 AM (UTC).
    %
    %   [ epoch_in_julian_age ] = current_epoch( julian_date, time_stamp )
    %
    %   INPUT
    %       julian_date    julian day;
    %       time_stamp     number of seconds from start of current julian day.
    %
    %   OUPUT
    %       epoch_in_julian_age   - date from 01 Jan 2000 12:00 AM (UTC) in Julian Age.
    %
    
    % 2451545.0 - 2000-Jan-01 12:00:00.00 (UTC)
    jdFromDate      = 2451545.0;
    julianAgeInDays = 36525;
    secondsInDay    = 86400;
    moscowUtcOffset = 10800;
    
    epoch_in_julian_age = (julian_date + (time_stamp - moscowUtcOffset) / secondsInDay - jdFromDate) / julianAgeInDays;
end
