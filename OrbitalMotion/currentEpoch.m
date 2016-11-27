function [ epochInJulianAge ] = currentEpoch( julianDate, timeStamp )
    % currentEpoch. Convert current Julian Date and provided time stamp to the ephemeresis time in Julian age.
    %   Provided time stamp converted to the UTC+3 time zone and then converted to the fractional day part.
    %   Provided Julian date with converted time stamp converted to the Julian age from 01 Jan 2000 12:00 AM (UTC).
    %
    %   [ epochInJulianAge ] = currentEpoch( julianDate, timeStamp )
    %
    %   INPUT
    %       julianDate    julian day;
    %       timeStamp     number of seconds from start of current julian day.
    %
    %   OUPUT
    %       epochInJulianAge   - date from 01 Jan 2000 12:00 AM (UTC) in Julian Age.
    %
    
    %   2451545.0 - 2000-Jan-01 12:00:00.00 (UTC)
    jdFromDate      = 2451545.0;
    julianAgeInDays = 36525;
    secondsInDay    = 86400;
    moscowUtcOffset = 10800;
    
    epochInJulianAge = (julianDate + (timeStamp - moscowUtcOffset) /secondsInDay - jdFromDate) / julianAgeInDays;
end
