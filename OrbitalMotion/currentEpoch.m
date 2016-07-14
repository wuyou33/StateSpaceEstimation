function currentEpoch = currentEpoch( julianDate, previousEpoch )
% from 2000 year, 01 Jan 12:00 AM (UTC) to time of ephemeresis in Julian% age.
    currentEpoch = (julianDate + (previousEpoch - 10800) /86400 - 2451545.0) / 36525;
end
