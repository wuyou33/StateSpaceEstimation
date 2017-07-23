classdef TimeExt < handle
    % TimeExt. Class which alllow to easy manipulate with date time variables.
    
    properties (Access = private)
        startTime;                      % 'HH:MM:SS.FFF'
        duration;                       % 'HH:MM:SS.FFF'
        sampleTime;                     % [sec]
        date;                           % day-month-year
        refreshSunMoonInfluenceTime;    % [sec]
        startSecond;                    % [sec]
        endSecond;                      % [sec]
    end
    
    properties (Access = private, Constant)
        dateFormat = 'dd.mm.yyyy';
        timeFormat = 'HH:MM:SS.FFF';
        referenceYear = 1996;
    end
    
    properties (Dependent, Access = public)
        StartTime;
        Duration;
        SampleTime;
        Date;
        SimulationNumber;
        TotalSeconds;
        JD;
        Time;
        RelTime;
        DigitTime;
        StartSecond;
        EndSecond;
        RefreshSunMoonInfluenceTime;
    end
    
    methods (Access = public)
        function obj = TimeExt(startTime, duration, sampleTime, date, refreshSunMoonInfluenceTime)
            %   Constructor
            %       startTime   - "hh.mm.ss.sss" (string);
            %       duration    - "hh.mm.ss.sss" (string);
            %       sampleTime  - sample time (number);
            %       date        - current date "yyy.mm.dd" (string);
            %       refreshSunMoonInfluenceTime - time in second when information about Moon & Sun influence should be refreshed (number).
            
            obj.startTime                   = startTime;
            obj.duration                    = duration;
            obj.sampleTime                  = sampleTime;
            obj.date                        = date;
            obj.refreshSunMoonInfluenceTime = refreshSunMoonInfluenceTime;
            
            dateArr = getDateArr(startTime, TimeExt.timeFormat,  date);
            obj.startSecond = TimeExt.tFun( dateArr(4:6) );
            
            deltaDateTime = datevec(duration, TimeExt.timeFormat);
            dayAdd = deltaDateTime(3) - 1;
            obj.endSecond = dayAdd * 3600 * 24 + TimeExt.tFun(deltaDateTime(4:6));
        end
        
        function sample = evalSampleFromTime(this, t)
            sample = find(abs(this.Time - t) < (this.sampleTime / 10000) );
        end
        
        function range = getSimulationRange(this, timeScale)
            %   Calculate sample array for whole simulation interval specifed in timeScale
            narginchk(2, 2);
            
            if ~isa(timeScale, 'TimeExt'); error(' [ TimeExt:getRange ] timeScale should be instance of the TimeExt'); end
            
            batchSize = this.SampleTime \ timeScale.SampleTime;
            range = 1 : batchSize : this.SimulationNumber;
        end
        
        function unixEpoch = timeStampToUnixEpoch(this, timeStampSec)
            startOfUnixEpcoh = datenum('01-Jan-1970');
            currentTime = datenum(this.GetDateArr()) + TimeExt.secToDate(timeStampSec);
            unixEpoch = int32( floor((currentTime - startOfUnixEpcoh) * 86400) );
        end
    end
    
    methods
        function val = get.StartTime(this)
            val = this.startTime;
        end
        
        function val = get.Duration(this)
            val = this.duration;
        end
        
        function val = get.SampleTime(this)
            val = this.sampleTime;
        end
        
        function val = get.Date(this)
            val = this.date;
        end
        
        function val = get.SimulationNumber(this)
            val = length(this.Time);
        end
        
        function val = get.TotalSeconds(this)
            val = this.EndSecond - this.StartSecond;
        end
        
        function val = get.Time(this)
            val = this.StartSecond : this.sampleTime : this.EndSecond;
        end
        
        function val = get.RelTime(this)
            val = this.Time;
        end
        
        function val = get.DigitTime(this)
            val = 1 : this.SimulationNumber;
        end
        
        function val = get.StartSecond(this)
            val = this.startSecond;
        end
        
        function val = get.EndSecond(this)
            val = this.endSecond;
        end
        
        function val = get.JD(this)
            dateArr = this.GetDateArr();
            val = juliandate(dateArr);
        end
        
        function val = get.RefreshSunMoonInfluenceTime(this)
            val = this.refreshSunMoonInfluenceTime;
        end
    end
    
    methods (Access = private)
        function val = GetDateArr(this)
            val = getDateArr(this.startTime, this.timeFormat,  this.date);
        end
    end
    
    methods (Access = private, Static)
        function dayNum = secToDate(secNum)
            % 1 day = 24 (hour) * 60 (minute) * 60 (sec)
            dayNum = secNum / 24 / 60 / 60;
        end
        
        function val = tFun (x)
            val = x(1)*60*60 + x(2)*60 + x(3);
        end
    end
    
end
