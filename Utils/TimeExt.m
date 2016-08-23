classdef TimeExt < handle
    %TIMEEXT Summary of this class goes here
    
    properties (Access = private)
        startTime;                      % 'HH:MM:SS.FFF'
        duration;                       % 'HH:MM:SS.FFF'
        sampleTime;                     % [sec]
        date;                           % day-month-year
        refreshSunMoonInfluenceTime;    % [sec]
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
            dateArr = this.GetDateArr();
            val = this.tFun( dateArr(4:6) );
        end
        
        function val = get.EndSecond(this)
            deltaDateTime = datevec(this.duration, this.timeFormat);
            dayAdd = deltaDateTime(3) - 1;
            val = dayAdd * 3600 * 24 + this.tFun(deltaDateTime(4:6));
        end
        
        function val = get.JD(this)
            dateArr = this.GetDateArr();
            % n4 - slot number by 4 years (leap between each) that contains the current date
            n4 = fix( (dateArr(1)-this.referenceYear) / 4 ) + 1;
            
            % nt - number of days within a 4-year leap period
            nt = datenum( dateArr(1:3) ) - datenum([this.referenceYear + (n4-1)*4 1 1])+1;
            
            val = 1461*(n4 - 1) + nt + 2450082.5;
        end
        
        function val = get.RefreshSunMoonInfluenceTime(this)
            val = this.refreshSunMoonInfluenceTime;
        end
    end
    
    methods (Access = private)
        function val = GetDateArr(this)
            val = datevec(this.startTime, this.timeFormat);
            
            fieldsDate = fieldnames(this.date);
            for i = 1:size(fieldsDate, 1)
                val(4-i)=this.date.( char( fieldsDate(i) ) );
            end            
        end
        
        function val = tFun (~, x)
            val = x(1)*60*60 + x(2)*60 + x(3);
        end
    end
end
