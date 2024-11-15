function [varargout] = rsac(varargin)

    %RSAC    Read SAC binary files.
    %    RSAC('sacfile') reads in a SAC (seismic analysis code) binary
    %    format file into a 3-column vector.
    %    Column 1 contains time values.
    %    Column 2 contains amplitude values.
    %    Column 3 contains all SAC header information.
    %    Default byte order is big-endian.  M-file can be set to default
    %    little-endian byte order.
    %
    %    usage:  output = rsac('sacfile')
    %
    %    Examples:
    %
    %    KATH = rsac('KATH.R');
    %    plot(KATH(:,1),KATH(:,2))
    %
    %    [SQRL, AAK] = rsac('SQRL.R','AAK.R');
    %
    %    by Michael Thorne (4/2004)   mthorne@asu.edu

    for nrecs=1:nargin

        sacfile = varargin{nrecs};

        %----------------------------------------------------------------------
        % Default byte-order
        % endian  = 'big-endian' byte order (e.g., UNIX)
        %         = 'little-endian' byte order (e.g., LINUX)

        % First try using little-endian
        fid = fopen(sacfile,'r','ieee-le');

        % read in single precision real header variables:
        %----------------------------------------------------------------------
        for i=1:70
            hd(i,1) = fread(fid,1,'single');
        end

        % read in single precision integer header variables:
        %----------------------------------------------------------------------
        for i=71:105
            hd(i,1) = fread(fid,1,'int32');
        end

        % Check header version = 6 and issue warning
        %----------------------------------------------------------------------
        % If the header version is not NVHDR == 6 then the sacfile is likely of
        % the opposite byte order.  This will give h(77) some ridiculously
        % large number.  NVHDR can also be 4 or 5.  In this case it is an old
        % SAC file and rsac cannot read this file in.  To correct, read the SAC
        % file into the newest verson of SAC and w over.

        if hd(77)==4 || hd(77)==5
            message = strcat(['NVHDR = 4 or 5. File: "' sacfile '" may be from an old version of SAC.']); 
            error(message)
        elseif hd(77)~=6
    %         message = strcat(['Current rsac byte order: "' endian '". File: "' sacfile '" may be of opposite byte-order.']);
    %         error(message)

            % Try opening SAC file with big-endian byte order
            fid = fopen(sacfile,'r','ieee-be');

            % read in single precision real header variables:
            %------------------------------------------------------------------
            for i=1:70
                hd(i,1) = fread(fid,1,'single');
            end

            % read in single precision integer header variables:
            %------------------------------------------------------------------
            for i=71:105
                hd(i,1) = fread(fid,1,'int32');
            end

        end

        % read in logical header variables
        %----------------------------------------------------------------------
        for i=106:110
            hd(i,1) = fread(fid,1,'int32');
        end

        % read in character header variables
        %----------------------------------------------------------------------
        for i=111:302
            hd(i,1) = (fread(fid,1,'char'))';
        end

        % read in amplitudes
        %----------------------------------------------------------------------
        YARRAY = fread(fid,'single');

        if hd(106)==1
            XARRAY = single((linspace(hd(6),hd(7),hd(80)))'); 
        else
            error('LEVEN must = 1; SAC file not evenly spaced')
        end 

        hd(307,1) = 0;
        if length(YARRAY)~=length(XARRAY)
            nlength = max(length(YARRAY),length(XARRAY));
            YARRAY=zeros(nlength,1);
            XARRAY=zeros(nlength,1);
            hd(307,1) = 1;   %if 307 is 1, the dimension diesn't agree, and sac file is not read
        end
        % add header signature for testing files for SAC format
        %----------------------------------------------------------------------
        hd(303,1) = 77;
        hd(304,1) = 73;
        hd(305,1) = 75;
        hd(306,1) = 69;

        % arrange output files
        %----------------------------------------------------------------------
        OUTPUT(:,1) = XARRAY;
        OUTPUT(:,2) = YARRAY;
        OUTPUT(1:307,3) = hd(1:307);

        %pad xarray and yarray with NaN if smaller than header field
        if hd(80)<307
            OUTPUT((hd(80)+1):307,1) = NaN;
            OUTPUT((hd(80)+1):307,2) = NaN;
        end

        fclose(fid);

        varargout{nrecs} = OUTPUT;

    end
    
end