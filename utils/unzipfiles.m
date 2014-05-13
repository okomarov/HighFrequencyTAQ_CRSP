function outfiles = unzipfiles(zipFilename, files, outputDirectory)
%   UNZIPFILES Extract the contents of a zip file.
%
%   UNZIPFILES(ZIPFILE)
%       extracts the contents of a zip file into the current
%      directory.
%
%   UNZIPFILES(ZIPFILE,FILES)
%       extracts the files specified in FILES, a string or a
%       a cell array of strings (one or more), from the zipped file
%       into the current directory.
%
%   UNZIPFILES(ZIPFILE,FILES,OUTPUTDIR)
%       extracts the files specified in FILES, a string or a
%       a cell array of strings, from the zipped file
%       into the directory OUTPUTDIR
%
%   Examples:
%           unzipfiles('junk.zip')
%               will extract the entire contents of junk into the current
%              current directory.
%          unzipfiles('junk.zip','temp1.mat')
%               will extract "temp1.mat" into the current
%              current directory.
%          unzipfiles('junk.zip','temp1.mat','c:\temp')
%               will extract "temp1.mat" into c:\temp
%          unzipfiles('junk.zip',[],'c:\temp')
%               will extract the entire contents of junk into c:\temp
%
%   See also ZIP, UNZIP.


%   This file is an extention of  The MathWorks, Inc's UNZIP function
%
%   Author :    Tal Pasi
%   Date:     15-Oct-2004

% Author: Oleg Komarov (o.komarov11@bcsprime.com)
% Tested on R2014a Win7 64bit
% 13 May 2014 - Refactored and corrected name matching 

import java.io.*;
import java.util.zip.ZipFile;
import com.mathworks.mlwidgets.io.InterruptibleStreamCopier;

% Argument parsing
error(nargchk(1,3,nargin));
if nargin < 3 
    outputDirectory = pwd;
elseif ~exist(outputDirectory,'dir')
    error('Directory "%s" does not exist.',outputDirectory)
end

if nargin < 2 || isempty(files)
    unzip(zipFilename, outputDirectory);
    return
end

outfiles = {};

% Open the Zip file.
if ~exist(zipFilename,'file')
    error('File "%s" does not exist.',zipFilename);
end
try
    zipFile = ZipFile(zipFilename);
catch
    error('Error opening zip file "%s".',zipFilename);
end

% This InterruptibleStreamCopier is unsupported and may change without notice.
interruptibleStreamCopier = InterruptibleStreamCopier.getInterruptibleStreamCopier;

enumeration = zipFile.entries;
while enumeration.hasMoreElements
    zipEntry   = enumeration.nextElement;
    outputName = fullfile(outputDirectory,files(strcmpi(zipEntry,files)));
    if ~isempty(outputName)
        file = java.io.File(outputName);
        parent = File(file.getParent);
        parent.mkdirs;
        try
            fileOutputStream = java.io.FileOutputStream(file);
        catch
            warning('Could not create "%s".',outputName);
        end
        % Extract entry to output stream.
        inputStream = zipFile.getInputStream(zipEntry);
        interruptibleStreamCopier.copyStream(inputStream,fileOutputStream);
        if nargout == 1
            outfiles = [outfiles, char(zipEntry)]; %#ok<AGROW>
        end
        % Close streams.
        fileOutputStream.close;
        inputStream.close;
    end
end
% Close zip.
zipFile.close;
end