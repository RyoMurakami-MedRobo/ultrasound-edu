function setup_must(varargin)
%SETUP_MUST  Download MUST (Matlab UltraSound Toolbox) and add it to the path.
%
%   SETUP_MUST()
%   SETUP_MUST('InstallDir', '~/Documents/MATLAB/MUST')
%   SETUP_MUST('Force', true)
%
%   This tool works without MUST (SIM_ENGINE falls back to an analytic mock
%   backend), but the MUST backend gives physically accurate RF data and is
%   required to reproduce the validation results in validation/ and in the
%   accompanying paper.
%
%   MUST is distributed by its authors under the GNU LGPLv3 license from
%   https://www.biomecardio.com/MUST/ and is NOT bundled with this
%   repository (a different license, and a large binary payload that does
%   not belong in a small teaching-tool repo). This function downloads the
%   official zip and adds it to the current MATLAB path; run it once per
%   session (or add the resulting addpath line to your startup.m to make it
%   permanent).
%
%   Name-value options:
%     'InstallDir'  Where to unzip MUST. Default: ~/Documents/MATLAB/MUST
%     'Force'       Re-download even if InstallDir already contains MUST
%                   functions. Default: false.
%     'AddToPath'   Add InstallDir to the MATLAB path after installing.
%                   Default: true. Uses addpath (session-only); call
%                   `savepath` yourself if you want it to persist.
%
%   See also SIM_ENGINE, VALIDATION/VALIDATE_MUST_VS_MOCK.

p = inputParser();
p.addParameter('InstallDir', fullfile(char(java.lang.System.getProperty('user.home')), ...
                                       'Documents', 'MATLAB', 'MUST'));
p.addParameter('Force', false);
p.addParameter('AddToPath', true);
p.parse(varargin{:});
installDir = p.Results.InstallDir;
force      = p.Results.Force;
addToPath  = p.Results.AddToPath;

alreadyInstalled = exist(fullfile(installDir, 'simus.m'), 'file') == 2;
if alreadyInstalled && ~force
    fprintf('MUST already found at %s (pass ''Force'', true to re-download).\n', installDir);
else
    url = 'https://www.biomecardio.com/MUST/functions/MUST.zip';
    fprintf('Downloading MUST from %s ...\n', url);
    tmpZip = [tempname() '.zip'];
    try
        websave(tmpZip, url);
    catch ME
        error('setup_must:download', ...
            ['Could not download MUST automatically (%s).\n' ...
             'Download it by hand from https://www.biomecardio.com/MUST/ and\n' ...
             'unzip it into %s, or pass a different InstallDir.'], ME.message, installDir);
    end

    if ~exist(installDir, 'dir'), mkdir(installDir); end
    fprintf('Extracting to %s ...\n', installDir);
    unzip(tmpZip, installDir);
    delete(tmpZip);
end

if addToPath
    addpath(genpath(installDir));
    fprintf('Added %s to the MATLAB path for this session.\n', installDir);
end

ok = ~isempty(which('simus')) && ~isempty(which('txdelay')) && ~isempty(which('das'));
if ok
    fprintf('MUST is ready: simus() -> %s\n', which('simus'));
else
    warning('setup_must:notFound', ...
        'MUST was installed but simus()/txdelay()/das() are not all on the path. Check %s.', installDir);
end
end
