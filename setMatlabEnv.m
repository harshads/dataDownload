function setMatlabEnv

    %%% Set the Matlab environment variable
    pth = getenv('PATH');
    if isempty(strfind(pth, '/Volumes/ToolsMac/bin'))
        setenv('PATH', [pth ':/Volumes/ToolsMac/bin/']);
    end

    pth = getenv('PATH');
    if isempty(strfind(pth, '/Volumes/ToolsMac/fsl-5.0.5/bin/'))
        setenv('PATH', [pth ':/Volumes/ToolsMac/fsl-5.0.5/bin/']);
    end

    setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
