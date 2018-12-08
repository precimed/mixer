if isunix
    bgmg_shared_library = 'lib/libbgmg.so';
elseif ispc
    bgmg_shared_library = 'bin/bgmg.dll';
else
    disp('Platform not supported')
end
bgmg_shared_library_header = 'bgmg_matlab.h';
logfile = 'test_mixer_plugin.bgmglib.log';
BGMG_cpp.unload(); 
BGMG_cpp.load(bgmg_shared_library, bgmg_shared_library_header);
BGMG_cpp.init_log(logfile);
BGMG_cpp.log('mixer plugin setup success');
