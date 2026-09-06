function finalPath = acc01_write_session_checkpoint(finalPath, checkpoint)
%ACC01_WRITE_SESSION_CHECKPOINT Publish one checkpoint by tmp.mat rename.

finalPath = string(finalPath);
[folder, name, ext] = fileparts(finalPath);
if ~isfolder(folder), mkdir(folder); end
tmpPath = fullfile(folder, name + ".tmp" + ext);
if isfile(finalPath)
    error("Bmovie:ACC01:CheckpointExists", "Refusing to overwrite complete checkpoint: %s", finalPath);
end
if isfile(tmpPath), delete(tmpPath); end
save(tmpPath, "checkpoint", "-v7.3");
movefile(tmpPath, finalPath, "f");
end
