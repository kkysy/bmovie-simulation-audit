function checkpoint = run_acc00_powercurve_profile_v2(root, gridIndex, seedIndex)
% Authorized v2 profiling wrapper; v1 checkpoints are never overwritten.
arguments
    root (1,1) string
    gridIndex (1,1) double {mustBeMember(gridIndex,[1 6 7])}
    seedIndex (1,1) double {mustBeInteger,mustBePositive} = 1
end
checkpoint=run_acc00_powercurve_world(root,gridIndex,seedIndex,"v2");
fprintf("P0-2 v2 profiling %s seed %d complete (%.3f s; peak %.0f bytes).\n",checkpoint.grid_id,checkpoint.seed,checkpoint.timing.elapsed_s,checkpoint.timing.peak_memory_bytes);
end
