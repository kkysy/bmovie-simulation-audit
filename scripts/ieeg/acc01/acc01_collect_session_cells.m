function sessions = acc01_collect_session_cells(sessionIds, loader)
%ACC01_COLLECT_SESSION_CELLS Preserve heterogeneous loader outputs by cell.
sessionIds = string(sessionIds(:));
sessions = cell(numel(sessionIds), 1);
for i = 1:numel(sessionIds)
    sessions{i} = loader(sessionIds(i));
end
end
