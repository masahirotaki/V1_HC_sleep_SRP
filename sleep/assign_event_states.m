function labels = assign_event_states(times, states)
% ASSIGN_EVENT_STATES  Label event times by the behavioural state they fall in.
% Same precedence as Sleep_analysis_Ellie (later states overwrite earlier ones).
% Events outside every state interval are labelled "unclassified".

labels = repmat("unclassified", size(times));
for name = ["SWS", "REM", "quietWake", "movement"]
    iv = states.(name);
    for k = 1:size(iv, 1)
        labels(times >= iv(k, 1) & times <= iv(k, 2)) = name;
    end
end
end
