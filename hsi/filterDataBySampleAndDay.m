function [filteredData, selectedType, selectedDay] = filterDataBySampleAndDay(dataTable)
    % Select sample type (and day, when relevant)
    sampleTypes = unique(dataTable.SampleName);
    [selectedTypeIndex, tf] = listdlg('PromptString','Select sample type:', 'SelectionMode','single', 'ListString', sampleTypes);
    if ~tf, return; end
    selectedType = sampleTypes{selectedTypeIndex};

    dayList = unique(dataTable.Day(strcmp(dataTable.SampleName, selectedType)));
    if length(dayList)>1
        [selectedDayIndex, tf] = listdlg('PromptString','Select sample acquisition day:', 'SelectionMode','single', 'ListString', dayList);
        if ~tf, return; end
        selectedDay = dayList{selectedDayIndex};
    else
        selectedDay = dayList{1};
    end

    % Filter selected dataset
    filteredData = dataTable(strcmp(dataTable.SampleName, selectedType) & strcmp(dataTable.Day, selectedDay), :);
end