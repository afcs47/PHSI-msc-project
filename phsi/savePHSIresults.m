function savePHSIresults(outputFolder, mainDatasetName, fusionChoice, polMethod, hsiMethod)
    %Save currently opened images
    
    if nargin < 5
      hsiMethod = '';
    end
    
    % Build filename based on chosen methods
    if ~isempty(hsiMethod)
        methodTag = sprintf('%s_%s_%s', fusionChoice, polMethod, hsiMethod);
    else
        methodTag = sprintf('%s_%s', fusionChoice, polMethod);
    end
    methodTag = strrep(methodTag, ' ', '_'); % remove spaces
    
    timestamp = datestr(now,'yyyy-mm-dd_HHMMSS');
    saveTag = sprintf('%s_%s_%s', mainDatasetName, methodTag, timestamp);
    
    outputFolder = fullfile(outputFolder, saveTag);
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    % Save all open figures
    figHandles = findall(0, 'Type', 'figure');
%     for i = 1:numel(figHandles) %figures saved as numbers
%         % png
%         figNamePng = sprintf('%s_Fig%d.png', saveTag, i);
%         saveas(figHandles(i), fullfile(outputFolder, figNamePng));
%     
%         % fig
%         figNameFig = sprintf('%s_Fig%d.fig', saveTag, i);
%         savefig(figHandles(i), fullfile(outputFolder, figNameFig));
%     end

    for i = 1:numel(figHandles)
        % Try to get figure title (use Name if no title is available)
        ax = get(figHandles(i), 'CurrentAxes');
        if ~isempty(ax)
            titleObj = get(ax, 'Title');
            figTitle = get(titleObj, 'String');
            if iscell(figTitle) % if multiline cell
                figTitle = strjoin(figTitle, '_');
            end
        else
            figTitle = get(figHandles(i), 'Name');
        end

        % Fallback if title is empty
        if isempty(figTitle)
            figTitle = sprintf('Figure_%d', i);
        end

        % Clean filename (remove invalid chars)
        figTitle = regexprep(figTitle, '[^\w\d-_]', '_');

        % png
        figNamePng = sprintf('%s.png', figTitle);
        saveas(figHandles(i), fullfile(outputFolder, figNamePng));

        % fig
        figNameFig = sprintf('%s.fig', figTitle);
        savefig(figHandles(i), fullfile(outputFolder, figNameFig));
    end
    
    fprintf('Saved %d figures to %s\n', numel(figHandles), outputFolder);
   
end