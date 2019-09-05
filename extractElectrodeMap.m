function elecMap = extractElectrodeMap(cds)
    % Function to extract the electrode map from the CDS.
    %
    % Outputs a table with columns:
    %   monkey: monkey name (string)
    %   array: name of array (string, ex. 'LeftS1Area2')
    %   chan: channel number (int)
    %   electrode label (ex. 'elec78')
    %   rowNum: index of row, labeled 1-10 from bottom to top (from pad side with wire bundle on the right)
    %   colNum: index of col, labeled 1-10 from left to right (from pad side with wire bundle on the right)

    elecMap = struct2table(cds.units);
    elecMap = elecMap(:,{...
        'monkey',...
        'array',...
        'chan',...
        'label',...
        'rowNum',...
        'colNum'});
    elecMap = unique(elecMap);


    elecMap.Properties.VariableDescriptions = {...
        'monkey name',...
        'name of array',...
        'channel number',...
        'electrode label',...
        'index of row, labeled 1-10 from bottom to top (from pad side with wire bundle on the right)',...
        'index of col, labeled 1-10 from left to right (from pad side with wire bundle on the right)'};
