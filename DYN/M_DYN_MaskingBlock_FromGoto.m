clear 'gotoTagStack'
blocks = get_param(gcs, 'Blocks');
colorPalette = colormap('colorcube');
close(gcf);

for ii = 1:length(blocks)
    if (contains(blocks{ii}, 'Goto'))
%         idNum = str2double(regexp(blocks{ii}, '\d*', 'match'));
% 		if (isempty(idNum))
% 			idNum = 0;
% 		end

        % Identify block name
        blockName = [gcs, '/', blocks{ii}];
        
        % Identify goto tag
        gotoTag = get_param(blockName, 'GotoTag');
        if exist('gotoTagStack', 'var')
            if any(matches(cellstr(gotoTagStack), gotoTag))
                idNum = find(matches(cellstr(gotoTagStack), gotoTag));
            else
                gotoTagStack = [gotoTagStack, gotoTag];
                % Identify ID number
                idNum = length(gotoTagStack);
            end
        else
            gotoTagStack = {gotoTag};
            
            % Identify ID number
            idNum = length(gotoTagStack);
            
        end

        % Set Background/Foreground Colors
        colorNum = mod(idNum * 18 + 1, 256);
        colorValue = colorPalette(colorNum ,:);
        colorValueInverse = [1, 1, 1] - colorValue;
		set_param(blockName, 'BackgroundColor', ['[', num2str(colorValue), ']']);
        set_param(blockName, 'ForegroundColor', ['[', num2str(colorValueInverse), ']']);
        set_param(blockName, 'FontSize', '11', 'FontWeight', 'Bold');
        
        % Adjust Size
        prevPos = get_param(blockName, 'Position');
        newPos = [prevPos(1), prevPos(2), prevPos(1)+100, prevPos(2)+20];
        set_param(blockName, 'Position', ['[', num2str(newPos), ']']);

    elseif (contains(blocks{ii}, 'From'))
% 		idNum = str2double(regexp(blocks{ii}, '\d*', 'match'));
% 		if (isempty(idNum))
% 			idNum = 0;
%         end
        
        % Identify block name
        blockName = [gcs, '/', blocks{ii}];
        
        % Identify goto tag
        gotoTag = get_param(blockName, 'GotoTag');
        if exist('gotoTagStack', 'var')
            if any(matches(cellstr(gotoTagStack), gotoTag))
                idNum = find(matches(cellstr(gotoTagStack), gotoTag));
            else
                gotoTagStack = [gotoTagStack, gotoTag];
                % Identify ID number
                idNum = length(gotoTagStack);
            end
        else
            gotoTagStack = {gotoTag};
            
            % Identify ID number
            idNum = length(gotoTagStack);
        end
        
        % Set Background/Foreground Colors
        colorNum = mod(idNum * 18 + 1, 256);
        colorValue = colorPalette(colorNum ,:);
        colorValueInverse = [1, 1, 1] - colorValue;
		set_param(blockName, 'BackgroundColor', ['[', num2str(colorValue), ']']);
        set_param(blockName, 'ForegroundColor', ['[', num2str(colorValueInverse), ']']);
        set_param(blockName, 'FontSize', '11', 'FontWeight', 'Bold');
        
        % Adjust Size
        prevPos = get_param(blockName, 'Position');
        newPos = [prevPos(1), prevPos(2), prevPos(1)+100, prevPos(2)+20];
        set_param(blockName, 'Position', ['[', num2str(newPos), ']']);
        
    end


end