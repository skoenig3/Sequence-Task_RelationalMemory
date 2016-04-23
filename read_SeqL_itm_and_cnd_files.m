function [itmlist,sequence_items,sequence_locations] = read_SeqL_itm_and_cnd_files(itmfile,monkey)
% writteen by Seth Konig August, 2014
% Function imports item file and and grabs condition file to determine which items
% are associated with which condition (itmlist) since conditions are randomly
% organized. Further function determines the number of sequences and which
% items (the largets item #) are assoicated with which sequence.

ITMFile = ['C:\Users\seth.koenig\Documents\MATLAB\Sequence Task_RelationalMemory\Item and Condition Files\' itmfile];
if strcmpi(monkey,'TT') && strcmpi(itmfile(1:end-4),'SeqL11')
    %run accintally on the wrong condition file
    CNDFile = ['C:\Users\seth.koenig\Documents\MATLAB\Sequence Task_RelationalMemory\Item and Condition Files\SeqL01.cnd'];
else
    CNDFile = ['C:\Users\seth.koenig\Documents\MATLAB\Sequence Task_RelationalMemory\Item and Condition Files\' itmfile(1:end-4) '.cnd'];
end


itmfil=[];
h =fopen(ITMFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(itmfil)
        if length(tline)>size(itmfil,2)
            tline=tline(1:size(itmfil,2));
        end
    end
    tline = [tline ones(1,(size(itmfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        itmfil=[itmfil; tline];
    else
        break
    end
end
fclose(h);

cndfil=[];
h=fopen(CNDFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(cndfil)
        if length(tline)>size(cndfil,2)
            tline=tline(1:size(cndfil,2));
        end
    end
    tline = [tline ones(1,(size(cndfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        cndfil=[cndfil; tline];
    else
        break
    end
end
fclose(h);

first_img_item = [];
for i = 9:size(itmfil,1); %first 9 are hearder, background, and color change
    str = textscan(itmfil(i,:),'%s');
    if ~isempty(strfind(str{1}{end},'.bmp')) && isempty(first_img_item)
        first_img_item  = str2num((str{1}{1}));
        break
    end
end

items_per_seq = 8; %4 for actual shapes and 4 for displaying fixation window
total_seq_items = first_img_item-1-3;%3 items for color change
sequence_items = [3+items_per_seq*(1:total_seq_items/items_per_seq)];


itmlist = zeros(size(cndfil,1)-1,1);
for i = 2:size(cndfil,1);
    str = textscan(cndfil(i,:),'%d');
    itmlist(i-1) = str{1}(end);
end

sequence_locations = cell(1,length(sequence_items));
for seq = 1:length(sequence_items)
    items = 9+(seq-1)*items_per_seq+(1:4);
    for i = 1:length(items)
        str = textscan(itmfil(items(i),:),'%d');
        sequence_locations{seq}(:,i) =  24*double(str{1}(4:5))+[400;300]; %turn into pixel coordinates
    end
end

