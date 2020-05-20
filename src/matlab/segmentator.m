function segmentator()

rootDir = fullfile('D:\', 'Home', 'Dimitris', ...
    'OneDrive - University College London', 'Data', ...
    'David', 'Midbrain Sequencing Data', 'fov');

my_dir = dir(rootDir);
my_dir = my_dir(~ismember({my_dir.name},{'.','..'}));
% my_dir = {my_dir.folder};
for i = 1: size(my_dir, 1)
    folderName = my_dir(i).folder;
    subfolderName = my_dir(i).name;
    img_name = [subfolderName, '.tif'];
    filePath = fullfile(folderName, subfolderName, 'img', img_name);
    
    fprintf('%s: Segmenting %s \n', datestr(now), img_name)
    try
        Dapi = imread(filePath);
    catch me
        fprintf('%s \n', me.message)
        continue
    end
        
    if max(Dapi(:)) == 0
        fprintf('%s: Image %s is totally black. Skipping segmentation\n', datestr(now), img_name)
        CellMap = Dapi;
        Boundaries = Dapi;
    else
        [CellMap, Boundaries] = my_segment_dapi(Dapi);
    end
    
    [~, myFilename, ~] = fileparts(img_name);
    label_image_path = save_csv(CellMap, folderName, subfolderName, myFilename, 'label_image');
    fprintf('%s: Label Image saved at %s \n', datestr(now), label_image_path)
    
    boundaries_path = save_csv(Boundaries, folderName, subfolderName, myFilename, 'boundaries');
    fprintf('%s: Boundaries saved at %s \n', datestr(now), boundaries_path)
    
    close all

end
end


function out = save_csv(arr, folderName, subfolderName, myFilename, id)
    % id is a string. 'label_image' or 'boundaries'
    save_path =  fullfile(folderName, subfolderName, id, [id, '_', myFilename, '.csv']);
    if ~exist(fileparts(save_path), 'dir')
       mkdir(fileparts(save_path))
    end
    writetable(array2table(arr), save_path, 'Delimiter', ',', 'WriteVariableNames', false)
    
    out = save_path;
end


