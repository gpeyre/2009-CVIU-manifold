% batch mode


name_list = {'irreg_stones', 'violettes', 'nuts', 'bonbons', ...
    'monet', 'mousse', ...
};



name_list = {
    'cherries',...
    'pieuvre',...
    'thai_art',...
    'town',...
};

name_list = {
    'rafia', ...
    'peoples',...
    'reptil_skin', ...
    'mures',...
    'sea_grid', ...
    'crochetcropti',...
    'empilements', ...
    'radishes',...
    'warped_grid', ...
    'olives',...
    'caustiques',...
    'ondulations',...
    'pasta',...
    'frenchfries',...
    'mures',...
    'framboise', ...
    'brick',...
    'tomatoes',...
    'stones',...
    'dunes', ...
};

name_list = {
      'apples', ...
      'caustiques', ...
      'cherries', ...
      'chocolate', ...
      'nuts', ...
      'thai_art', ...
      'town', ...
      'yellow_peppers', ...
    'pointilliste', ...
    'text', ...
    'corral', ...
    'zigzag', ...  
};

name_list = {
      'cherries', ...
    'pointilliste', ...
    'radishes',...
    'mures',...
};


for i=1:length(name_list)
    clear name;
    name = name_list{i};
    disp(['--- Synthesizing ' name ' ---']);
    test_synthesis;
end