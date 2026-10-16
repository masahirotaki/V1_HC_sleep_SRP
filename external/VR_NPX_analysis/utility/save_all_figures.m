function save_all_figures(save_path,filename,varargin)
% Save figures in .fig and .png format
% Default values
p = inputParser;
addParameter(p,'ContentType','vector',@isstr)
addParameter(p,'SVG_option',[])

% assign parameters (either defaults or given)
parse(p,varargin{:});
ContentType = p.Results.ContentType;
SVG_option = p.Results.SVG_option;



figlist = findobj('type','figure');

for i = 1 : numel(figlist)
    if exist(filename,'var')
        name = filename;
    else 
        name = get(figlist(i),'Name');
        if contains(name,'/') 
            name(strfind(name,'/')) = ';';
        elseif contains(name,'\')
            name(strfind(name,'\')) = ';';
        elseif contains(name,':')
            name(strfind(name,':')) = '-';
        end
        if isempty(name)
           disp('Figure has no name!')
           name = ['Figure_' num2str(i)];
           if exist(name,'file')
                  name = ['Figure_' num2str(i) '_A'];
           end
        end
    end 
    if exist([save_path,'\figures'], 'dir')
        saveas(figlist(i),[save_path,'\figures\',name,'.fig']);
%         saveas(figlist(i),[save_path,'\figures\png_figs\',name,'.png']);
%         saveas(figlist(i),[save_path,'\figures\png_figs\',name,'.pdf']);
        exportgraphics(figlist(i),[save_path,'\figures\png_figs\',name,'.pdf'],'ContentType',ContentType)
        % exportgraphics(figlist(i),[save_path,'\figures\png_figs\',name,'.pdf'],'ContentType','image')
    else 
        disp('not saved in figures folder!')
        saveas(figlist(i),[save_path,'\',name,'.fig']);
%         saveas(figlist(i),[save_path,'\',name,'.png']);
%         saveas(figlist(i),[save_path,'\',name,'.pdf']);

         % exportgraphics(figlist(i),[save_path,'\',name,'.pdf'],'ContentType','image')
        if ~isempty(SVG_option)
            print(figlist(i), [save_path,'\',name,'.svg'], '-dsvg', '-painters');
            exportgraphics(figlist(i),[save_path,'\',name,'.pdf'],'ContentType','image')
        else
            exportgraphics(figlist(i),[save_path,'\',name,'.pdf'],'ContentType',ContentType)
        end
    end
    close  
end

end