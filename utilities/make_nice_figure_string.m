function figure_string=make_nice_figure_string(string)
%make a string which is much nicer for labelling figures by:
%1. replacing '_' with ' ' which stops tex interpreting as superscript
%
%INPUT
%string - ugly string
%OUTPUT
%figure_string - nice string
%
%Author 
%Paddy Slator (p.slator@ucl.ac.uk)


%check if this is a cell of strings
if iscell(string)
    %if it is a cell, call the function for each string in the cell
    figure_string=cell(size(string));
    for i=1:length(string)
       figure_string{i}=make_nice_figure_string(string{i});
    end
else
    %remove x from beginning of string (matlab.lang.makeValidName adds x if the string begins with a number)
    if strcmp(string(1),'x')
        string(1)=[];
    end
    
    string_length=length(string);
    figure_string=string;
        
    for i=1:string_length
        if strcmp(figure_string(i),'_')
            if i>1
                if ~strcmp(figure_string(i-1),' ')
                    %replace underscore with space
                    figure_string(i)=' ';
                    %string_length=string_length+1;
                    %replace underscore with space
                    %figure_string(i:string_length)=[' ' figure_string(i+1:end)];
                end
            end
        end
    end
end

end