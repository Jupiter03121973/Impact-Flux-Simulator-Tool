function DelPlot = Delete_Plot(opt_list)  
    bool = false;
    hfig=uifigure('Name','Delete Plot','Position',[100 100 232 140]);
    gl = uigridlayout(hfig,[3 2]);
    gl.RowHeight = {22,22,22,'1x'};
    gl.ColumnWidth = {100,100};

    %set defaults
    opt=opt_list{1};

    %create GUI           
    dropdown1=uidropdown(gl,'Items',opt_list,'Value',opt);
    deletebutton = uibutton(gl,'Text', 'Delete','ButtonPushedFcn', @deleteplotbutton);
    cancelbutton=uibutton(gl,'Text', 'Cancel','ButtonPushedFcn', @cancelfun);
    label1=uilabel(gl,'Text','Plot to delete');

    dropdown1.Layout.Row = 1;
    dropdown1.Layout.Column = 2;
    deletebutton.Layout.Row = 4;
    deletebutton.Layout.Column = 1;
    cancelbutton.Layout.Row = 4;
    cancelbutton.Layout.Column = 2;
    label1.Layout.Row = 1;
    label1.Layout.Column = 1;
    
    uiwait(hfig);
    if bool == 1
        DelPlot = dropdown1.Value;
    else
        DelPlot = 0;
    end
    if isvalid(hfig)
        close(hfig)
    end
    
    function cancelfun(src,event)
        bool = 0;
        uiresume(hfig)
    end

    function deleteplotbutton(src,event)
        bool = 1;
        uiresume(hfig)
    end
end