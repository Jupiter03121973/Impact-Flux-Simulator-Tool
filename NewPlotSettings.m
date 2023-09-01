function [xax,yax,XScale,YScale] = NewPlotSettings(opt_list)  
    bool = false;
    hfig=uifigure('Name','Plot Creation','Position',[100 100 232 240]);
    gl = uigridlayout(hfig,[7 2]);
    gl.RowHeight = {22,22,22,22,22,22,22,'1x'};
    gl.ColumnWidth = {100,100};

    %set defaults
    opt=opt_list{1};

    ScaleList = {'linear' 'log'};
    ScaleOpt = ScaleList{1};

    %create GUI           
    dropdown1=uidropdown(gl,'Items',opt_list,'Value',opt);
    dropdown2=uidropdown(gl,'Items',opt_list,'Value',opt);
    ddScale1=uidropdown(gl,'Items',ScaleList,'Value',ScaleOpt);
    ddScale2=uidropdown(gl,'Items',ScaleList,'Value',ScaleOpt);
    addbutton = uibutton(gl,'Text', 'Add','ButtonPushedFcn', @addplotbutton);
    cancelbutton=uibutton(gl,'Text', 'Cancel','ButtonPushedFcn', @cancelfun);
    label1=uilabel(gl,'Text','Variabel for x-Axis');
    label2=uilabel(gl,'Text','Variabel for y-Axis');
    ScaleLabel1=uilabel(gl,'Text','Scale for x-Axis');
    ScaleLabel2=uilabel(gl,'Text','Scale for y-Axis');

    dropdown1.Layout.Row = 1;
    dropdown1.Layout.Column = 2;
    dropdown2.Layout.Row = 2;
    dropdown2.Layout.Column = 2;
    ddScale1.Layout.Row = 4;
    ddScale1.Layout.Column = 2;
    ddScale2.Layout.Row = 5;
    ddScale2.Layout.Column = 2;
    addbutton.Layout.Row = 7;
    addbutton.Layout.Column = 1;
    cancelbutton.Layout.Row = 7;
    cancelbutton.Layout.Column = 2;
    label1.Layout.Row = 1;
    label1.Layout.Column = 1;
    label2.Layout.Row = 2;
    label2.Layout.Column = 1;
    ScaleLabel1.Layout.Row = 4;
    ScaleLabel1.Layout.Column = 1;
    ScaleLabel2.Layout.Row = 5;
    ScaleLabel2.Layout.Column = 1;

    uiwait(hfig);
    if bool == 1
        xax = dropdown1.Value;
        yax = dropdown2.Value;
        XScale = ddScale1.Value;
        YScale = ddScale2.Value;
    else
        xax = 0;
        yax = 0;
        XScale = 0;
        YScale = 0;
    end
    if isvalid(hfig)
        close(hfig)
    end
    
    function cancelfun(src,event)
        bool = 0;
        uiresume(hfig)
    end

    function addplotbutton(src,event)
        bool = 1;
        uiresume(hfig)
    end
end