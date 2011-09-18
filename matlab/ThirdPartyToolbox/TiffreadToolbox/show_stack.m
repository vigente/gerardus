function hImage = show_stack( stack )

% function show_stack( mov )
%
% Display movie frames with a slider and autoscaling
%
% F. Nedelec, Jan 2009


hImage = show_image(image(stack, 1));
hAxes  = gca;

hFig = gcf;
fPos = get(hFig, 'Position');
set(hFig,  'Position', fPos + [0 0 0 24], 'Name', 'Images');
set(hFig,  'WindowButtonMotionFcn', {@callback_mouse});
set(hAxes, 'Units', 'pixels', 'Position', [1 25 fPos(3:4)]);

imx  = size(stack, ndims(stack));
indx = 1;

hSlider = uicontrol(hFig, 'Style', 'slider',...
    'Position',[ 20 0 fPos(3)-50 22 ],...
    'Min', 1, 'Max', imx, 'Value', 1,...
    'SliderStep', [ 1/imx, 2/imx ],...
    'Callback',{@callback_slider});

drawnow;


    function callback_slider(hObject, eventData)
        indx = round( get( hObject, 'Value' ) );
        show_image(image(stack,indx), hImage);
        set(hFig, 'Name', sprintf('Image %i', indx));
        drawnow;
    end

    function callback_mouse(hObject, eventData)
        i = round( get( hSlider, 'Value' ) );
        if i ~= indx
            indx = i;
            show_image(image(stack, indx), hImage);
            set(hFig, 'Name', sprintf('Image %i', indx));
            drawnow;
        end
    end

    function pix = image(stack, indx)
        pix = [];
        if ndims(stack) <= 2
            %compatibility with tiffread:
            if isfield(stack, 'data')  &&  ~iscell(stack(indx).data)
                pix = stack(indx).data;
            elseif isfield(stack, 'data') 
                pix = stack(indx).data{1};
            end
        elseif ndims(stack) == 3
            pix = stack(:,:,indx);
        end
    end

end

