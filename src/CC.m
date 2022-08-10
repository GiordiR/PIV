function [x, y, dx, dy] = CC (image1, image2, w_size, step, subpx)

    % Author: Riccardo Giordani
    %
    % Correlation analysis between two PIV images (Convolution)
    %
    %
    % arguments (input):
    %   image1 - first time PIV image
    %
    %   image2 - Second time PIV image
    %
    %   w_size - Interrogation window size
    %
    %   step - Distance between interrogation windows center points
    %
    %   subpx - Subpixel analysis mode (1=activated, 0=deactivated)
    %
    % arguments (output):
    %   x - x coordinate 
    %
    %   y - y coordinate 
    %
    %   dx - x displacement
    %
    %   dy - y displacement
    %

    

    image1 = double(image1);
    image2 = double(image2);


    % Windows center point limits
    ymin = 1+(ceil(w_size/2)); % NB! +1 per evitare errori di indici
    xmin = 1+(ceil(w_size/2));
    ymax = step*(floor(size(image1,1)/step))-(w_size-1)+(ceil(w_size/2)); % NB! -1 per evitare errori di indici 
    xmax = step*(floor(size(image1,2)/step))-(w_size-1)+(ceil(w_size/2));

    % Number of windows in the domain
    ycount = floor((ymax-ymin)/step+1);
    xcount = floor((xmax-xmin)/step+1);

    % Matrix expansion to allow search at borders
    % (adding w_size/2 at the matrix)
    image1 = padarray(image1,[ceil(w_size/2) ceil(w_size/2)], 0);
    image2 = padarray(image2,[ceil(w_size/2) ceil(w_size/2)], 0);



    % Matrix initialization
    x = zeros(ycount,xcount);
    y = zeros(ycount,xcount);
    dx = zeros(ycount,xcount);
    dy = zeros(ycount,xcount);

    % Matrix index counters initialization
    index_x = 0;
    index_x_real = 0;
    index_y = 0;


    for j = ymin:step:ymax % vertical loop
        index_y = index_y+1;
        for i = xmin:step:xmax % horizontal loop
            index_x = index_x+1;
            if index_x_real < xcount
                index_x_real = index_x_real+1;
            else
                index_x_real = 1;
            end

            startpoint=[i j];

            % Image position
            A = image1(j:j+w_size-1, i:i+w_size-1); % Interrogation windows overlap (w_size-step)
            B = image2(ceil(j-w_size/2):ceil(j+1.5*w_size-1), ceil(i-w_size/2):ceil(i+1.5*w_size-1));


            % Mean value subtraction to obtain an even number of positive
            % and negative values
            if mode(A,'all')==mean(A)  
                xpeak = w_size/2;
                ypeak = w_size/2;
            else
                A = A-mean(mean(A));
                B = B-mean(mean(B));
            

                % Discrete convolution (= correlation with conj)
                Conv = conv2(B,rot90(conj(A),2),'valid'); %      'valid' - returns only those parts of the convolution that are computed without the zero-padded edges              

                % Correlation normalization (peaks at 255)
                Conv = ((Conv-min(min(Conv)))/(max(max(Conv))-min(min(Conv))))*255;

                % Correlation peaks search
                [ypeak,xpeak] = find(Conv==255);

                % Select the first peak (in case of multiple peaks)
                if size(xpeak,1) > 1 
                  xpeak = xpeak(1:1); 
                end
                if size(ypeak,1) > 1
                    ypeak = ypeak(1:1);
                end
                
                % Subpixel level 
                if subpx==1
                    % 3x3 matrix extraction from correlation
                    if ypeak<size(Conv,1) && xpeak<size(Conv,2) && ypeak>1 && xpeak>1
                        Conv_sub = Conv(ypeak-1:ypeak+1,xpeak-1:xpeak+1);
                        if Conv_sub(2,1)==Conv_sub(2,3) && Conv_sub(1,2)==Conv_sub(3,2)
                            % centered peak 
                        else
                            if (rem(w_size,2) == 0) 
                                offset=1;
                            else
                                offset=0.5;
                            end
                            try
                                [displ] = subpx3gauss (Conv,w_size,xpeak,ypeak,offset);
                            catch
                                dispX = xpeak - (w_size/2);
                                dispY = ypeak - (w_size/2);
                                displ = [dispX, dispY]; 
                            end
                        end
                    end
                end
             end
            dispX = xpeak - (w_size/2);
            dispY = ypeak - (w_size/2);
            displ = [dispX, dispY];
            
            x(index_y,index_x_real) = startpoint(1) + w_size/2;
            y(index_y,:) = startpoint(1,2) + w_size/2;
            dx(index_y,index_x_real) = displ(1);
            dy(index_y,index_x_real) = displ(2);
        end
    end

    x = x-ceil(w_size/2);
    y = y-ceil(w_size/2);


    dx(dx>w_size/1.5) = NaN;
    dy(dx>w_size/1.5) = NaN;
    dy(dy>w_size/1.5) = NaN;
    dx(dy>w_size/1.5) = NaN;

end


function [displ] = subpx3gauss (Conv,w_size,x,y,offset)

    yphi0 = log(Conv(y,x));
    yphim1 = log(Conv(y-1,x));
    yphip1 = log(Conv(y+1,x));
    ypeak = y + (yphim1-yphip1)/(2*yphim1-4*yphi0+2*yphip1);

    xphi0 = log(Conv(y,x));
    xphim1 = log(Conv(y,x-1));
    xphip1 = log(Conv(y,x+1));
    xpeak = x + (xphim1-xphip1)/(2*xphim1-4*xphi0+2*xphip1);

    SubpxX = xpeak-(w_size/2)-offset;
    SubpzY = ypeak-(w_size/2)-offset;
    displ = [SubpxX, SubpzY];
    
end