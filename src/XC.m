function [x, y, dx, dy] = XC(image1,image2,w_size, wc_step,subpx)

    % Author: Riccardo Giordani
    %
    % Correlation analysis between two PIV images (Cross-correlation)
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



    %     image1 = double(image1);
    %     image2 = double(image2);


    % Windows center point limits
    ymin = 1+(ceil(w_size/2)); % NB! +1 per evitare errori di indici
    xmin = 1+(ceil(w_size/2));
%     ymax = wc_step*(floor(size(image1,1)/wc_step))-(w_size-1)+(ceil(w_size/2)); % NB! -1 per evitare errori di indici 
%     xmax = wc_step*(floor(size(image1,2)/wc_step))-(w_size-1)+(ceil(w_size/2));
    ymax = size(image1,1) - w_size/2 -1;
    xmax = size(image1,2) - w_size/2 -1;

    % Number of windows in the domain
    ycount = floor((ymax-ymin)/wc_step+1);
    xcount = floor((xmax-xmin)/wc_step+1);

    % Matrix initialization
    x = zeros(ycount,xcount);
    y = zeros(ycount,xcount);
    dx = zeros(ycount,xcount);
    dy = zeros(ycount,xcount);

    % Matrix index counters initialization
    ix = 0;
    ixx = 0;
    iy = 0;


    for j = ymin:wc_step:ymax % vertical loop  
        iy = iy+1;
        for i = xmin:wc_step:xmax % horizontal loop     
            ix = ix+1;
            if ixx < xcount
                ixx = ixx+1;
            else
                ixx = 1;
            end

            startpoint=[i j];

            % Image position
            A = image1(j-w_size/2:j+w_size/2, i-w_size/2:i+w_size/2); 
            B = image2(j-w_size/2:j+w_size/2, i-w_size/2:i+w_size/2);
            
            xc_corr = (size(A,2) + size(B,2))/2;
            yc_corr = (size(A,1) + size(B,1))/2;

            % Check if it we are inside a glaring part of the image (where we cannot
            % perform correlations)
            if mode(A,'all')==mean(A)
                xpeak = (size(A,2) + size(B,2))/2;
                ypeak = (size(A,1) + size(B,1))/2;
            else
                % Mean value subtraction to obtain an even number of positive
                % and negative values
                A = A-mean(mean(A));
                B = B-mean(mean(B));

                % Cross-correlation
%                 corr = xcorr2(B,A);
                corr = normxcorr2(A,B);
                
                [ypeak, xpeak] = find(corr==max(corr(:)));

            
            
                % Select the first peak (in case of multiple peaks)
                if size(xpeak,1) > 1 
                  xpeak = xpeak(1); 
                end
                if size(ypeak,1) > 1
                    ypeak = ypeak(1);
                end


                % Subpixel level 
                % 3x3 matrix extraction from correlation
                if ypeak<size(corr,1) && xpeak<size(corr,2) && ypeak>1 && xpeak>1
                    corr_sub = corr(ypeak-1:ypeak+1,xpeak-1:xpeak+1);
                    if corr_sub(2,1)==corr_sub(2,3) && corr_sub(1,2)==corr_sub(3,2)
                        % centered peak 
                    else
                        try
                            if subpx==1
                                ypeak_old = ypeak;
                                xpeak_old = xpeak;
                                [ypeak,xpeak] = subpx3gauss (corr_sub,xpeak,ypeak);
                                if isreal(ypeak) && isreal(xpeak)
                                else
                                    ypeak = ypeak_old;
                                    xpeak = xpeak_old;
                                end
                            end
                        catch
                        end
                    end
                else
                end
            end
            dispX = xpeak - xc_corr;
            dispY = ypeak - yc_corr;
            displ = [dispX, dispY]; 

%             x(iy,ixx) = startpoint(1) + w_size/2;
%             y(iy,:) = startpoint(1,2) + w_size/2;
            x(iy,ixx) = startpoint(1) + wc_step;
            y(iy,:) = startpoint(1,2) + wc_step;
            dx(iy,ixx) = displ(1);
            dy(iy,ixx) = displ(2);
        end
    end


    x = x-ceil(w_size/2);
    y = y-ceil(w_size/2);

    % Check if displacements are bigger than window dimension
    % (interrogation window)
    dx(dx>w_size/1.5) = NaN;
    dy(dx>w_size/1.5) = NaN;
    dy(dy>w_size/1.5) = NaN;
    dx(dy>w_size/1.5) = NaN;

end


function [ypeakG,xpeakG] = subpx3gauss (corr_sub,xpeak,ypeak)

        xcenter = ceil(size(corr_sub,2)/2);
        ycenter = ceil(size(corr_sub,1)/2);
        
        yphi0 = log(corr_sub(ycenter,xcenter));
        yphim1 = log(corr_sub(ycenter-1,xcenter));
        yphip1 = log(corr_sub(ycenter+1,xcenter));
        ypeakG = ypeak + 0.5*(yphim1-yphip1)/(yphim1-2*yphi0+yphip1);
        
        xphi0 = log(corr_sub(ycenter,xcenter));
        xphim1 = log(corr_sub(ycenter,xcenter-1));
        xphip1 = log(corr_sub(ycenter,xcenter+1));
        xpeakG = xpeak + 0.5*(xphim1-xphip1)/(xphim1-2*xphi0+xphip1);
end