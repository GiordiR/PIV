function [x, y, dx, dy] = XC(image1,image2,w_size, wc_step,subpx)


%     image1 = double(image1);
%     image2 = double(image2);


    % Posizione limite dei centri delle finestre
    ymin = 1+(ceil(w_size/2)); % NB! +1 per evitare errori di indici
    xmin = 1+(ceil(w_size/2));
%     ymax = wc_step*(floor(size(image1,1)/wc_step))-(w_size-1)+(ceil(w_size/2)); % NB! -1 per evitare errori di indici 
%     xmax = wc_step*(floor(size(image1,2)/wc_step))-(w_size-1)+(ceil(w_size/2));
    ymax = size(image1,1) - w_size/2 -1;
    xmax = size(image1,2) - w_size/2 -1;

    % Numero di finestre nel dominio
    ycount = floor((ymax-ymin)/wc_step+1);
    xcount = floor((xmax-xmin)/wc_step+1);

    % Inizializzazione matrici
    x = zeros(ycount,xcount);
    y = zeros(ycount,xcount);
    dx = zeros(ycount,xcount);
    dy = zeros(ycount,xcount);

    % Inizializzazione conteggi per posizione nella matrice finale
    ix = 0;
    ixx = 0;
    iy = 0;


    for j = ymin:wc_step:ymax % loop verticale   
        iy = iy+1;
        for i = xmin:wc_step:xmax % loop orizzontale     
            ix = ix+1;
            if ixx < xcount
                ixx = ixx+1;
            else
                ixx = 1;
            end

            startpoint=[i j];

            % Porzioni delle immagini
            A = image1(j-w_size/2:j+w_size/2, i-w_size/2:i+w_size/2); 
            B = image2(j-w_size/2:j+w_size/2, i-w_size/2:i+w_size/2);
            
            xc_corr = (size(A,2) + size(B,2))/2;
            yc_corr = (size(A,1) + size(B,1))/2;

            % Verifico di non essere in una zona con riflessi spuri in cui non
            % posso fare la correlazione
            if mode(A,'all')==mean(A)
                xpeak = (size(A,2) + size(B,2))/2;
                ypeak = (size(A,1) + size(B,1))/2;
            else
                % Sottrazione del valore medio in modo da avere a grandi linee un
                % numero uguale di valori negativi e positivi
                A = A-mean(mean(A));
                B = B-mean(mean(B));

                % Cross-correlazione
%                 corr = xcorr2(B,A);
                corr = normxcorr2(A,B);
                
                [ypeak, xpeak] = find(corr==max(corr(:)));

            
            
                % Seleziona solo il primo picco (in caso di picchi multipli)
                if size(xpeak,1) > 1 
                  xpeak = xpeak(1); 
                end
                if size(ypeak,1) > 1
                    ypeak = ypeak(1);
                end


                % Verifica a livello subpixel
                % Estrazione matrice 3x3 nella correlazione
                if ypeak<size(corr,1) && xpeak<size(corr,2) && ypeak>1 && xpeak>1
                    corr_sub = corr(ypeak-1:ypeak+1,xpeak-1:xpeak+1);
                    if corr_sub(2,1)==corr_sub(2,3) && corr_sub(1,2)==corr_sub(3,2)
                        % picco centrato 
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

    % Verifico che gli spostamenti non superino la dimensione della finestra
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