function [x, y, dx, dy] = CC (image1,image2,w_size, step, subpx)


    image1 = double(image1);
    image2 = double(image2);


    % Posizione limite dei centri delle finestre
    ymin = 1+(ceil(w_size/2)); % NB! +1 per evitare errori di indici
    xmin = 1+(ceil(w_size/2));
    ymax = step*(floor(size(image1,1)/step))-(w_size-1)+(ceil(w_size/2)); % NB! -1 per evitare errori di indici 
    xmax = step*(floor(size(image1,2)/step))-(w_size-1)+(ceil(w_size/2));

    % Numero di finestre nel dominio
    ycount = floor((ymax-ymin)/step+1);
    xcount = floor((xmax-xmin)/step+1);

    % Espansione delle matrici delle immagini per permettere la ricerca ai
    % bordi (aggiunta di un contorno di dimensione w_size/2 alla matrice)
    image1 = padarray(image1,[ceil(w_size/2) ceil(w_size/2)], 0);
    image2 = padarray(image2,[ceil(w_size/2) ceil(w_size/2)], 0);



    % Inizializzazione matrici
    x = zeros(ycount,xcount);
    y = zeros(ycount,xcount);
    dx = zeros(ycount,xcount);
    dy = zeros(ycount,xcount);

    % Inizializzazione conteggi per posizione nella matrice
    index_x = 0;
    index_x_real = 0;
    index_y = 0;


    for j = ymin:step:ymax % loop verticale
        index_y = index_y+1;
        for i = xmin:step:xmax % loop orizzontale
            index_x = index_x+1;
            if index_x_real < xcount
                index_x_real = index_x_real+1;
            else
                index_x_real = 1;
            end

            startpoint=[i j];

            % Porzioni delle immagini
            A = image1(j:j+w_size-1, i:i+w_size-1); % Finestre di interrogazione si sovrappongono per (w_size-step)
            B = image2(ceil(j-w_size/2):ceil(j+1.5*w_size-1), ceil(i-w_size/2):ceil(i+1.5*w_size-1));


            % Sottrazione del valore medio in modo da avere a grandi linee un
            % numero uguale di valori negativi e positivi
            if mode(A,'all')==mean(A)  
                xpeak = w_size/2;
                ypeak = w_size/2;
            else
                A = A-mean(mean(A));
                B = B-mean(mean(B));
            

                % Convoluzione discreta (= correlazione usando conj)
                Conv = conv2(B,rot90(conj(A),2),'valid'); %      'valid' - returns only those parts of the convolution that are computed without the zero-padded edges              

                % Normalizzazione della correlazione (picchi a 255)
                Conv = ((Conv-min(min(Conv)))/(max(max(Conv))-min(min(Conv))))*255;

                % Ricerca picchi della correlazione
                [ypeak,xpeak] = find(Conv==255);

                % Seleziona solo il primo picco (in caso di picchi multipli)
                if size(xpeak,1) > 1 
                  xpeak = xpeak(1:1); 
                end
                if size(ypeak,1) > 1
                    ypeak = ypeak(1:1);
                end
                
                % Verifica a livello subpixel
                if subpx==1
                    % Estrazione matrice 3x3 nella correlazione
                    if ypeak<size(Conv,1) && xpeak<size(Conv,2) && ypeak>1 && xpeak>1
                        Conv_sub = Conv(ypeak-1:ypeak+1,xpeak-1:xpeak+1);
                        if Conv_sub(2,1)==Conv_sub(2,3) && Conv_sub(1,2)==Conv_sub(3,2)
                            % picco centrato  
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