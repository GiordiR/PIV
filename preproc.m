function image_out = preproc (image_in,order,clahe,clahe_window,hpf,hpf_size,minmax,minmax_size,int_capt,capt_scaling)
    
    k = 1;
    while k<=size(order,2)
        l=order(k);
        switch l
            case 1
                % CLAHE (Contrast Limited Adaptive Histogram Equalization)
                % Clahe opera una equalizzazione dell'istogramma in piccole regioni
                % dell'immagine (tiles). In casi di PIV con un esposizione non uniforme, il
                % contrasto di ogni porzione di immagine viene ottimizzato
                % indipendentemente. I tiles confinanti vengono sucessivamente ricollegati
                % mediante una intepolazione bilineare
                if clahe == 1
                    ntiles1 = round(size(image_in,1)/clahe_window);
                    ntiles2 = round(size(image_in,2)/clahe_window);
                    image_out = adapthisteq(image_in, 'NumTiles',[ntiles1 ntiles2], 'ClipLimit', 0.01, 'NBins', 256, 'Range', 'full', 'Distribution', 'uniform');
                else
                end
                
            case 2
                % HIGH PASS FILTER
                % Le informazioni a bassa frequenza del background possono essere rimosse
                % applicando un filtro passa alto che conserva le frequenze alte relative
                % alle particelle. Questo filtro è calcolato applicando un filtro passa
                % basso all'immagine (blurring) e sottraendolo all'immagine originale.
                if hpf == 1
                    h = fspecial('gaussian',hpf_size,hpf_size);
                    image = double(image_in-(imfilter(image_in,h,'replicate')));
                    image_out = image/max(max(image))*255;
                else
                end
                
            case 3
                % MIN/MAX FILTER
                if minmax==1
            %         image_in = im2double(image_in);
                    image_in = double(image_in);
                    [y,x] = size(image_in);

                    % Inizializzazione degli envelope
                    lower = zeros(y,x);
                    upper = zeros(y,x);

                    neighb = floor(minmax_size/2); 
                    image_modified = padarray(image_in,[neighb neighb],NaN);
                    for i=1+neighb:y-neighb
                        for j=1+neighb:x-neighb
                            if rem(minmax_size,2)==0
                                disp('Errore: minmax_size deve essere un numero dispari');
                                break
                            else
                                tile = image_modified(i-neighb:i+neighb, j-neighb:j+neighb);

                                Min = min(tile,[],'all');
                                Max = max(tile,[],'all');

                                lower(i-neighb,j-neighb) = Min;
                                upper(i-neighb,j-neighb) = Max;
                            end
                        end         
                    end
                    % Moving average filter
                    kernel = ones(minmax_size,1) / minmax_size;
                    lower = filter(kernel, 1, lower);
                    upper = filter(kernel, 1, upper);

                    image_out = 255*(image_in - lower)./(upper - lower);
                else
                end
                
            case 4
                % INTENSITY CAPPING
                % Imposta i valori di intensità superiori ad un certo limite al limite
                % stesso. Il limite viene calcolato come
                %           I_lim = median(I) + n*std(I)
                % dove 0.5 < n < 2
                if int_capt==1

                    image_in = double(image_in);

                    I_lim = median(image_in,'all') + capt_scaling*std(image_in,0,'all','omitnan');

                    for i=1:size(image_in,1)
                        for j=1:size(image_in,2)
                            if image_in(i,j)>I_lim
                                image_in(i,j) = I_lim;
                            end
                        end
                    end
                    image_out = image_in;       
                else
                end
        end
    k = k+1;    
    end        
    

    % Filtri disattivati    
    if clahe==0 && hpf==0 && minmax==0 && int_capt==0
        image_out = image_in;
    else
    end

    image_out=uint8(image_out);
end