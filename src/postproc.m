function [u_valid,v_valid,V_valid] = postproc(u,v,V,method,interp,parameter)

    % Author: Riccardo Giordani
    %
    % Post-processing correlation analysis
    %
    %
    % arguments (input):
    %   u - x velocity field
    %
    %   v - y velocity field
    %
    %   method - Velocity post-processing method 
    %           (0=manual limits, 1=Vector difference test, 2=std dev test, 3=Normalized median test)
    %
    %   interp - Missing velocity interpolation
    %           (0=deactivated, 1=activated)
    %
    %   parameter - Parameter struct for test method
    %
    % arguments (output):
    %   u_valid - x coordinate velocities after post-processing
    %
    %   v_valid - y coordinate velocities after post-processing
    %
    %   V_valid - Velocity field after post-processing
    %

    % Outliers removal
    u_old = u;
    v_old = v;

    neighb = floor(parameter.tile_size/2);
    
    if method.on_uv == 1
        switch method.choice
            case 0
                % Manual velocity limits 
                umin = parameter.u_lim(1); 
                umax = parameter.u_lim(2); 
                vmin = parameter.v_lim(1); 
                vmax = parameter.v_lim(2); 

                % Check compatibility with limits
                u(u<umin) = NaN;
                u(u>umax) = NaN;
                v(v<vmin) = NaN;
                v(v>vmax) = NaN;

            case 1
                % Vector Difference Test
                if parameter.u_diff==0 && parameter.v_diff==0
                else
                    for i=neighb+1:size(u,1)-neighb
                        for j=neighb+1:size(u,2)-neighb

                            % Cell selection
                            tile_u = u(i-neighb:i+neighb,j-neighb:j+neighb);
                            tile_v = v(i-neighb:i+neighb,j-neighb:j+neighb);

                            count = 0;
                            for m=1:parameter.tile_size
                                for n=1:parameter.tile_size
                                    % Test
                                    u_diff = abs(tile_u(m,n)-u(i,j));
                                    v_diff = abs(tile_v(m,n)-v(i,j));

                                    if u_diff>parameter.u_diff && v_diff>parameter.v_diff               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        count = count + 1;
                                    else
                                    end
                                end
                            end
                            if count>=(parameter.tile_size*parameter.tile_size-1)/2
                                u(i,j) = NaN;
                                v(i,j) = NaN;
                            end
                        end
                    end
                end
            case 2
                % Standard deviation filter
                     for i=neighb+1:size(u,1)-neighb
                        for j=neighb+1:size(u,2)-neighb

                            % Cell selection
                            tile_u = u(i-neighb:i+neighb,j-neighb:j+neighb);
                            tile_v = v(i-neighb:i+neighb,j-neighb:j+neighb);

                            mean_u = mean(tile_u,'all');
                            mean_v = mean(tile_v,'all');

                            std_u = std(tile_u,0,'all','omitnan');
                            std_v = std(tile_v,0,'all','omitnan');

                            u_lower = mean_u - parameter.tile_size*std_u;
                            u_upper = mean_u + parameter.tile_size*std_u;
                            v_lower = mean_v - parameter.tile_size*std_v;
                            v_upper = mean_v + parameter.tile_size*std_v;

                            if u(i,j)<u_lower && v(i,j)<v_lower
                                u(i,j) = NaN;
                                v(i,j) = NaN;
                            end
                            if u(i,j)>u_upper && v(i,j)>v_upper
                                u(i,j) = NaN;
                                v(i,j) = NaN;
                            end

                        end
                    end 
                 
            case 3
                % Normalized Median filter
                if parameter.builtin==0
                    for i=neighb+1:size(u,1)-neighb
                        for j=neighb+1:size(u,2)-neighb

                            % Cell selection
                            tile_u = u(i-neighb:i+neighb,j-neighb:j+neighb);
                            tile_v = v(i-neighb:i+neighb,j-neighb:j+neighb);

                            % Reference velocity for the cell
                            u_ref = median(tile_u,'all');
                            v_ref = median(tile_v,'all');

                            % Compute residuals
                            residual_u = tile_u - u_ref;
                            residual_v = tile_v - v_ref;

                            % Residual normalization
                            r_norm_u = abs(u(i,j) - u_ref)/(median(residual_u,'all')+parameter.epsilon);
                            r_norm_v = abs(v(i,j) - v_ref)/(median(residual_v,'all')+parameter.epsilon);

                            if r_norm_u>parameter.res_toll
                                u(i,j) = NaN;
                            end   
                            if r_norm_v>parameter.res_toll
                                v(i,j) = NaN;
                            end 
                        end
                    end
                elseif parameter.builtin==1
                    u = medfilt2(u);
                    v = medfilt2(v);
                end
               
        end        
        u_valid = u;
        v_valid = v;
        V_valid = sqrt(u_valid.^2+v_valid.^2);
        
        % Accurancy measure (Outliers percentage)
        n_u_outliers = nnz(isnan(u));
        n_v_outliers = nnz(isnan(v));
        
        n_u_tot = size(u,1)*size(u,2);
        n_v_tot = size(v,1)*size(v,2);
        
        error_u = n_u_outliers/n_u_tot*100;
        error_v = n_v_outliers/n_v_tot*100;
        
        fprintf('u outlier percentage: %f %% \n',error_u);
        fprintf('v outlier percentage: %f %% \n',error_v);   
        
        % Missing velocity interpolation
        if interp.uv==1
            u_valid = inpaint_nans(u_valid,4);
            v_valid = inpaint_nans(v_valid,4);
            V_valid = sqrt(u_valid.^2 + v_valid.^2);
        end
        
    elseif method.on_uv == 0       
        switch method.choice
            case 0
                % Manual velocity limits
                Vmin = parameter.V_lim(1);
                Vmax = parameter.V_lim(2);

                % Check limit compatibility
                V(V<Vmin) = NaN;
                V(V>Vmax) = NaN;

            case 1
                % Vector Difference Test
                if parameter.V_diff==0
                else
                    for i=neighb+1:size(V,1)-neighb
                        for j=neighb+1:size(V,2)-neighb

                            % Cell selection
                            tile_V = V(i-neighb:i+neighb,j-neighb:j+neighb);

                            count = 0;
                            for m=1:parameter.tile_size
                                for n=1:parameter.tile_size
                                    % Test
                                    V_diff = abs(tile_V(m,n)-V(i,j));
                                    if V_diff>parameter.V_diff
                                        count = count + 1;
                                    else
                                    end
                                end
                            end
                            if count>=(parameter.tile_size*parameter.tile_size-1)/2
                                V(i,j) = NaN;
                            end
                        end
                    end
                end
            case 2
                % Standard deviation filter
                     for i=neighb+1:size(V,1)-neighb
                        for j=neighb+1:size(V,2)-neighb

                            % Cell selection
                            tile_V = V(i-neighb:i+neighb,j-neighb:j+neighb);

                            mean_V = mean(tile_V,'all');

                            std_V = std(tile_V,0,'all','omitnan');

                            V_lower = mean_V - parameter.tile_size*std_V;
                            V_upper = mean_V + parameter.tile_size*std_V;

                            if V(i,j)<V_lower
                                V(i,j) = NaN;
                            end
                            if V(i,j)>V_upper
                                V(i,j) = NaN;
                            end                       

                        end
                    end
                  
            case 3
                % Normalized Median filter
                if parameter.builtin == 0
                   for i=neighb+1:size(V,1)-neighb
                        for j=neighb+1:size(V,2)-neighb

                            % Cell selection
                            tile_V = V(i-neighb:i+neighb,j-neighb:j+neighb);

                            % Reference velocity for the cell
                            V_ref = median(tile_V,'all');

                            % Compute residual
                            residual_V = tile_V - V_ref;

                            % Residual normalization
                            r_norm_V = abs(V(i,j) - V_ref)/(median(residual_V,'all')+parameter.epsilon);

                            if r_norm_V>parameter.res_toll
                                V(i,j) = NaN;
                            end   
                        end
                   end 
                elseif parameter.builtin==1
                    V = medfilt2(V);                  
                end
               
        end
        u_valid = u_old;
        v_valid = v_old;
        V_valid = V;
        
        % Accuracy measure (Outlier percentage)
        n_V_outliers = nnz(isnan(V));
        
        n_V_tot = size(V,1)*size(V,2);
        
        error_V = n_V_outliers/n_V_tot*100;
        
        fprintf('V outlier percentage: %f %% \n',error_V);
        
        % Missing velocity interpolation
        if interp.V==1
        V_valid = inpaint_nans(V_valid,4);
        end        
    else
        disp('Error! method.on_uv must be 0 or 1!');
    end
       

end

