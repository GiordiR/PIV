function image_out = rm_airfoil(image_in)

    % Author: Riccardo Giordani
    %
    % Raw method to eliminate airfoil from image
    % 
    %
    % arguments (input):
    %   image_in - Input image
    %  
    % arguments (output):
    %   image_out - Pre-processed image
    %


    % Windows dimention to seek consecutive white pixels
    D = 5;

    for j=1:size(image_in,2)
        for i=1:size(image_in,1)-D
            
            A = image_in(i:i+D,j);
            if mode(A)==255
                image_in(i:size(image_in,1),j) = NaN;
            end
        end
    end
    
    image_out = image_in;        

end