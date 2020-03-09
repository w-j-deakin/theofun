function [angD,angR] = findAngle2D(vec1,vec2)

cosTheta = dot(vec1,vec2)/(norm(vec1)*norm(vec2));

angD = acosd(cosTheta);
angR = acos(cosTheta);


end

