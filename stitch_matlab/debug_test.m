Dw = zeros(imgNum, 3);
Dthe = zeros(imgNum,  1);

for  i = 1:imgNum
    DR = hat_cams_waved(i).R;
    norm_w = acos((trace(DR) - 1)/2);
    Dw(i, :) = (1/(2 * sin(norm_w))) .* [DR(3, 2) - DR(2, 3), DR(1, 3) - DR(3,1), DR(2,1)-DR(1,2)];
    Dthe(i) = norm_w;
end

figure;
hold all;

for i = 1:imgNum
    plot3([0, Dw(i,1)], [0, Dw(i, 2)], [0, Dw(i,3)]);
end
