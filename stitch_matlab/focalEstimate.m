function [f1, f1Ok, f2, f2Ok] = focalEstimate(H)

% x1 = Hx2

f1Ok = true;
d1 = H(3, 1) * H(3, 2);
v1 = -(H(1, 1) * H(1, 2) + H(2, 1) * H(2, 2)) / d1;
d2 = (H(3, 2) - H(3, 1)) * (H(3,2) + H(3, 1));
v2 = (H(1, 1)^2 + H(2, 1)^2 - H(1, 2)^2 - H(2,2)^2) / d2;
if (v1 <= 0 && v2 <= 0)
    f1Ok = false;
    f1 = nan;
else
    if (v1 > 0 && v2 > 0)
        if (abs(d1) > abs(d2))
            f1 = sqrt(v1);
        else
            f1 = sqrt(v2);
        end
    else
        if (v1 > 0)
            f1 = sqrt(v1);
        else
            f1 = sqrt(v2);
        end
    end
end

f2Ok = true;
d1 = H(1, 1) * H(2, 1) + H(1 ,2) * H(2, 2);
v1 = -H(1, 3) * H(2, 3) / d1;
d2 = H(1, 1)^2 + H(1, 2)^2 - H(2, 1)^2 - H(2, 2)^2;
v2 = (H(2, 3)^2 - H(1, 3)^2) / d2;
if (v1 <= 0 && v2 <= 0)
    f2Ok = false;
    f2 = nan;
else
    if (v1 > 0 && v2 > 0)
        if (abs(d1) > abs(d2))
            f2 = sqrt(v1);
        else
            f2 = sqrt(v2);
        end
    else
        if (v1 > 0)
            f2 = sqrt(v1);
        else
            f2 = sqrt(v2);
        end
    end
end