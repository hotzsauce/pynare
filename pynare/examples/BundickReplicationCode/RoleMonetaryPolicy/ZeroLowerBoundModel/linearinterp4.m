function [zi] = linearinterp4(nx1,x1grid,x1step,x1i,nx2,x2grid,x2step,x2i,nx3,x3grid,x3step,x3i,nx4,x4grid,x4step,x4i,z)


location1 = min(nx1-1,max(1,floor((x1i-x1grid(1))/x1step)+1));
w1(2) = (x1i-x1grid(location1)) / (x1grid(location1+1)-x1grid(location1));
w1(1) = 1-w1(2);

location2 = min(nx2-1,max(1,floor((x2i-x2grid(1))/x2step)+1));
w2(2) = (x2i-x2grid(location2)) / (x2grid(location2+1)-x2grid(location2));
w2(1) = 1-w2(2);

location3 = min(nx3-1,max(1,floor((x3i-x3grid(1))/x3step)+1));
w3(2) = (x3i-x3grid(location3)) / (x3grid(location3+1)-x3grid(location3));
w3(1) = 1-w3(2);

location4 = min(nx4-1,max(1,floor((x4i-x4grid(1))/x4step)+1));
w4(2) = (x4i-x4grid(location4)) / (x4grid(location4+1)-x4grid(location4));
w4(1) = 1-w4(2);

zi = 0;

for m1 = 0:1
    for m2 = 0:1
        for m3 = 0:1
            for m4 = 0:1
            zi = zi + w1(m1+1)*w2(m2+1)*w3(m3+1)*w4(m4+1)*z(location1+m1,location2+m2,location3+m3,location4+m4);
            end
        end
    end
end
