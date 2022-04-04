function kruh = gen_circle(r)
    [y,x] = meshgrid(1:r*2+1,1:r*2+1);
    kruh = zeros(size(x));
    s = r+1;
    kruh(sqrt( (x - s).^2  + (y - s).^2 ) < (r+0.5) & sqrt((x - s).^2  + (y - s).^2 ) > (r-0.5) ) = 1;
end