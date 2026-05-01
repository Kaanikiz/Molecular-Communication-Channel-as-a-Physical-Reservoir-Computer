function xyz = randomPointsInSphere(numPts, radius)

    v = randn(numPts, 3);
    v = v ./ vecnorm(v, 2, 2);

    r = radius * rand(numPts, 1).^(1/3);
    xyz = v .* r;
end
