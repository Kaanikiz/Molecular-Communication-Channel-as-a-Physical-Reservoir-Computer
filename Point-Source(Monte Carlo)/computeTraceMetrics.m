function metrics = computeTraceMetrics(b_ref, b_test, t)

    b_ref = b_ref(:);
    b_test = b_test(:);
    t = t(:);

    valid = isfinite(b_ref) & isfinite(b_test);
    b_ref = b_ref(valid);
    b_test = b_test(valid);
    t = t(valid);

    b_ref_centered = b_ref - mean(b_ref);
    b_test_centered = b_test - mean(b_test);

    metrics.corr = sum(b_ref_centered .* b_test_centered) / ...
        sqrt(sum(b_ref_centered.^2) * sum(b_test_centered.^2));

    metrics.nrmse = sqrt(mean((b_ref - b_test).^2)) / ...
        (max(b_ref) - min(b_ref));

    [~, idx_ref] = max(b_ref);
    [~, idx_test] = max(b_test);

    metrics.peakTimeError = abs(t(idx_ref) - t(idx_test));
end
