function tab = tmppp()

tab = zeros(4096 * 4, 1);
for c = 0 : 31
    for r = 0 : 31
        i = c * 32 + r; n1 = i * 4;
        if c < 31, tab(n1 + 0 + 1) = i + 32; else tab(n1 + 0 + 1) = i; end
        if r < 31, tab(n1 + 1 + 1) = i + 1; else tab(n1 + 1 + 1) = i; end
        if c < 31 && r < 31, tab(n1 + 2 + 1) = i + 32 + 1; else tab(n1 + 2 + 1) = i; end
        if c > 0 && r < 31, tab(n1 + 3 + 1) = i - 32 + 1; else tab(n1 + 3 + 1) = i; end
    end
end

end