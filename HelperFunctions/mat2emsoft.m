function mat = mat2emsoft(Material)
    if any(strcmpi(Material, {'ni', 'nickel'}))
        mat = 'nickel';
    else
        mat = Material;
    end
end