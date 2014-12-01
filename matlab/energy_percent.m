function energy_ratio = energy_percent(coeff)
    coeff_count = size(coeff, 1);
    total_energy = coeff'*coeff;
    energy_ratio = zeros(coeff_count, 1);
    cur_energy = 0;
    for i = 1:coeff_count
        cur_energy = cur_energy + coeff(i,1)*coeff(i,1);
        energy_ratio(i) = cur_energy / total_energy;
    end
end