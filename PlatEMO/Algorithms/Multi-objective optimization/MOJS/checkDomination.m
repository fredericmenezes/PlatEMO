%% This function checks the domination inside the population.
function domi_vector = checkDomination(fitness)
    Np = size(fitness,1);
    if Np>2
        domi_vector = zeros(Np,1);
        all_perm = nchoosek(1:Np,2);    % Possible permutations
        all_perm = [all_perm; [all_perm(:,2) all_perm(:,1)]];

        d = dominates(fitness(all_perm(:,1),:),fitness(all_perm(:,2),:));
        dominated_particles = unique(all_perm(d==1,2));
        domi_vector(dominated_particles) = 1;
    else
        domi_vector=ones(Np,1);
    end
end