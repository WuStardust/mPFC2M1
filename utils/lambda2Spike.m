% generate spike as a Bernoulli random variable
function spike = lambda2Spike(lambda)
    randSeq = rand(size(lambda));
    spike = double(randSeq < lambda);
end
