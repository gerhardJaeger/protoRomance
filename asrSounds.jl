using RCall;
using DataFrames;
using CSV;
using Statistics;
using StatsFuns;
using Random;
Random.seed!(12345);
R"library(ape)";

const N = 100;

const posteriorFull = CSV.read("output/soundRates.log", delim='\t');
const posterior = posteriorFull[posteriorFull.Iteration .>
                                div(maximum(posteriorFull.Iteration), 2), :];
const treetrace = CSV.read("romance.posterior.tree", delim='\t')[1];
const data = CSV.read("romanceAlignments.csv")[:,2:end];
const taxa = CSV.read("romanceMultiMtx.csv")[:,1];
const states = "358CELNSTZabcdefghijklmnopqrstuvwxyz0";

const nChars = size(data)[2];

#__convention__
#- nodes are numbered in postorder:
#   - 1..n: tips
#   - n+1..2n-1: internal nodes, mother before daughters

const k = 37;


asr = zeros((N, nChars, k));



for g in 1:N
    treeIndex = rand(1:size(treetrace)[1])
    rateIndex = rand(1:size(posterior)[1])
    rate1 = posterior[rateIndex, Symbol("nrates[1]")]
    rate2 = posterior[rateIndex, Symbol("nrates[2]")]
    rate3 = posterior[rateIndex, Symbol("nrates[3]")]
    q = zeros((k, k)) .+ rate3
    q[:, k] = repeat([rate1], k)
    q[k, :] = repeat([rate2], k)
    for i in 1:k
        q[i, i] -= sum(q[i, :])
    end
    treeS = treetrace[treeIndex];
    treeR = R"read.tree(text=$treeS)";
    tree = convert(Dict, treeR);
    nInternal = tree["Nnode"];
    n = length(tree["tip.label"]);
    nNodes = nInternal+n;
    branchLengths = zeros((nNodes, nNodes));
    for i in 1:length(tree["edge.length"])
        b = tree["edge"][i, :]
        l = tree["edge.length"][i]
        branchLengths[b[1], b[2]]=l
    end;
    for c in 1:nChars
        loglik = zeros((nNodes, k))
        rate = posterior[rateIndex, Symbol("siteRates[" * string(c) * "]")]
        for nd in 1:n
            nm = taxa[nd]
            s = data[taxa .== nm, c][1]
            if s != "-"
                si = findfirst(s, states)[1]
                loglik[nd, :] .-= Inf
                nm = tree["tip.label"][nd]
                loglik[nd, si] = 0
            end
        end
        for nd in nNodes:-1:(n+1)
            daughters = tree["edge"][tree["edge"][:, 1] .== nd, 2]
            ll = branchLengths[nd, daughters]
            for j in 1:length(daughters)
                d = daughters[j]
                loglik[nd, :] += mapslices(logsumexp,
                                           log.(rate * exp(ll[j] .* q')) .+
                                           loglik[d, :], dims=1)'
            end
        end
        logAsr = loglik[n+1, :]
        x = exp.(logAsr)
        x /= sum(x)
        asr[g, c, :] = x
    end
    print(g)
    print("\n")
end;

bestChars = mapslices(argmax, mapslices(mean, asr, dims=1)[1, :, :], dims=2);

characters = map(string, names(data));

soundASR = [states[c] for c in bestChars]


recon = DataFrame(concept=[split(c, ":")[1] for c in characters],
                  reconstruction = soundASR[:,1]);

concepts = unique(recon.concept);

protoRomance = DataFrame(concept=concepts,
                         reconstruction=[replace(join(recon.reconstruction[recon.concept .== c]), "0" => "") for c in concepts]) |> CSV.write("protoRomance.csv")


