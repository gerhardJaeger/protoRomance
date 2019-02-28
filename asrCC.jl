using RCall;
using DataFrames;
using CSV;
using Statistics;
using StatsFuns;
using Random;
Random.seed!(12345);
R"library(ape)";

const N = 1000;

const posterior = CSV.read("output/ccRates.log", delim='\t');
const treetrace = CSV.read("romance.posterior.tree", delim='\t')[1];
const data = CSV.read("romanceMultiMtx.csv")[:,2:end];
const taxa = CSV.read("romanceMultiMtx.csv")[:,1];
const states = "abcdefghijk";


#__convention__
#- nodes are numbered in postorder:
#   - 1..n: tips
#   - n+1..2n-1: internal nodes, mother before daughters

const k = 11;
q = ones((k, k));
for i in 1:k
    q[i, i] = 1-k
end


asr = zeros((N, 40, k));



for g in 1:N
    treeIndex = rand(1:size(treetrace)[1]);
    rateIndex = rand(1:size(posterior)[1]);
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
    for c in 1:40
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

bestChars = mapslices(argmax, mapslices(mean, asr, dims=1)[1, :, :], dims=2) .- 1;

concepts = map(string, names(data));

ccASR = [concepts[c] * ":" * string(bestChars[c])
         for c in 1:40];

DataFrame(concept=concepts, cc = ccASR) |> CSV.write("asrCC.csv");
