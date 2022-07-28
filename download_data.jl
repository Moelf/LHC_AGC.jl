using JSON3

const N = 100

const a = JSON3.read(read("./ntuples.json"));

for key in keys(a)
    Threads.@threads for n in first(a[key][:nominal][:files], N)
        url = n[:path]
        fname = last(split(url, '/'))
        out = joinpath("/data/jiling/Analysis_Grand_Challenge/", fname)
        isfile(out) && continue
        download(url, out)
    end
end
