#include <cstdint>
#include <cstddef>
#include <cmath>
#include <cstdio>
#include <unordered_map>

extern "C" double computeLosslessEntropy_32bits_fast(void* data, size_t nbEle)
{
    std::unordered_map<std::uint32_t, size_t> count;
    const std::uint32_t* values = static_cast<const std::uint32_t*>(data);

    size_t totalSteps = nbEle < 1000 ? nbEle : 1000;
    size_t stepSize = nbEle / totalSteps;
    if (stepSize == 0) stepSize = 1;

    for (size_t i = 0; i < nbEle; i++) {
        count[values[i]]++;
        if ((i + 1) % stepSize == 0) {
            size_t step = (i + 1) / stepSize;
            if (step <= totalSteps)
                fprintf(stderr, "\rcomputing 32bit entropy: %zu/%zu", step, totalSteps);
        }
    }
    fprintf(stderr, "\rcomputing 32bit entropy: %zu/%zu\n", totalSteps, totalSteps);

    double h = 0;
    for (auto it = count.begin(); it != count.end(); ++it)
        h += it->second * std::log2(static_cast<double>(it->second));
    h /= nbEle;
    h = std::log2(static_cast<double>(nbEle)) - h;

    return h;
}
