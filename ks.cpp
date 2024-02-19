#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <random>

#include <vector>
double ksTestStatistic(const std::vector<double>& sample1, const std::vector<double>& sample2) {
    // Combine and sort the two samples
    std::vector<double> combined = sample1;
    combined.insert(combined.end(), sample2.begin(), sample2.end());
    std::sort(combined.begin(), combined.end());

    // Calculate the empirical distribution functions
    const int totalSize = combined.size();
    std::vector<double> cdf1(totalSize), cdf2(totalSize);

    for (int i = 0; i < totalSize; ++i) {
        cdf1[i] = static_cast<double>(i + 1) / totalSize;
        cdf2[i] = static_cast<double>(std::lower_bound(combined.begin(), combined.end(), combined[i]) - combined.begin()) / totalSize;
    }

    // Calculate the K-S test statistic
    double ksStatistic = 0.0;
    for (int i = 0; i < totalSize; ++i) {
        double diff1 = std::abs(cdf1[i] - cdf2[i]);
        double diff2 = std::abs(cdf2[i] - cdf1[i]);
        ksStatistic = std::max(ksStatistic, std::max(diff1, diff2));
    }

    return ksStatistic;
}

double ksTestPValue(double ksStatistic, int n, int m) {
    // Use the K-S distribution or approximation to calculate the p-value
    // In this example, we use a simple simulation approach for illustration

    const int numSimulations = 10000;
    int countExceed = 0;

    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    for (int i = 0; i < numSimulations; ++i) {
        std::vector<double> simulatedSample1(n);
        std::vector<double> simulatedSample2(m);

        // Generate random samples
        for (double& value : simulatedSample1) {
            value = distribution(generator);
        }

        for (double& value : simulatedSample2) {
            value = distribution(generator);
        }

        // Calculate the K-S test statistic for the simulated samples
        double simulatedKS = ksTestStatistic(simulatedSample1, simulatedSample2);

        // Check if the simulated test statistic exceeds the observed test statistic
        if (simulatedKS >= ksStatistic) {
            countExceed++;
        }
    }

    // Calculate the p-value based on the simulations
    double pValue = static_cast<double>(countExceed) / numSimulations;

    return pValue;
}

int main() {
    // Example data (replace with your own samples)
    std::random_device rd;
    std::default_random_engine generator(rd());
    double mean1 = 0.0;
    double stddev1 = 1.0;
    std::normal_distribution<double> distribution1(mean1, stddev1);

    double mean2 = 0.0;
    double stddev2 = 1.0;
    std::normal_distribution<double> distribution2(mean1, stddev1);

    int n=1000;
    std::vector<double> samples1;

    for (int i = 0; i < n; ++i) {
        double sample = distribution1(generator);
        samples1.push_back(sample);
    }

    std::vector<double> samples2;

    for (int i = 0; i < n; ++i) {
        double sample = distribution2(generator);
        samples2.push_back(sample);
    }




    // Perform the K-S test
    double ksStatistic = ksTestStatistic(samples1, samples2);

    // Calculate the p-value
    double pValue = ksTestPValue(ksStatistic, samples1.size(), samples2.size());


    // Output the results
    std::cout << "K-S Test Statistic: " << ksStatistic << std::endl;
    std::cout << "P-value: " << pValue << std::endl;

    return 0;
}

