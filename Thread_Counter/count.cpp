#include <execution>
#include <vector>
#include <atomic>
#include <thread>
#include <iostream>

int main() {
    std::vector<int> data(10'000'000, 1);
    std::atomic<int> thread_count{0};
    thread_local bool counted = false;

    std::for_each(std::execution::par_unseq, data.begin(), data.end(), [&](int& x) {
        if (!counted) {
            counted = true;
            thread_count++;
        }
        x += 1;
    });

    std::cout << "Estimated threads used: " << thread_count.load() << "\n";
}