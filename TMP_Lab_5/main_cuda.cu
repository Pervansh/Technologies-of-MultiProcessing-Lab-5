#include <chrono>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// #define USE_DOUBLE_PREC // Использование двойной точности в вычислениях вместо одинарной

#include "body.h"
#include "file.h"
#include "random.h"
#include "runge_kutta_cuda.cu"

std::filesystem::path PATH_TO_RESULTS = "results/random_m1";

int main(int argc, char **argv) {
    double t_0 = 0;
    double t_n = 2;
    double tau = 0.1;
    int steps = static_cast<int>((t_n - t_0) / tau);

    std::srand(static_cast<unsigned>(
        std::time(nullptr))); // Инициализация генератора случайных чисел

    // Список чисел тел, для которых выполняется программа
    std::vector<size_t> numbers_of_bodies = {100, 200, 300};

    // Цикл по количеству тел
    for (size_t n : numbers_of_bodies) {
        // Генерация случайных тел
        std::vector<Body> bodies = generate_random_bodies(n);

        size_t total_bodies = bodies.size();

        // Запуск метода Рунге-Кутты
        const double time_start = MPI_Wtime();
        runge_kutta_2_parallel_improved(
            bodies, bodies_local, t_0, t_n, tau, counts, offsets,
            PATH_TO_RESULTS / ("random_test_" + std::to_string(size) + "_" +
                               std::to_string(n) + ".txt"));
        const double time = MPI_Wtime() - time_start;

        if (rank == 0) {
            printf("processes: %d\nn: %zu\ntime: %f\nstep time: %f\n", size, n,
                   time, time / steps);
        }
    }

    // Освобождение ресурсов
    free_mpi_types();
    MPI_Finalize();

    return 0;
}
