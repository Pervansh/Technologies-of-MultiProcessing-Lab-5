#pragma once
#include <cstdio>
#include <string>
#include <vector>

#include <crt/device_functions.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "body.h"
#include "file.h"

__device__ void get_acceleration(const Body *bodies, int& n, const int& cudaBSshared, Vector3D *accelerations) {
    Vector3D difference;
    Vector3D acceleration;
    numType denom, difference_norm;

    int i = blockIdx.x * blockDim.x + threadIdx.x;

    int localI = threadIdx.x;

    auto radius_vector_i = bodies[i].radius_vector;

    if (i < n) {
        extern __shared__ Vector3D sharedPoses[];
        extern __shared__ numType  sharedMasses[];

        numType res = 0;
        for (int bk = 0; bk < n; bk += cudaBSshared) {
            sharedPoses[localI]  = bodies[bk + localI].radius_vector;
            sharedMasses[localI] = bodies[bk + localI].mass;

            __syncthreads(); 

            for (int localK = 0; localK < cudaBSshared; ++localK) {
                difference = sharedPoses[localK] - radius_vector_i;
                difference_norm = difference.norm();
                denom = difference_norm * difference_norm * difference_norm;
                denom = denom > EPS ? denom : EPS;
                acceleration += difference * sharedMasses[localK] * (1. / denom);
            }
                // res += sharedPoses[localI * (cudaBSshared + 1) + localK] * b[localK * cudaBSshared + localJ];

            __syncthreads();
        }

        accelerations[i] = acceleration;
    }

    for (const auto &body : bodies) {
        difference = body.radius_vector - radius_vector_i;
        difference_norm = difference.norm();
        denom = difference_norm * difference_norm * difference_norm;
        denom = denom > EPS ? denom : EPS;
        acceleration += difference * body.mass * (1. / denom);
    }

    acceleration = acceleration * G;
}

// Параллельный метод Рунге-Кутты второго порядка
__host__ void runge_kutta_2_parallel(std::vector<Body> bodies,
                            std::vector<Body> bodies_local, double t_0,
                            double t_n, double tau,
                            const std::vector<int> &counts,
                            const std::vector<int> &offsets,
                            const std::string &filename) {
    std::ofstream file(filename);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int local_size = counts[rank];

    std::vector<Body> bodies_local_buffer(bodies_local);
    std::vector<Vector3D> acceleration_prev;
    Vector3D acceleration_buffer, acceleration;

    for (double t = t_0; t < t_n + tau / 2.; t += tau) {
        // Процесс 0 записывает результаты в файл
        if (rank == 0) {
            for (size_t i = 0; i < bodies.size(); ++i) {
                write_step_to_file(t, bodies[i].radius_vector, filename);
            }
        }

        // Первый шаг метода Рунге-Кутты
        // printf("[INFO] Make first step from t = %f\n", t);
        for (size_t i = 0; i < local_size; ++i) {
            // printf("[INFO] Rank %d: \n", rank);
            get_acceleration(bodies, bodies_local[i].radius_vector,
                             acceleration_buffer);
            // printf("[INFO] Acceleration = {%f, %f, %f}\n",
            // acceleration_buffer.vector[0], acceleration_buffer.vector[1],
            // acceleration_buffer.vector[2]);

            bodies_local_buffer[i].radius_vector +=
                bodies_local[i].velocity * tau * 0.5;
            bodies_local_buffer[i].velocity += acceleration_buffer * tau * 0.5;

            // printf("[INFO] Rank %d: Body %zu: r = {%f, %f, %f}\n", rank, i +
            // 1, bodies_local_buffer[i].radius_vector.vector[0],
            // bodies_local_buffer[i].radius_vector.vector[1],
            // bodies_local_buffer[i].radius_vector.vector[2]);
        }

        // Синхронизируем данные между процессами
        MPI_Allgatherv(bodies_local_buffer.data(), local_size, mpi_body_type,
                       bodies.data(), counts.data(), offsets.data(),
                       mpi_body_type, MPI_COMM_WORLD);

        // printf("[INFO] Make second step from t = %f\n", t);
        for (size_t i = 0; i < local_size; ++i) {
            // printf("[INFO] Rank %d: \n", rank);
            get_acceleration(bodies, bodies_local_buffer[i].radius_vector,
                             acceleration);
            // printf("[INFO] Acceleration = {%f, %f, %f}\n",
            // acceleration.vector[0], acceleration.vector[1],
            // acceleration.vector[2]);
            bodies_local[i].radius_vector +=
                bodies_local_buffer[i].velocity * tau;
            bodies_local[i].velocity += acceleration * tau;

            // printf("[INFO] Rank %d: Body %zu: r = {%f, %f, %f}\n", rank, i +
            // 1, bodies_local[i].radius_vector.vector[0],
            // bodies_local[i].radius_vector.vector[1],
            // bodies_local[i].radius_vector.vector[2]);
        }

        // Синхронизируем данные между процессами
        MPI_Allgatherv(bodies_local.data(), local_size, mpi_body_type,
                       bodies.data(), counts.data(), offsets.data(),
                       mpi_body_type, MPI_COMM_WORLD);

        bodies_local_buffer = bodies_local;
    }

    // printf("[INFO] Rank %d: Finished solving system\n", rank);
}

// Параллельный метод Рунге-Кутты второго порядка
__host__ void runge_kutta_2_parallel_improved(std::vector<Body> bodies,
                                     std::vector<Body> bodies_local, double t_0,
                                     double t_n, double tau,
                                     const std::string &filename) {
    std::ofstream file(filename);

    std::vector<Body> bodies_local_buffer(bodies_local);
    std::vector<Vector3D> acceleration_prev;
    Vector3D acceleration_buffer, acceleration;

    for (double t = t_0; t < t_n + tau / 2.; t += tau) {
        // Процесс 0 записывает результаты в файл
        for (size_t i = 0; i < bodies.size(); ++i) {
            write_step_to_file(t, bodies[i].radius_vector, filename);
        }

        // Первый шаг метода Рунге-Кутты
        // printf("[INFO] Make first step from t = %f\n", t);
        for (size_t i = 0; i < local_size; ++i) {
            // printf("[INFO] Rank %d: \n", rank);
            get_acceleration(bodies, bodies_local[i].radius_vector,
                             acceleration_buffer);
            // printf("[INFO] Acceleration = {%f, %f, %f}\n",
            // acceleration_buffer.vector[0], acceleration_buffer.vector[1],
            // acceleration_buffer.vector[2]);

            bodies_local_buffer[i].radius_vector +=
                bodies_local[i].velocity * tau * 0.5;
            bodies_local_buffer[i].velocity += acceleration_buffer * tau * 0.5;

            // printf("[INFO] Rank %d: Body %zu: r = {%f, %f, %f}\n", rank, i +
            // 1, bodies_local_buffer[i].radius_vector.vector[0],
            // bodies_local_buffer[i].radius_vector.vector[1],
            // bodies_local_buffer[i].radius_vector.vector[2]);
        }

        // Синхронизируем данные между процессами
        MPI_Allgatherv(bodies_local_buffer.data(), local_size,
                       mpi_body_radius_vector_type_part, bodies.data(),
                       counts.data(), offsets.data(),
                       mpi_body_radius_vector_type_part, MPI_COMM_WORLD);

        // printf("[INFO] Make second step from t = %f\n", t);
        for (size_t i = 0; i < local_size; ++i) {
            // printf("[INFO] Rank %d: \n", rank);
            get_acceleration(bodies, bodies_local_buffer[i].radius_vector,
                             acceleration);
            // printf("[INFO] Acceleration = {%f, %f, %f}\n",
            // acceleration.vector[0], acceleration.vector[1],
            // acceleration.vector[2]);
            bodies_local[i].radius_vector +=
                bodies_local_buffer[i].velocity * tau;
            bodies_local[i].velocity += acceleration * tau;

            // printf("[INFO] Rank %d: Body %zu: r = {%f, %f, %f}\n", rank, i +
            // 1, bodies_local[i].radius_vector.vector[0],
            // bodies_local[i].radius_vector.vector[1],
            // bodies_local[i].radius_vector.vector[2]);
        }

        // Синхронизируем данные между процессами
        MPI_Allgatherv(bodies_local.data(), local_size,
                       mpi_body_radius_vector_type_part, bodies.data(),
                       counts.data(), offsets.data(),
                       mpi_body_radius_vector_type_part, MPI_COMM_WORLD);

        bodies_local_buffer = bodies_local;
    }

    // printf("[INFO] Rank %d: Finished solving system\n", rank);
}
