#include <array>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

#include "body.h"

// Функция генерации случайного числа в диапазоне
double generate_random(double min, double max) {
    return min + static_cast<double>(std::rand()) / RAND_MAX * (max - min);
}

// Функция генерации случайного тела
Body generate_random_body() {
    double mass = generate_random(1e9, 1e10);
    Vector3D position({generate_random(0, 100.0), generate_random(0, 100.0),
                       generate_random(0, 100.0)});
    Vector3D velocity({generate_random(-10.0, 10.0),
                       generate_random(-10.0, 10.0),
                       generate_random(-10.0, 10.0)});
    return Body(mass, position, velocity);
}

// Функция для генерации вектора из n тел
std::vector<Body> generate_random_bodies(size_t n) {
    std::vector<Body> bodies;
    bodies.reserve(n); // Резервируем место в памяти

    for (size_t i = 0; i < n; ++i) {
        bodies.push_back(generate_random_body());
    }

    return bodies;
}