#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <array>

#ifdef USE_DOUBLE_PREC
    using numType = double;
#else
    using numType = float;
#endif

const numType G = 6.67e-11;
const numType EPS = 1e-9;

// Вектор в трехмерном пространстве
struct Vector3D {
    std::array<numType, 3> vector;

    // Конструкторы
    Vector3D() : vector{0.0, 0.0, 0.0} {} // Инициализируем по умолчанию с нулями

    Vector3D(const std::array<numType, 3>& other) : vector(other) {}

    Vector3D(const Vector3D& other) : vector(other.vector) {}

    // Метод для получения массива
    const std::array<numType, 3>& getArray() const {
        return vector;
    }

    // Перегрузка оператора умножения 
    Vector3D operator*(numType scalar) const {
        std::array<numType, 3> new_vector = {scalar * vector[0], scalar * vector[1], scalar * vector[2]};
        return Vector3D(new_vector);
    }

    // Перегрузка оператора умножения 
    Vector3D operator*(Vector3D& other) const {
        std::array<numType, 3> new_vector = {other.vector[0] * vector[0], other.vector[1] * vector[1], other.vector[2] * vector[2]};
        return Vector3D(new_vector);
    }

    // Перегрузка оператора сложения
    Vector3D operator+(const Vector3D& other) const {
        std::array<numType, 3> new_vector = {vector[0] + other.vector[0], vector[1] + other.vector[1], vector[2] + other.vector[2]};
        return Vector3D(new_vector);
    }

    // Перегрузка оператора вычетания
    Vector3D operator-(const Vector3D& other) const {
        std::array<numType, 3> new_vector = {vector[0] - other.vector[0], vector[1] - other.vector[1], vector[2] - other.vector[2]};
        return Vector3D(new_vector);
    }

    // Перегрузка оператора +=
    Vector3D& operator+=(const Vector3D& other) {
        vector[0] += other.vector[0];
        vector[1] += other.vector[1];
        vector[2] += other.vector[2];
        return *this;
    }

    // Перегрузка оператора -=
    Vector3D& operator-=(const Vector3D& other) {
        vector[0] -= other.vector[0];
        vector[1] -= other.vector[1];
        vector[2] -= other.vector[2];
        return *this;
    }

    Vector3D& operator=(const std::initializer_list<numType>& list) {
        auto it = list.begin();
        vector[0] = *it++;
        vector[1] = *it++;
        vector[2] = *it;
        return *this;
    }

    // Норма вектора
    double norm() const {
        return std::sqrt(vector[0] * vector[0] +
                         vector[1] * vector[1] +
                         vector[2] * vector[2]);
    }
};


// Структура тела
struct Body {
    double mass;
    Vector3D radius_vector, velocity;

    Body() : mass(0.0), radius_vector(), velocity() {}
    
    Body(numType mass_, const std::array<numType, 3>& pos, const std::array<numType, 3>& vel) : mass(mass_), radius_vector(pos), velocity(vel) {}

    Body(double mass_, const Vector3D& pos, const Vector3D& vel)
        : mass(mass_), radius_vector(pos), velocity(vel) {}

    Body(const Body& other)
        : mass(other.mass),
          radius_vector(other.radius_vector),
          velocity(other.velocity) {}

    Body& operator=(const Body& other) {
        mass = other.mass;
        radius_vector = other.radius_vector;
        velocity = other.velocity;

        return *this;
    }
};
