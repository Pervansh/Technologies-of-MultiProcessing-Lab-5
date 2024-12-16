#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include "body.h"

// Функция для чтения данных из файла
std::vector<Body> read_from_file(const std::string& input_file) {
    std::vector<Body> bodies; 
    
    // Сообщение о начале чтения файла
    printf("[INFO] Start reading file %s to initialize system\n", input_file.c_str());

    // Открытие файла
    std::ifstream file(input_file);
    if (!file.is_open()) {
        std::cerr << "[ERROR] Failed to open file " << input_file << std::endl;
        throw std::runtime_error("File opening failed");
    }

    // Чтение количества тел
    int n;
    file >> n;

    // Чтение данных тел
    double x, y, z, vx, vy, vz, mass;
    for (int i = 0; i < n; ++i) {
        file >> mass >> x >> y >> z >> vx >> vy >> vz;

        Body body;
        body.radius_vector = Vector3D({x, y, z});
        body.velocity = Vector3D({vx, vy, vz});
        body.mass = mass;
        bodies.push_back(body);  // Добавление тела в систему
    }

    // Сообщение об успешном завершении чтения
    printf("[INFO] Successfully read file %s\n", input_file.c_str());
    file.close();

    // Возвращение вектора тел
    printf("[INFO] System size = %zu\n", bodies.size());
    return bodies;
}


// Дозаписывает состояние системыв файл
void write_step_to_file(double t, Vector3D& vector, const std::string& filename) {
    std::ofstream file(filename, std::ios::app);
    if (!file) {
        throw std::runtime_error("[ERROR] Failed to open file " + filename);
    }
    file << t << " " << vector.vector[0] << " " << vector.vector[1] << " " << vector.vector[2] << std::endl;
    file.close();
}
