#include <iostream>
#include <cstdlib>

int main()
{
    setlocale(LC_ALL, "rus");
    std::cout << "Запускается визуализация..." << std::endl;
    system("python visualize.py");

    return 0;
}