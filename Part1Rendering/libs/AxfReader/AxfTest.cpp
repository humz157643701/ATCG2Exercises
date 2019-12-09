#include "AxfReader.hpp"

#include <string>
#include <iostream>
#include <iomanip>


int main()
{
    std::cout << "start AxfReader test" << std::endl << std::endl;

    std::string axfFilename = "../happy_birthday.axf";

    // textures are stored in a map from texture names to float vectors holding the raw pixel values
    AxfReader::TextureType textures;

    // texture dimensions are stored separately
    AxfReader::TextureDimType texture_dims;

    if (!AxfReader::readAxF(axfFilename, textures, texture_dims))
    {
        std::cerr << "error reading AxF " << axfFilename << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << std::endl << "print shapes:" << std::endl;
    for (auto& pair : texture_dims)
    {
        std::cout << "key: " << std::setw(15) << std::left << pair.first
                  << " shape (h,w,c): (" << pair.second.height << ","
                                         << pair.second.width << ","
                                         << pair.second.num_channels << ")" << std::endl;
    }

    std::cout << std::endl << "print first 6 floats:" << std::endl;
    for (auto& pair : textures)
    {
        std::cout << "key: " << std::setw(15) << std::left << pair.first
                  << " values: " << pair.second[0];

        for (size_t i = 1; i < 6 && i < pair.second.size(); i++)
        {
            std::cout << ", " << pair.second[i];
        }
        std::cout << std::endl;
    }

    std::cout << std::endl << "test successful!" << std::endl;

    return EXIT_SUCCESS;
}
