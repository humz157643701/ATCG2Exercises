/**************************************************************************
 * Copyright 2019 University of Bonn
 *
 * authors:
 *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
 *  - Lukas Bode <lbode@cs.uni-bonn.de>
 *
 * file creation date: 2019-10-20
 *
 *************************************************************************

 Read X-Rite's AxF material format, which is based on HDF5.

 This is a custom and unofficial implementation!
 For better reliability, please refer to X-Rite's official SDK:
 https://www.xrite.com/axf

*/

#include <cstdint>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

namespace AxfReader
{

// simple struct storing texture dimensions
struct dims_t
{
    uint32_t height;
    uint32_t width;
    uint32_t num_channels;
};

// auxiliary data types for storing textures
typedef std::map<std::string, std::vector<float> > TextureType;
typedef std::map<std::string, dims_t > TextureDimType;

// simple data type for an AxF resource definition
typedef std::pair<std::string, std::string> AxFResourceType;

// read an array from HDF5
template<typename T> bool
read_hdf5_array(HighFive::DataSet& dataset, std::vector<T>& data, dims_t& dims)
{
    std::vector<size_t> dims_hdf5 = dataset.getSpace().getDimensions();
    if (dataset.getDataType() != HighFive::AtomicType<T>())
    {
        std::cerr << "AxfReader: wrong data type" << std::endl;
        return false;
    }
    size_t num_elements = 1;
    for (auto dim: dims_hdf5)
    {
        num_elements *= dim;
    }
    data.resize(num_elements);
    dataset.read(&data[0]);

    // handle single value cases like uniform fresnel
    if (dims_hdf5.size() == 1)
    {
        dims.height = dims_hdf5[0];
        dims.width = 1;
        dims.num_channels = 1;
    }
    else
    {
        dims.height = dims_hdf5[0];
        dims.width = dims_hdf5[1];
        dims.num_channels = dims_hdf5[2];
    }

    return true;
}

// read AxF from file name, in case multiple materials are contained, material_index specifies which one to load
bool readAxF(std::string filename,
             TextureType& textures,
             TextureDimType& texture_dims,
             uint32_t material_index = 0)
{

    // AxF resource definitions: texture name and corresponding path in AxF
    // inspect other AxF definitions using hdfview under /com.xrite.Materials/material_name/com.xrite.Resources/
    std::vector<std::pair<std::string, std::string> > axf_defs;
    axf_defs.clear();
    axf_defs.push_back(AxFResourceType("diffuse", "DiffuseModel/Color"));
    axf_defs.push_back(AxFResourceType("normal", "DiffuseModel/Normal"));
    axf_defs.push_back(AxFResourceType("specular", "SpecularModel/Color"));
    axf_defs.push_back(AxFResourceType("lobes", "SpecularModel/Lobes"));
    axf_defs.push_back(AxFResourceType("aniso", "SpecularModel/AnisotropicRotation"));
    axf_defs.push_back(AxFResourceType("fresnel", "SpecularModel/Fresnel"));
    axf_defs.push_back(AxFResourceType("displacement", "DisplacementFilter/Height"));
    axf_defs.push_back(AxFResourceType("transparency", "TransparencyFilter/Alpha"));
    //axf_defs.push_back(AxFResourceType("clearcoat_ior", "ClearCoat/IOR"));
    //axf_defs.push_back(AxFResourceType("clearcoat_normal", "ClearCoat/Normal"));
    // ...

    // dimensions of material
    std::string material_name;
    try
    {
        HighFive::File file(filename, HighFive::File::ReadOnly);

        HighFive::Group group;
        group = file.getGroup("/com.xrite.Materials");
        std::vector<std::string> material_list = group.listObjectNames();
        material_name = material_list[material_index];
        std::cout << "AxfReader: read materials: " << material_name << std::endl;
        group = file.getGroup("/com.xrite.Materials/" + material_name + "/com.xrite.Resources/");
        std::vector<std::string> resources_list = group.listObjectNames();

        uint32_t version_major, version_minor, version_revision;
        file.getGroup("/").getAttribute("axf.version.major").read(version_major);
        file.getGroup("/").getAttribute("axf.version.minor").read(version_minor);
        file.getGroup("/").getAttribute("axf.version.revision").read(version_revision);

        for (size_t ri = 0; ri < resources_list.size(); ri++)
        {
            // iterates over resource folders like DiffuseModel, DisplacementFilter, ...
            std::string resource_base = resources_list[ri].substr(0, resources_list[ri].find_first_of("/"));
            for (size_t ii = 0; ii < axf_defs.size(); ii++)
            {
                std::string object_path = axf_defs[ii].second;
                // read all resources under the corresponding folder
                if (object_path.size() >= resource_base.size() &&
                    object_path.substr(0, resource_base.size()) == resource_base)
                {
                    std::string object = object_path.substr(resource_base.size() + 1);
                    std::vector<std::string> objects = file.getGroup("/com.xrite.Materials/" +
                        material_name + "/com.xrite.Resources/" + resource_base).listObjectNames();
                    if (std::find(objects.begin(), objects.end(), object) == objects.end())
                    {
                        continue;
                    }
                    std::string texture_name = axf_defs[ii].first;
                    std::string dataset_path = "/com.xrite.Materials/" +
                            material_name + "/com.xrite.Resources/" +
                            resource_base + "/" + object + "/Data";
                    HighFive::DataSet dataset = file.getDataSet(dataset_path);
                    std::vector<float> pixels;
                    dims_t dims;
                    if (!read_hdf5_array(dataset, pixels, dims))
                    {
                        std::cerr << "AxfReader: error reading " << dataset_path << std::endl;
                        return false;
                    }
                    textures[texture_name] = pixels;
                    texture_dims[texture_name] = dims;
                    continue;
                }
            }
        }
    }
    catch (std::exception err)
    {
        std::cerr << "AxfReader: " << err.what() << std::endl;
        return false;
    }

    // duplicate lobe channel for isotropic SVBRDFs so we can use the same anisotropic model in the shader
    if (textures.find("lobes") != textures.end() && texture_dims["lobes"].num_channels == 1)
    {
        std::vector<float> pix_orig = textures["lobes"];
        std::vector<float> lobe_pixels(pix_orig.size() * 2);
        std::copy(pix_orig.begin(), pix_orig.end(), &lobe_pixels[0]);
        std::copy(pix_orig.begin(), pix_orig.end(), &lobe_pixels[pix_orig.size()]);
        textures["lobes"] = lobe_pixels;
        texture_dims["lobes"].num_channels = 2;
    }

    // also create a map for anisotropic rotation (initialized with 0s)
    if (textures.find("aniso") == textures.end())
    {
        textures["aniso"] = std::vector<float>(texture_dims["diffuse"].width * texture_dims["diffuse"].height, 0.);
        texture_dims["aniso"].width = texture_dims["diffuse"].width;
        texture_dims["aniso"].height = texture_dims["diffuse"].height;
        texture_dims["aniso"].num_channels = 1;
    }
    std::cout << "AxfReader: done" << std::endl;
    return true;
}

}
