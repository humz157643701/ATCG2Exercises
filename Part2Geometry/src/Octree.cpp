//
// Created by alex on 12.12.19.
//

#include <Octree.h>

//namespace students {

size_t clz(uint64_t x) {
    int n = 64;
    unsigned y;

    y = x >> 32;
    if (y != 0) {
        n = n - 32;
        x = y;
    }
    y = x >> 16;
    if (y != 0) {
        n = n - 16;
        x = y;
    }
    y = x >> 8;
    if (y != 0) {
        n = n - 8;
        x = y;
    }
    y = x >> 4;
    if (y != 0) {
        n = n - 4;
        x = y;
    }
    y = x >> 2;
    if (y != 0) {
        n = n - 2;
        x = y;
    }
    y = x >> 1;
    if (y != 0) return n - 2;
    return n - x;
}

//----------------------------------------------------------------------------------------------------------------------

void Octree::clear() {
    indices.clear();
    OctreeBase::clear();
}

void Octree::build(size_t leaf_size) {
    init(leaf_size);
    create_bfs();
    storage.shrink_to_fit();
}

Octree::ResultSet Octree::query_radius(const Vec3 &point, float radius, size_t depth) const {
    ResultSet resultSet;
    query_radius(0, root_aabb, point, radius, resultSet, depth);
    return resultSet;
}

Octree::ResultSet Octree::query_knn(const Vec3 &point, int k, size_t depth) const {
    ResultSet resultSet(k);
    query_knn(0, root_aabb, point, resultSet, depth);
    return resultSet;
}

void Octree::init(size_t leaf_size) {
    this->leaf_size = leaf_size;

    size_t N = points.size();
    indices = std::vector<size_t>(N);
    root_aabb = AABB();

    for (size_t i = 0; i < N; ++i) {
        root_aabb.grow(points[i]);
        indices[i] = i;
    }
}

void Octree::create_bfs() {
    storage.push_back(Node(0, 1, 0, 0, OctantIndices(0, points.size() - 1, points.size())));

    std::queue<AABB> queue;
    queue.emplace(root_aabb);

    for (size_t i = 0; i < storage.size() && !queue.empty(); ++i) {
        AABB &aabb = queue.front();
        depth = std::max(tree_depth(storage[i].loc_code), depth);

        if (storage[i].data.size > leaf_size && tree_depth(storage[i].loc_code) < 21) {
            std::vector<size_t> buckets[8];

            for (size_t j = 0; j < storage[i].data.size; ++j) {
                size_t idx = indices[storage[i].data.start + j];
                uint8_t mortonCode = morton_code(points[idx], aabb.center);
                buckets[mortonCode].push_back(idx);
            }

            bool firsttime = true;
            size_t child_start = storage[i].data.start;

            Vec3 child_extent = aabb.halfsize / 2;
            for (uint8_t octant = 0; octant < 8; ++octant) {
                if (buckets[octant].empty()) continue;
                if (firsttime) {
                    firsttime = false;
                    storage[i].first_child_index = storage.size();
                }

                set_child_exists(storage[i].config, octant);


                size_t child_size = buckets[octant].size();
                size_t child_end = child_start + child_size - 1;
                storage.push_back(Node(0, child_loc_code(storage[i].loc_code, octant), 0, i,
                                       OctantIndices(child_start, child_end, child_size)));
                for (size_t l = 0; l < child_size; ++l) {
                    indices[child_start + l] = buckets[octant][l];
                }
                buckets[octant].clear();
                child_start += child_size;
                queue.emplace(child_box(aabb, octant, child_extent));
            }
        }
        queue.pop();
    }
}

void
Octree::query_radius(size_t index, const AABB &aabb, const Vec3 &point, float radius, ResultSet &resultSet,
                     size_t depth) const {
    const Node &node = storage[index];

    if (aabb.is_inside_sphere(point, radius)) {
        for (size_t i = node.data.start; i < node.data.end; ++i) {
            resultSet.add_point(indices[i]);
        }
        return;
    }

    if (node.config == 0 || tree_depth(node.loc_code) == depth) {
        for (size_t i = node.data.start; i < node.data.end; ++i) {
            size_t idx = indices[i];
            float dist = (points[idx] - point).length();
            if (dist < radius) resultSet.add_point(idx);
        }
        return;
    }

    Vec3 child_extent = aabb.halfsize / 2;
    uint8_t offset = 0;
    for (uint8_t i = 0; i < 8; ++i) {
        if (!child_exists(node.config, i)) continue;
        AABB childBox = child_box(aabb, i, child_extent);
        if (!childBox.intersect(point, radius)) continue;
        query_radius(node.first_child_index + offset, childBox, point, radius, resultSet, depth);
        ++offset;
    }
}

void
Octree::query_knn(size_t index, const AABB &aabb, const Vec3 &point, ResultSet &resultSet, size_t depth) const {
    //descend to leaf
    //first order children by distance to query point

    const Node &node = storage[index];
    //collect up to k points
    if (node.config == 0 || tree_depth(node.loc_code) == depth) {
        if (node.data.size > 0) {
            for (size_t i = node.data.start; i <= node.data.end; ++i) {
                size_t idx = indices[i];
                resultSet.add_point(idx, (point - points[idx]).length());
            }
            resultSet.sort();
            return;
        }
    }

    uint8_t morton = morton_code(point, aabb.center);
    Vec3 child_extent = aabb.halfsize / 2;
    uint8_t offset = 0;
    for (uint8_t i = 0; i < 8; ++i) {
        if (morton == i) break;
        if (!child_exists(node.config, i)) continue;
        ++offset;
    }
    if (child_exists(node.config, morton)) {
        query_knn(node.first_child_index + offset, child_box(aabb, morton, child_extent), point, resultSet, depth);
    }

    offset = 0;
    for (uint8_t i = 0; i < 8; ++i) {
        if (!child_exists(node.config, i)) continue;
        if (morton == i) {
            ++offset;
            continue;
        };
        AABB childBox = child_box(aabb, i, child_extent);

        if (!resultSet.is_full() || childBox.intersect(point, resultSet.worst_dist)) {
            query_knn(node.first_child_index + offset, childBox, point, resultSet, depth);
        }
        ++offset;
    }
}

//}