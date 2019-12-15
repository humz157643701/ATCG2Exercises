//
// Created by alex on 12.12.19.
//

#ifndef STUDENTS_OCTREE_H
#define STUDENTS_OCTREE_H

/*
 * Use like this:
 * start by calling constructor:
 *
 *      Octree octree(points);
 *      build octree with e.g. leaf size = 5
 *      octree.build(5);
 *
 *  if you want radius query for all points in a sphere around p with radius r call:
 *
 *      ResultSet results = octree.query_radius(p, r);
 *
 *  result.idx_dist_pair contains the indices (and dist is 0 for all points due to speed).
 *
 *  if you want knn query for k nearest neghbors of p call:
 *
 *      ResultSet results = octree.query_knn(p, k);
 *
 *  result.idx_dist_pair contains the indices (and dist is the actual distance for all points to the query point).
 *
 * */

#include <vector>
#include <limits>
#include <functional>
#include <queue>
#include <cassert>
#include <cstdint>
#include <cmath>

//namespace students {

template<typename T>
inline size_t BIT(T pos) { return 1 << pos; }

template<typename T, typename U>
inline T DROP(T var, U pos) { return var >> pos; }

template<typename T, typename U>
inline bool CHECK_BIT(T var, U pos) { return DROP(var, pos) & 1; }

template<typename T, typename U>
inline void SET_BIT(T &var, U pos) { var |= BIT(pos); }

size_t clz(uint64_t var);


//----------------------------------------------------------------------------------------------------------------------

class Vec3 {
public:
    Vec3() : Vec3(Vec3::Zero()) {}

    Vec3(float x, float y, float z) : data{x, y, z} {}

    static Vec3 Zero() { return Vec3(0, 0, 0); }

    static Vec3 Lowest() {
        return Vec3(std::numeric_limits<float>::lowest(),
                    std::numeric_limits<float>::lowest(),
                    std::numeric_limits<float>::lowest());
    }

    float operator[](size_t i) const {
        return data[i];
    }

    float &operator[](size_t i) {
        return data[i];
    }

    Vec3 operator+(const Vec3 &o) const {
        return Vec3(data[0] + o[0], data[1] + o[1], data[2] + o[2]);
    }

    Vec3 operator-(const Vec3 &o) const {
        return Vec3(data[0] - o[0], data[1] - o[1], data[2] - o[2]);
    }

    Vec3 operator*(const Vec3 &o) const {
        return Vec3(data[0] * o[0], data[1] * o[1], data[2] * o[2]);
    }

    Vec3 operator/(const Vec3 &o) const {
        assert(o[0] != 0);
        assert(o[1] != 0);
        assert(o[2] != 0);
        return Vec3(data[0] / o[0], data[1] / o[1], data[2] / o[2]);
    }

    friend inline Vec3 operator+(const Vec3 &v, float s) {
        return Vec3(v[0] + s, v[1] + s, v[2] + s);
    }

    friend inline Vec3 operator+(float s, const Vec3 &v) {
        return v + s;
    }

    friend inline Vec3 operator-(const Vec3 &v, float s) {
        return Vec3(v[0] - s, v[1] - s, v[2] - s);
    }

    friend inline Vec3 operator-(float s, const Vec3 &v) {
        return v - s;
    }

    friend inline Vec3 operator*(const Vec3 &v, float s) {
        return Vec3(v[0] * s, v[1] * s, v[2] * s);
    }

    friend inline Vec3 operator*(float s, const Vec3 &v) {
        return v * s;
    }

    friend inline Vec3 operator/(const Vec3 &v, float s) {
        assert(s != 0);
        return Vec3(v[0] / s, v[1] / s, v[2] / s);
    }

    friend inline Vec3 operator/(float s, const Vec3 &v) {
        return v / s;
    }

    float squared_length() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }

    float length() const {
        return std::sqrt(squared_length());
    }

    Vec3 cwise_abs() const {
        return Vec3(std::abs(data[0]), std::abs(data[1]), std::abs(data[2]));
    }

private:
    float data[3];
};


//----------------------------------------------------------------------------------------------------------------------

struct AABB {
    AABB() : AABB(Vec3::Zero(), Vec3::Lowest()) {}

    AABB(const Vec3 &center, const Vec3 &halfsize) : center(center), halfsize(halfsize) {}

    void grow(const Vec3 &point) {
        Vec3 min = center - halfsize;
        Vec3 max = center + halfsize;

        min[0] = std::min(min[0], point[0]);
        min[1] = std::min(min[1], point[1]);
        min[2] = std::min(min[2], point[2]);

        max[0] = std::max(max[0], point[0]);
        max[1] = std::max(max[1], point[1]);
        max[2] = std::max(max[2], point[2]);

        halfsize = (max - min) / 2;
        center = min + halfsize;
    }

    bool is_inside_sphere(const Vec3 &point, float radius) const {
        return ((point - center).cwise_abs() + halfsize).squared_length() < radius * radius;
    }

    bool intersect(const Vec3 &point, float radius) const {
        Vec3 XYZ = (point - center).cwise_abs() - halfsize;

        if (XYZ[0] > radius || XYZ[1] > radius || XYZ[2] > radius) return false;
        int num_less_extent = (XYZ[0] < 0) + (XYZ[1] < 0) + (XYZ[2] < 0);
        if (num_less_extent > 1) return true;
        XYZ[0] = std::max<float>(XYZ[0], 0);
        XYZ[1] = std::max<float>(XYZ[1], 0);
        XYZ[2] = std::max<float>(XYZ[2], 0);

        return XYZ.squared_length() < radius * radius;
    }

    Vec3 center, halfsize;
};
//----------------------------------------------------------------------------------------------------------------------

template<class UserData>
struct OctreeBase {
    struct Node {
        Node(uint8_t config, uint64_t loc_code, size_t first_child_index, size_t parent_index, UserData data) :
                config(config),
                loc_code(loc_code),
                first_child_index(first_child_index),
                parent_index(parent_index),
                data(std::move(data)) {}

        uint8_t config;
        uint64_t loc_code;
        size_t first_child_index, parent_index;
        UserData data;
    };

    using container_t = std::vector<Node>;

    AABB root_aabb;
    container_t storage;

    struct ResultSet {
        using points_item_t  = std::pair<Vec3, float>;
        using points_container_t = std::vector<points_item_t>;
        using index_item_t  = std::pair<size_t, float>;
        using index_container_t = std::vector<index_item_t>;
        using iterator_t = index_container_t::iterator;
        using const_iterator_t = index_container_t::const_iterator;

        points_container_t points_dist_pair;
        index_container_t idx_dist_pair;
        float worst_dist;

        ResultSet() : worst_dist(std::numeric_limits<float>::max()), k(0) {}

        explicit ResultSet(size_t k) : worst_dist(0), k(k) {}

        bool add_point(const Vec3 &point, float dist) {
            if (!is_full() || dist < worst_dist) {
                worst_dist = std::max(worst_dist, dist);
                points_dist_pair.emplace_back(point, dist);
                return true;
            }
            return false;
        }

        bool add_point(size_t idx, float dist) {
            if (!is_full() || dist < worst_dist) {
                worst_dist = std::max(worst_dist, dist);
                idx_dist_pair.emplace_back(idx, dist);
                return true;
            }
            return false;
        }

        bool add_point(size_t idx) {
            idx_dist_pair.emplace_back(idx, 0);
            return true;
        }

        void sort() {
            if (!idx_dist_pair.empty()) {
                std::sort(idx_dist_pair.begin(), idx_dist_pair.end(),
                          [](const std::pair<size_t, float> &lhs, const std::pair<size_t, float> &rhs) {
                              return lhs.second < rhs.second;
                          });
                if (is_full()) {
                    idx_dist_pair.resize(k);
                    worst_dist = idx_dist_pair.back().second;
                }
            }
            if (!points_dist_pair.empty()) {
                std::sort(points_dist_pair.begin(), points_dist_pair.end(),
                          [](const std::pair<Vec3, float> &lhs, const std::pair<Vec3, float> &rhs) {
                              return lhs.second < rhs.second;
                          });
                if (is_full()) {
                    points_dist_pair.resize(k);
                    worst_dist = points_dist_pair.back().second;
                }
            }
        }

        bool is_full() const {
            return idx_dist_pair.size() >= k || points_dist_pair.size() >= k;
        }

        iterator_t begin() { return idx_dist_pair.begin(); }

        iterator_t end() { return idx_dist_pair.end(); }

        const_iterator_t begin() const { return idx_dist_pair.begin(); }

        const_iterator_t end() const { return idx_dist_pair.end(); }

    private:
        size_t k;
    };

    OctreeBase() = default;

    virtual ~OctreeBase() = default;

    virtual void clear() {
        root_aabb = AABB();
        storage.clear();
    }

    //function signature: std::function<bool(size_t index, Node &node, AABB &aabb)>
    void traverse_bfs(const std::function<bool(size_t, Node &, AABB &)> &function) {
        std::queue<AABB> queue;
        queue.emplace(root_aabb);

        for (size_t i = 0; i < storage.size() && !queue.empty(); ++i) {
            AABB &aabb = queue.front();
            if (function(i, storage[i], aabb)) {
                Vec3 child_extent = aabb.halfsize / 2;
                for (uint8_t octant = 0; octant < 8; ++octant) {
                    if (!child_exists(storage[i].config, octant)) continue;
                    queue.emplace(child_box(aabb, octant, child_extent));
                }
            }
            queue.pop();
        }
    }

    //function signature: std::function<bool(size_t index, const Node &node, const AABB &aabb)>
    void traverse_bfs(const std::function<bool(size_t, const Node &, const AABB &)> &function) const {
        std::queue<AABB> queue;
        queue.emplace(root_aabb);

        for (size_t i = 0; i < storage.size() && !queue.empty(); ++i) {
            const AABB &aabb = queue.front();
            if (function(i, storage[i], aabb)) {
                Vec3 child_extent = aabb.halfsize / 2;
                for (uint8_t octant = 0; octant < 8; ++octant) {
                    if (!child_exists(storage[i].config, octant)) continue;
                    queue.emplace(child_box(aabb, octant, child_extent));
                }
            }
            queue.pop();
        }
    }

    //function signature: std::function<bool(size_t index, Node &node, AABB &aabb)>
    void traverse_dfs(const std::function<bool(size_t, Node &, AABB &)> &function) {
        traverse_dfs(0, root_aabb, function);
    }

    //function signature: std::function<bool(size_t index, const Node &node, const AABB &aabb)>
    void traverse_dfs(const std::function<bool(size_t, const Node &, const AABB &)> &function) const {
        traverse_dfs(0, root_aabb, function);
    }

protected:

    //function signature: std::function<bool(size_t index, Node &node, AABB &aabb)>
    void traverse_dfs(size_t index, AABB &aabb, const std::function<bool(size_t, Node &, AABB &)> &function) {
        if (!function(index, storage[index], aabb)) return;

        Vec3 child_extent = aabb.halfsize / 2;
        size_t offset = 0;
        for (uint8_t i = 0; i < 8; ++i) {
            if (!child_exists(storage[index].config, i)) continue;
            traverse_dfs(storage[index].first_child_index + offset, child_box(aabb, i, child_extent), function); // WHYYYYYY??? Please! Never! Never ever use recursion in production code!
            ++offset;
        }
    }

    //function signature: std::function<bool(size_t index, const Node &node, const AABB &aabb)>
    void traverse_dfs(size_t index, const AABB &aabb,
                      const std::function<bool(size_t, const Node &, const AABB &)> &function) const {
        if (!function(index, storage[index], aabb)) return;

        Vec3 child_extent = aabb.halfsize / 2;
        size_t offset = 0;
        for (uint8_t i = 0; i < 8; ++i) {
            if (!child_exists(storage[index].config, i)) continue;
            traverse_dfs(storage[index].first_child_index + offset, child_box(aabb, i, child_extent), function);
            ++offset;
        }
    }

    static uint8_t morton_code(const Vec3 &point, const Vec3 &center) {
        uint8_t mortonCode = 0;
        if (point[0] > center[0]) mortonCode |= 1 << 0;
        if (point[1] > center[1]) mortonCode |= 1 << 1;
        if (point[2] > center[2]) mortonCode |= 1 << 2;
        return mortonCode;
    }

    static uint8_t morton_code(uint64_t loc_code, size_t depth) {
        return (((1 << 3) - 1) & (loc_code >> (depth * 3)));
    }

    static Vec3 center(const AABB &root_box, uint64_t loc_code) {
        int depth = tree_depth(loc_code);
        AABB box = root_box;
        for (int i = depth - 1; i >= 0; --i) {
            box = child_box(box, morton_code(loc_code, i), box.halfsize / 2);
        }
        return box.center;
    }

    static Vec3 child_center(const AABB &aabb, uint8_t i) {
        Vec3 child_center;
        child_center[0] = aabb.center[0] + (CHECK_BIT(i, 0) - 0.5) * aabb.halfsize[0];
        child_center[1] = aabb.center[1] + (CHECK_BIT(i, 1) - 0.5) * aabb.halfsize[1];
        child_center[2] = aabb.center[2] + (CHECK_BIT(i, 2) - 0.5) * aabb.halfsize[2];
        return child_center;
    }

    static AABB child_box(const AABB &aabb, uint8_t i, const Vec3 &child_extent) {
        return AABB(child_center(aabb, i), child_extent);
    }

    static Vec3 parent_center(const AABB &aabb, uint8_t i) {
        Vec3 parent_center;
        parent_center[0] = aabb.center[0] + CHECK_BIT(i, 0) * aabb.halfsize[0];
        parent_center[1] = aabb.center[1] + CHECK_BIT(i, 1) * aabb.halfsize[1];
        parent_center[2] = aabb.center[2] + CHECK_BIT(i, 2) * aabb.halfsize[2];
        return parent_center;
    }

    static AABB parent_box(const AABB &aabb, uint8_t i) {
        return AABB(parent_center(aabb, i), aabb.halfsize * 2);
    }

    static uint64_t child_loc_code(uint64_t loc_code, uint8_t i) {
        return loc_code << 3 | i;
    }

    static uint64_t parent_loc_code(uint64_t loc_code) {
        return loc_code >> 3;
    }

    static void set_child_exists(uint8_t &config, uint8_t i) {
        SET_BIT(config, i);
    }

    static bool child_exists(uint8_t config, uint8_t i) {
        return CHECK_BIT(config, i);
    }

    static size_t tree_depth(uint64_t loc_code) {
        return (63 - clz(loc_code)) / 3;
    }
};

//----------------------------------------------------------------------------------------------------------------------

struct OctantIndices {
    OctantIndices(size_t start, size_t end, size_t size) : start(start), end(end), size(size), idx(0) {}

    size_t start, end, size, idx;
};

struct Octree : public OctreeBase<OctantIndices> {
    const std::vector<Vec3> &points;
    size_t leaf_size;
    std::vector<size_t> indices;
    size_t depth;

    explicit Octree(const std::vector<Vec3> &points) :  points(points), leaf_size(1), depth(0) {};

    ~Octree() override = default;

    void clear() override;

    void build(size_t leaf_size = 1);

    ResultSet query_radius(const Vec3 &point, float radius, size_t depth = 21) const;

    ResultSet query_knn(const Vec3 &point, int k, size_t depth = 21) const;

protected:
    void init(size_t leaf_size);

    virtual void create_bfs();

    virtual void
    query_radius(size_t index, const AABB &aabb, const Vec3 &point, float radius, ResultSet &resultSet,
                 size_t depth) const;

    virtual void query_knn(size_t index, const AABB &aabb, const Vec3 &point, ResultSet &resultSet, size_t depth) const;

};

//}
#endif //BCG_OCTREE_H
