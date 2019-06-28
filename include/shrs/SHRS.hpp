#ifndef SHRS_HPP
#define SHRS_HPP

// std
#include <string>
#include <vector>
#include <utility>
#include <unordered_map>
#include <unordered_set>

// Eigen
#include <Eigen/Dense>

// PCL
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>

#include <GridPoint.hpp>

namespace fsd {

using Eigen::Vector2d;
using std::vector;
using std::unordered_map;
using std::unordered_set;
using std::begin;
using std::end;
using std::endl;

/**
 * @brief The SHRS class spatial hashing radius search,
 * some data structure
 * for fast radius search in sparse pointcloud.
 * currently implemented only for 2D
 * and only for constant search_radius as
 * the cell_size is search_radius * 2.,
 * and no search circle rasterizing is implemented
 * yet but only getting the 1 - 4 cells which
 * this search circle can intersect
 * @tparam non-virtual dataset interface
 */
template<typename DatasetAdaptor>
class SHRS { 
public:
    using key_t = typename DatasetAdaptor::key_t;
    using coords = GridPoint::coordinates;
    using hash_func_t = std::hash<key_t>;
    using cell_t = unordered_set<key_t, hash_func_t>;
    using search_map_t = unordered_map<key_t, cell_t, hash_func_t>;

    /**
      *@brief constructor; if adj_matrix is passed, its calculated
      */
    SHRS(DatasetAdaptor dataset_adaptor, double search_radius_,
         unordered_map<key_t, vector<key_t>> *adj_matrix = nullptr) :
        dataset_adaptor(dataset_adaptor), cell_size(search_radius_ * 2.),
        search_radius(search_radius_),
        search_radius_2(pow(search_radius_, 2.)) {
        // curently only search radius of cell_size
        // implemented
        cell_keys_to_search.reserve(4);
        search_map.max_load_factor(0.5);
        if(adj_matrix) {
            calculate_adj_matrix(*adj_matrix);
        }
    }

    SHRS &set_initial_cell_size(size_t value) {
        initial_cell_size = value;
        return *this;
    }

    /**
     * @brief find_neighbours find
     * all neigbours within set radius
     * @param key key of point
     * @return number of neighbours
     */
    size_t find_neighbours(key_t point_key,
                           vector<key_t> &neighbours) {
        neighbours.reserve(initial_neighbours_size);
        return on_find_neighbours(point_key, [&](key_t key) {
            neighbours.push_back(key);
        });
    }

    /**
     * @brief on_find_neighbours found neigbours
     * with callback to eventually avoid copy
     * @param point_key
     * @param func
     * @return
     */
    size_t on_find_neighbours(key_t point_key,
                              std::function<void(key_t)> func) {
        auto vec_pnt = dataset_adaptor.get_point(point_key);
        auto cell_coords =
                GridPoint::to_coordinates(vec_pnt, cell_size);
        auto cell_key = GridPoint::pack_key(cell_coords);
        if(search_map.find(cell_key) == end(search_map)) {
            // this point doesn't exist in search_map
            return 0;
        }
        size_t neighbours_found = 0;
        get_cell_keys_in_range(vec_pnt, cell_keys_to_search);
        // now search all cells ( maximum 4 cells are searched)
        for(auto cell_key_to_search : cell_keys_to_search) {
            if(search_map.find(cell_key_to_search) == end(search_map)) {
                continue;
            }
            auto &point_keys = search_map.at(cell_key_to_search);
            for(auto other_point_key : point_keys) {
                if(point_key == other_point_key) {
                    continue;
                }
                auto other_point =
                        dataset_adaptor.get_point(other_point_key);
                double distance2 = (vec_pnt - other_point).squaredNorm();
                if(distance2 < search_radius_2) {
                    func(other_point_key);
                    neighbours_found++;
                }
            }
        }
        return neighbours_found;
    }


    /**
     * @brief find_neighbours_vec find neigbours to a point
     * not in dataset
     * @param vec_pnt
     * @param neighbours
     * @return num of neigbours
     */
    size_t find_neighbours_vec(Vector2d const &vec_pnt,
                               vector<key_t> &neighbours) {
        neighbours.reserve(initial_neighbours_size);
        return on_find_neighbours_vec(vec_pnt, [&](key_t key) {
            neighbours.push_back(key);
        });
    }


    /**
     * @brief find_neighbours_vec find neigbours to a point
     * not in dataset
     * @param vec_pnt
     * @param neighbours
     * @return num of neigbours
     */
    size_t on_find_neighbours_vec(Vector2d const &vec_pnt,
                                  std::function<void(key_t)> func) {
        auto cell_coords =
                GridPoint::to_coordinates(vec_pnt, cell_size);
        size_t neighbours_found = 0;
        get_cell_keys_in_range(vec_pnt, cell_keys_to_search);
        // now search all cells ( maximum 4 cells are searched)
        for(auto cell_key_to_search : cell_keys_to_search) {
            if(search_map.find(cell_key_to_search) == end(search_map)) {
                continue;
            }
            auto &point_keys = search_map.at(cell_key_to_search);
            for(auto other_point_key : point_keys) {
                auto other_point =
                        dataset_adaptor.get_point(other_point_key);
                double distance2 = (vec_pnt - other_point).squaredNorm();
                if(distance2 < search_radius_2) {
                    func(other_point_key);
                    neighbours_found++;
                }
            }
        }
        return neighbours_found;
    }

    /**
     * @brief find_nn_vec finds like find_neighbours_vec, but
     * nearest neigbour, which is in search radius
     * @param vec_pnt
     * @param neighbours
     * @return
     */
    size_t find_nn_vec(Vector2d const &vec_pnt,
                       vector<key_t> &neighbours) {
        auto cell_coords =
                GridPoint::to_coordinates(vec_pnt, cell_size);
        size_t neighbours_found = 0;
        neighbours.resize(1);
        double min_distance = std::numeric_limits<double>::infinity();
        get_cell_keys_in_range(vec_pnt, cell_keys_to_search);
        // now search all cells ( maximum 4 cells are searched)
        for(auto cell_key_to_search : cell_keys_to_search) {
            if(search_map.find(cell_key_to_search) == end(search_map)) {
                continue;
            }
            auto &point_keys = search_map.at(cell_key_to_search);
            for(auto other_point_key : point_keys) {
                if(points_index.find(other_point_key) == points_index.end()){
                    continue;
                }
                auto other_point =
                        dataset_adaptor.get_point(other_point_key);
                double distance2 = (vec_pnt - other_point).squaredNorm();
                if(distance2 < search_radius_2) {
                    if(distance2 < min_distance) {
                        min_distance = distance2;
                        neighbours[0] = other_point_key;
                        neighbours_found = 1;
                    }
                }
            }
        }
        if(!neighbours_found) {
            neighbours.clear();
        }
        return neighbours_found;
    }


    /**
      *
      * @brief inserts all points and calaculates adj_matrix
      */
    void calculate_adj_matrix(unordered_map<key_t, vector<key_t>> &adj_matrix) {
        insert_all();
        auto num_points = points_index.size();
        adj_matrix.max_load_factor(0.5);
        adj_matrix.reserve(num_points);
        dataset_adaptor.for_all_keys([&](key_t key) {
            find_neighbours(key, adj_matrix[key]);
        });
    }

    /**
     * @brief insert_point inserts the point in the seach_map,
     * the point is referenced by key, this means in must be first
     * have inserted in the underlying datastructure provided by
     * DatasetAdaptor
     * @param point_key
     * @return true if inserted
     */
    void insert_point(key_t point_key) {
        auto vec_pnt = dataset_adaptor.get_point(point_key);
        auto cell_coords =
                GridPoint::to_coordinates(vec_pnt, cell_size);
        auto cell_key = GridPoint::pack_key(cell_coords);
        if(search_map.find(cell_key) == end(search_map)) {
            auto set_to_insert = cell_t(initial_cell_size);
            set_to_insert.max_load_factor(.5f);
            search_map.emplace(cell_key, set_to_insert);
        }
        search_map.at(cell_key).insert(point_key);
        points_index.insert(point_key);
    }

    /**
     * @brief insert_all inserts all points currently in the
     * dataset
     */
    void insert_all() {
        points_index.reserve(dataset_adaptor.size());
        dataset_adaptor.for_all_keys([&](key_t key) {
            insert_point(key);
        });
    }

    /**
     * @brief delete_point deleted point must be available
     * while deleting, it can be deleted only
     * aftewards from the underlying datastructure provided by
     * DatasetAdaptor
     * @param point_key
     * @return true if deleted
     */
    bool delete_point(key_t point_key) {
        if(!dataset_adaptor.find_point(point_key)) {
            return false;
        }
        auto vec_pnt = dataset_adaptor.get_point(point_key);
        auto cell_key = GridPoint::to_key(vec_pnt, cell_size);
        auto &points = search_map.at(cell_key);
        points.erase(point_key);
        points_index.erase(point_key);
        if(points.empty()) {
            search_map.erase(cell_key);
        }
        return true;
    }

    /**
     * @brief get_cell_keys_in_range gets all cells on which
     * a circle with diameter cell_size and the center vec_pnt
     * lays
     * @param vec_pnt the point which is the center
     * @param OUT keys of these cells
     *
     */
    void get_cell_keys_in_range(Vector2d const &vec_pnt,
                                vector<key_t> &keys) {
        keys.clear();
        auto cell_coords =
                GridPoint::to_coordinates(vec_pnt, cell_size);
        keys.push_back(GridPoint::pack_key(cell_coords));
        // offset of this point from center
        double const half_cell_size = cell_size / 2.;
        double row_d = vec_pnt.x();
        double column_d = vec_pnt.y();
        double row_offset = fmod(row_d, cell_size);
        double column_offset = fmod(column_d, cell_size);
        auto fourth_point = cell_coords;
        if(row_offset < 0) {
            row_offset += cell_size;
        }
        if(column_offset < 0) {
            column_offset += cell_size;
        }
        // to jump from -1 to 1 evtl, as
        // 0 is not used as row and column
        if(row_offset < half_cell_size) {
            if(cell_coords.row == 1) {
                fourth_point.row -= 2;
                keys.push_back(GridPoint::pack_key(coords(cell_coords.row - 2,
                                                          cell_coords.column)));
            } else {
                fourth_point.row -= 1;
                keys.push_back(GridPoint::pack_key(coords(cell_coords.row - 1,
                                                          cell_coords.column)));
            }
        } else if(row_offset > half_cell_size) {
            if(cell_coords.row == -1) {
                fourth_point.row += 2;
                keys.push_back(GridPoint::pack_key(coords(cell_coords.row + 2,
                                                          cell_coords.column)));
            } else {
                fourth_point.row += 1;
                keys.push_back(GridPoint::pack_key(coords(cell_coords.row + 1,
                                                          cell_coords.column)));
            }
        }
        if(column_offset < half_cell_size) {
            if(cell_coords.column == 1) {
                fourth_point.column -= 2;
                keys.push_back(GridPoint::pack_key(coords(cell_coords.row,
                                                          cell_coords.column - 2)));
            } else {
                fourth_point.column -= 1;
                keys.push_back(GridPoint::pack_key(coords(cell_coords.row,
                                                          cell_coords.column - 1)));
            }
        } else if(column_offset > half_cell_size) {
            if(cell_coords.column == -1) {
                fourth_point.column += 2;
                keys.push_back(GridPoint::pack_key(coords(cell_coords.row,
                                                          cell_coords.column + 2)));
            } else {
                fourth_point.column += 1;
                keys.push_back(GridPoint::pack_key(coords(cell_coords.row,
                                                          cell_coords.column + 1)));
            }
        }

        if(keys.size() == 3) {
            keys.push_back(GridPoint::pack_key(fourth_point));
        }
    }

    /**
     * @brief get_points_index
     * @return
     */
    cell_t const &get_points_index() const {
        return points_index;
    }

    void set_initial_neighbours_size(size_t value) {
        initial_neighbours_size = value;
    }

    /**
     * @brief clear clears the whole shrs map, but not the
     * undelying data structure
     */
    void clear() {
        points_index.clear();
        search_map.clear();
    }

    void delete_from_index(key_t idx) {
        points_index.erase(idx);
    }

protected:
    /**
     * @brief initial_cell_size initial size to reserve
     * on a new cell when constructing
     */
    size_t initial_cell_size = 7;

    /**
     * @brief initial_neighbours_size initial size to
     * allocate for neigbours
     */
    size_t initial_neighbours_size = 5;

    /**
     * @brief search_map
     */
    search_map_t search_map;

    /**
     * @brief dataset_adaptor
     */
    DatasetAdaptor dataset_adaptor;

    /**
     * @brief cell_size
     */
    double cell_size;

    /**
     * @brief search_radius
     */
    double search_radius;

    /**
     * @brief search_radius_2 squared_search_radius
     */
    double search_radius_2;

    /**
     * @brief cells_to_search cells to search on
     * radius query, memebr to avoid realocation
     */
    vector<key_t> cell_keys_to_search;

    /**
     * @brief points_index which points are saved
     * in the shrs grid
     */
    cell_t points_index;

};

} // end namespace fsd

#endif
