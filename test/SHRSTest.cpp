#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

// std
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <utility>
#include <memory>
#include <chrono>
#include <random>
#include <catch2/catch.hpp>
#include <GridPoint.hpp>
#include <SHRS.hpp>

// PCL
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>

// Eigen
#include <Eigen/Dense>

using namespace fsd;

using std::string;
using std::cout;
using std::endl;
using std::vector;
using std::iota;
using std::numeric_limits;
using std::random_device;
using std::mt19937;
using std::uniform_real_distribution;
using std::uniform_int_distribution;
using std::unique_ptr;
using Eigen::Vector2d;

/*
 * @brief get_random_coordinates generates num number of unique and
 * uniform distributed random coordinates in the range of width and height,
 * in O(n) time with knuth-fisher-yates shuffling
 * can be optimized further for num > (width * heigh - 1) / 2
 * @param width
 * @param height
 * @param num num of random coords
 * @return vector with positions
 */
vector<GridPoint> generate_random_coordinates
(GridPoint area_min, GridPoint area_max, size_t num) {
    size_t width = std::abs(area_max.get_column() - area_min.get_column());
    size_t height = std::abs(area_max.get_row() - area_min.get_row());
    if(width == 0 || height == 0) {
        throw std::logic_error("width or height cannot be zero");
    }
    size_t max_val = width * height - 1;
    if(max_val <= num) {
        throw std::logic_error("number if random positions is bigger "
                               " or equal than max possible");
    }
    random_device rd;
    mt19937 eng(rd());
    vector<size_t> vals;
    vector<GridPoint> positions;
    vals.reserve(max_val + 1);
    positions.reserve(num);
    iota (begin(vals), end(vals), 0);
    while(num) {
        uniform_int_distribution<> u_distr(0, max_val);
        auto rnd_number = u_distr(eng);
        auto val_at_rnd = vals[rnd_number];
        vals[rnd_number] = vals[max_val];
        vals[max_val] = val_at_rnd;
        auto x = val_at_rnd % width;
        auto y = val_at_rnd / width;
        positions.emplace_back(GridPoint(x, y));
        max_val--;
        num--;
    }
    return positions;
}

/**
 *@brief adaptor for pointclouds
 */
template <typename point_t>
class PCLSHRSAdaptor {
public:
    using key_t = size_t;
    using data = pcl::PointCloud<point_t>;
    data const &pcl_ref;

    PCLSHRSAdaptor
    (data const &pcl_ref_) :
        pcl_ref(pcl_ref_) {
    }

    inline Vector2d get_point(size_t i_point) const {
        auto const &pnt = pcl_ref.at(i_point);
        return { pnt.x, pnt.y };
    }

    inline void for_all_keys(std::function<void(key_t)> func) const {
        for(size_t i = 0; i < pcl_ref.size(); i++) {
            func(i);
        }
    }

    inline bool find_point(size_t i_point) const {
        return i_point < pcl_ref.size();
    }

    inline size_t size() const {
        return pcl_ref.size();
    }
};

double bench_shrs(const double search_radius_pcl,
                  const Vector2d test_area_min_pcl,
                  const Vector2d test_area_max_pcl,
                  const size_t num_test_points_pcl,
                  const size_t points_to_search_pcl) {
    random_device rd;
    mt19937 eng(rd());

    using point_t = pcl::PointXYZ;
    pcl::PointCloud<point_t> test_pcl;

    SHRS<PCLSHRSAdaptor<point_t>> shrs_map(
                PCLSHRSAdaptor<point_t>(test_pcl),
                search_radius_pcl);

    bool gen_new_cloud = true;
    bool rw_cloud = false;
    pcl::PointCloud<pcl::PointXYZ>::Ptr read_cloud
            (new pcl::PointCloud<pcl::PointXYZ>);
    if (pcl::io::loadPCDFile<pcl::PointXYZ>
            ("/home/hsa/shrs_test_pcl.pcd", *read_cloud) == -1
            || gen_new_cloud) {
        // create new
        cout << "could not read cloud, creating new: ";
        uniform_real_distribution<> dis_x(test_area_min_pcl.x(),
                                          test_area_max_pcl.x());
        uniform_real_distribution<> dis_y(test_area_min_pcl.y(),
                                          test_area_max_pcl.y());
        // to save
        pcl::PointCloud<pcl::PointXYZ> save_cloud;
        save_cloud.width = 10;
        save_cloud.height = 100;
        save_cloud.is_dense = false;
        save_cloud.points.resize(num_test_points_pcl);
        for(size_t n = 0; n < num_test_points_pcl; n++) {
            pcl::PointXYZ tmp;
            auto rnd_x = dis_x(eng);
            auto rnd_y = dis_y(eng);
            tmp.x = rnd_x;
            tmp.y = rnd_y;
            tmp.z = 0.;
            save_cloud[n].x = tmp.x;
            save_cloud[n].y = tmp.y;
            save_cloud[n].z = 0;

            test_pcl.push_back(tmp);
        }
        if(rw_cloud) {
            pcl::io::savePCDFileASCII ("/home/hsa/shrs_test_pcl.pcd", save_cloud);
        }
    } else {
        cout << "read cloud from disk " << "\n""";
        test_pcl = *read_cloud;
    }

    auto begin_insert = std::chrono::steady_clock::now();

    shrs_map.insert_all();

    auto end_insert = std::chrono::steady_clock::now();
    auto insert_span = std::chrono::duration <double, std::milli>
            (end_insert - begin_insert).count();
    cout << "inserting " << num_test_points_pcl << " points took " << insert_span
         << " ms " << "\n";

    uniform_int_distribution<> u_distr_pcl(0, num_test_points_pcl - 1);
    double avg_neigbours = 0;
    double max_find = -std::numeric_limits<double>::infinity();
    double min_find =  std::numeric_limits<double>::infinity();
    double avg_find = 0;
    for(size_t iRndPoint = 0; iRndPoint < points_to_search_pcl; iRndPoint++) {
        bool random_idx = false;
        size_t rnd_index = 0;
        if(random_idx) {
            rnd_index = u_distr_pcl(eng);
        } else {
            rnd_index = iRndPoint;
        }
        auto rnd_pnt = test_pcl[rnd_index];
        vector<size_t> neighbours_shrs;
        vector<size_t> neighbours_reference;

        auto begin_find = std::chrono::steady_clock::now();

        shrs_map.find_neighbours(rnd_index,
                                 neighbours_shrs);

        auto end_find = std::chrono::steady_clock::now();
        double find_span = std::chrono::duration <double, std::milli>
                (end_find - begin_find).count();

        max_find = std::max(max_find, find_span);
        min_find = std::min(min_find, find_span);
        avg_find += find_span;
        avg_neigbours += (double)neighbours_shrs.size();

        sort(begin(neighbours_shrs), end(neighbours_shrs));

        for(size_t iOtherPoint = 0; iOtherPoint < num_test_points_pcl; iOtherPoint++) {
            if(iOtherPoint == rnd_index) {
                continue;
            }
            Vector2d  rnd_pnt_vec = { rnd_pnt.x, rnd_pnt.y };
            Vector2d  test_pcl_pnt_vec = { test_pcl[iOtherPoint].x,
                                           test_pcl[iOtherPoint].y };
            if((rnd_pnt_vec - test_pcl_pnt_vec).norm() <
                    (search_radius_pcl)) {
                neighbours_reference.push_back(iOtherPoint);
            }
        }
        sort(begin(neighbours_reference), end(neighbours_reference));
        if(neighbours_reference != neighbours_shrs) {
            cout << "found different neigbours ";
        }
        REQUIRE(neighbours_reference == neighbours_shrs);
    }

    avg_neigbours /= static_cast<double>(points_to_search_pcl);
    avg_find /= static_cast<double>(points_to_search_pcl);
    cout << "found " << avg_neigbours << " neighbours on average,"
                                         " one neighbour search took on average "
         << avg_find << " ms, max " << max_find << " ms and min "
         << min_find << " ms" << endl;


    // benchmark of adjacency matrix
    unordered_map<size_t, vector<size_t>> adj_matrix;
    auto begin_adj = std::chrono::steady_clock::now();

    unique_ptr<SHRS<PCLSHRSAdaptor<point_t>>> bench_map =
            unique_ptr<SHRS<PCLSHRSAdaptor<point_t>>>
                                                    (new SHRS<PCLSHRSAdaptor<point_t>>(PCLSHRSAdaptor<point_t>(test_pcl),
                                                                                       search_radius_pcl, &adj_matrix));

    auto end_adj = std::chrono::steady_clock::now();
    double adj_span = std::chrono::duration <double, std::milli>
            (end_adj - begin_adj).count();

    cout << "calculating adjacency list of " << num_test_points_pcl
         << " points with " << avg_neigbours
         << " neighbours on average took " << adj_span << " ms" << endl;
    return adj_span;
}

TEST_CASE("shrs insert and search neigbours with pcl adaptor") {

    vector<size_t> point_sizes_to_test = { 100, 1000, 10000, 100000 };
    vector<size_t> avg_neighbours_to_test = { 3, 7, 10, 20, 30, 50 };

    vector<vector<double>> times;
    auto const search_radius = 25.;
    for(auto size : point_sizes_to_test) {
        for(auto avg_neighbours : avg_neighbours_to_test) {
            auto area_side_len = std::sqrt(M_PI * std::pow(search_radius, 2.)
                                           * size
                                           / static_cast<double>(avg_neighbours));
            // to test also around the midpoint
            auto min_pnt = Vector2d(-(area_side_len / 2.),
                                    -(area_side_len / 2.));
            auto max_pnt = Vector2d(area_side_len / 2., area_side_len / 2.);

            cout << "min_pnt: x: " << min_pnt.x() << ", y: " << min_pnt.y()
                 << "max_pnt: x: " << max_pnt.x() << ", y: " << max_pnt.y();
            auto time_to_calc
                    = bench_shrs(search_radius, min_pnt, max_pnt, size,
                                 static_cast<size_t>(std::sqrt(size)));
            times.push_back({static_cast<double>(size),
                              static_cast<double>(avg_neighbours),
                              static_cast<double>(time_to_calc) });
        }
    }

    cout << "times: \n";

    for(auto const &time_struct : times) {
        cout << "size: " << time_struct[0] << ", "
             << "avg_neighbours: " << time_struct[1] << ", "
             << "time: " << time_struct[2] << " ms\n";
    }
}


