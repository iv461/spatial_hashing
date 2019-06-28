## Spatial Hashing Radius Search

A data structure for finding 2D points within a radius fast,
best suited for sparse pointcloud where the points are uniformly distributed

## Dependencies

This depends on Eigen and PCL library, for the example
dataset adaptor and catch2 which is included

## Building

This is a single-header library, a unit test executable is provided,
build it for example on Linux with:

```
mkdir build
cd build
cmake ..
make
```

And run it:

```
./shrs_test
```

## Usage

Just include the ```SHRS.hpp``` header in your project

### Dataset adaptor

Similar to nanoflann, you have to provide a dataset adaptor.

A example adaptor for PCL pointclouds:

```
template <typename point_t>
class PCLSHRSAdaptor {
public:
    // key_t needed for SHRS , must be defined
    using key_t = size_t;
    using data = pcl::PointCloud<point_t>;
    data const &pcl_ref;
    // the already existing pointcloud is passed
    PCLSHRSAdaptor
    (data const &pcl_ref_) :
        pcl_ref(pcl_ref_) {
    }
    // SHRS works with 2D points, this must return a Vector2d
    inline Vector2d get_point(size_t i_point) const {
        auto const &pnt = pcl_ref.at(i_point);
        return { pnt.x, pnt.y };
    }
    // needed for insert_all() and similar operations
    inline void for_all_keys(std::function<void(key_t)> func) const {
        for(auto i = 0; i < pcl_ref.size(); i++) {
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
```
### Constructing search data structure


```
using point_t = pcl::PointXYZ;
using shrs_t = SHRS<PCLSHRSAdaptor<point_t>>;
pcl::PointCloud<point_t> test_pcl;
auto const search_radius = .5;
shrs_t shrs_map(PCLSHRSAdaptor<point_t>(test_pcl),
                search_radius);
```
## Inserting and deleting points

Now you can fill the pointcloud with data and then insert all the points
in the search data structure:

``` 
// fill here the test_pcl ..

// now insert in search data structure
shrs_map.insert_all();

```

Points can also be inserted and deleted one by one in constant time complexity:

```
size_t some_index = 11;
shrs_map.insert_point(some_index);
shrs_map.delete_point(some_index);
```

## Examples

### Finding all points in radius for a point in dataset

```
// must exist in pointcloud
size_t search_point_idx = 7; 
vector<typename shrs_t::key_t> neighbours;
if(shrs_map.find_neighbours(search_point_idx, neighbours) {
    // found neighbour
    for(auto const n_keys : neighbours) {
       auto point_in_radius = test_pcl[n_keys];
       // use the point_in_radius ..
    }
}
```

### Finding nearest neighbour in radius of a point not in dataset

```
Vector2d search_point = { 70., 80. };
vector<typename shrs_t::key_t> neighbours;
if(shrs_map.find_nn_vec(search_point, neighbours)) {
    // found neighbour
    auto nn_index = neighbours.front();
    auto nn_point = test_pcl[nn_index];
    // use the nn_point ..
}
```


## Performance

Example output of the unit test measuring the performance:

```
could not read cloud, creating new: inserting 10000 points took 3.00308 ms 
found 7.02 neighbours on average, one neighbour search took on average 0.00214091 ms, max 0.015111 ms and min 0.000679 ms
calculating adjacency list of 10000 points with 7.02 neighbours on average took 13.4824 ms
```

Adjacency list calculation times 
for different sizes and average neighbours

```
size: 100, avg_neighbours: 3, time: 0.077742 ms
size: 100, avg_neighbours: 7, time: 0.089566 ms
size: 100, avg_neighbours: 10, time: 0.103969 ms
size: 100, avg_neighbours: 20, time: 0.149846 ms
size: 100, avg_neighbours: 30, time: 0.17453 ms
size: 100, avg_neighbours: 50, time: 0.200177 ms
size: 1000, avg_neighbours: 3, time: 0.934259 ms
size: 1000, avg_neighbours: 7, time: 1.09438 ms
size: 1000, avg_neighbours: 10, time: 1.41327 ms
size: 1000, avg_neighbours: 20, time: 2.16976 ms
size: 1000, avg_neighbours: 30, time: 2.76571 ms
size: 1000, avg_neighbours: 50, time: 4.04864 ms
size: 10000, avg_neighbours: 3, time: 9.99308 ms
size: 10000, avg_neighbours: 7, time: 15.0083 ms
size: 10000, avg_neighbours: 10, time: 16.0602 ms
size: 10000, avg_neighbours: 20, time: 24.8878 ms
size: 10000, avg_neighbours: 30, time: 34.01 ms
size: 10000, avg_neighbours: 50, time: 50.6205 ms
size: 100000, avg_neighbours: 3, time: 221.516 ms
size: 100000, avg_neighbours: 7, time: 331.998 ms
size: 100000, avg_neighbours: 10, time: 423.088 ms
size: 100000, avg_neighbours: 20, time: 697.058 ms
size: 100000, avg_neighbours: 30, time: 1007.68 ms
size: 100000, avg_neighbours: 50, time: 1527.43 ms


```

## License

MIT


