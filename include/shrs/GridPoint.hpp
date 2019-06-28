#ifndef PROJECT_GRIDPOINT_HPP
#define PROJECT_GRIDPOINT_HPP

#include <stdint.h>
#include <string>
#include <limits>

#include <Eigen/Dense>

namespace fsd {

using Eigen::Vector2d;
using std::string;
/**
 * @brief The GridPoint class representing a point in a grid by row, column, for
 * fast lookup
 */
class GridPoint {
public:
    /**
     * @brief the type of point this gridpoint is constructed from,
     * must have a (x, y)/coordinate constructor and x(), y()
     * / cordinate getter methods. and probably more ..
     *
    **/
    typedef Vector2d from_point_t;

    typedef int32_t coordinates_t;
    typedef int64_t key_t;
    /**
     * must be doulbe for packing, unpacking
     */
    static_assert(sizeof(coordinates_t) * 2 ==
                  sizeof(key_t), "key_t must be double "
                                 "the size of coordinates_t for"
                                 " packing and stuff");
    /**
     * @brief getter for row
     * @return row
     */
    coordinates_t get_row() const {
        return row;
    }

    /**
     * @brief getter for column
     * @return column
     */
    coordinates_t  get_column() const {
        return column;
    }

    /**
     * @brief hit_cnt how many points mapped to this cell
     */
    uint32_t hit_cnt = 1;

    /**
     * @brief The coordinates struct used for other grid things
     */
    struct coordinates {
        coordinates_t row, column;
        coordinates(coordinates_t row_,
                    coordinates_t column_) :
            row(row_), column(column_) {}
        coordinates() = default;
    };

    static_assert (sizeof (coordinates) == sizeof (key_t), "key size must equal "
                                                           "coordinates");

    /**
     * @brief to_key creates collision free key by packing row and column
     * in one integer
     * @return key
     */
    key_t to_key() const {
        coordinates coords;
        // copy instead of brace initializer to avoid some bugs
        // if the struct is redefined
        coords.row = row;
        coords.column = column;
        return pack_key(coords);
    }

    /**
     * @brief to_coordinates point to coordinates helper
     * @param point
     * @param cell_size
     * @param scale to use more entropy on the coordinates
     * and therefore have less collisions when saving in a hashtable,
     * the coordinates can be scaled the use the whole range
     * of coordinates_t type
     * @return
     */
    static inline coordinates to_coordinates(
            const from_point_t &point, double cell_size, double scale = 1.) {
        double x = point.x() * scale;
        if(x < 0.) {
            x -= cell_size;
        } else {
            // also zero
            x += cell_size;
        }
        double y = point.y() * scale;
        if(y < 0.) {
            y -= cell_size;
        } else {
            // also zero
            y += cell_size;
        }
        // cast is used to truncate towards zero
        coordinates coord;
        coord.row = static_cast<coordinates_t>(x / cell_size);
        coord.column = static_cast<coordinates_t>(y / cell_size);
        return coord;
    }

//#define MIXING8
    /**
     * @brief pack_key helper to pack key, this zigzac
     * packing is used
     * to perform better with std::hash identity hash
     * funtion, the lower bytes are packed,
     * as there is the most entropy
     * @param coords
     * @return
     */
    static inline key_t pack_key(coordinates coords) {
        key_t packed;
#ifdef MIXING8
        uint8_t *packed_ptr = ((uint8_t*)&packed);
        uint8_t *row_ptr = ((uint8_t*)&(coords.row));
        uint8_t *column_ptr = ((uint8_t*)&(coords.column));
        packed_ptr[0] = row_ptr[0];
        packed_ptr[1] = column_ptr[0];
        packed_ptr[2] = row_ptr[1];
        packed_ptr[3] = column_ptr[1];
        packed_ptr[4] = row_ptr[2];
        packed_ptr[5] = column_ptr[2];
        packed_ptr[6] = row_ptr[3];
        packed_ptr[7] = column_ptr[3];
#else
        ((uint32_t*)&packed)[0] = *((uint32_t*)&(coords.row));
        ((uint32_t*)&packed)[1] = *((uint32_t*)&(coords.column));
#endif
        return packed;
    }

    static inline coordinates unpack_key(key_t key) {
        coordinates coords;
#ifdef MIXING8
        uint8_t *packed_ptr = ((uint8_t*)&key);
        uint8_t *row_ptr = ((uint8_t*)&(coords.row));
        uint8_t *column_ptr = ((uint8_t*)&(coords.column));
        row_ptr[0] = packed_ptr[0];
        column_ptr[0] = packed_ptr[1];
        row_ptr[1] = packed_ptr[2];
        column_ptr[1] = packed_ptr[3];
        row_ptr[2] = packed_ptr[4];
        column_ptr[2] = packed_ptr[5];
        row_ptr[3] = packed_ptr[6];
        column_ptr[3] = packed_ptr[7];
#else
        *((uint32_t*)&(coords.row)) = ((uint32_t*)&key)[0];
        *((uint32_t*)&(coords.column)) = ((uint32_t*)&key)[1];
#endif
        return  coords;
    }

    /**
     * @brief to_key fast conversion to key. Allows checking in which
     * cell this point would map without first constructing the point
     * @param point
     * @param cell_size
     * @return key
     */
    static inline key_t to_key(const from_point_t &point,
                               double cell_size, double scale = 1.) {
        auto coords = to_coordinates(point, cell_size, scale);
        return pack_key(coords);
    }

    /**
     * @brief coordinates_would_overflow
     * @param point
     * @param cell_size
     * @return true if coordinates would overflow
     */
    static bool coordinates_would_overflow(
            const from_point_t &point, double cell_size) {
        static_assert (sizeof (coordinates_t) <= 4,
                       "coordinates_would_overflow check unfortunately"
                       "isn't implemented for coordinates_t type bigger than "
                       "4 byte");
        double x = point.x();
        if(x > 0.) {
            x += cell_size;
        } else if(x < 0.){
            x -= cell_size;
        }
        double y = point.y();
        if(y > 0.) {
            y += cell_size;
        } else if(y < 0.) {
            y -= cell_size;
        }

        double row = x / cell_size;
        double column = y / cell_size;
        bool row_overflow = row > std::numeric_limits<coordinates_t>::max()
                || row < std::numeric_limits<coordinates_t>::min();
        bool column_overflow = column > std::numeric_limits<coordinates_t>::max()
                || column < std::numeric_limits<coordinates_t>::min();
        return row_overflow || column_overflow;
    }

    /**
     * @brief to_point to_point creates a point with a datatype which is suitable
     * for doing geometric calculations
     * @param cell_size
     * @param scale scale by which the coordinates where multiplied,
     * so the inverse is taken to get the real coordinates
     * @return
     */
    from_point_t to_point(double cell_size, double scale = 1.) const {
        double tmp_x, tmp_y;
        auto scale_inv = 1. / scale;
        if (row < 0) {
            tmp_x = ((row * cell_size) + .5 * cell_size) * scale_inv;
        } else if(row > 0){
            tmp_x = (((row - 1) * cell_size) + .5 * cell_size) * scale_inv;
        } else {
            tmp_x = 0;
        }
        if (column < 0) {
            tmp_y = ((column * cell_size) + .5 * cell_size) * scale_inv;
        } else if(column > 0){
            tmp_y = (((column - 1) * cell_size) + .5 * cell_size) * scale_inv;
        } else {
            tmp_y = 0;
        }
        return from_point_t(tmp_x, tmp_y);

    }

    /**
     * @brief to_point static method to build a point
     * from key directly
     * @param key
     * @param cell_size
     * @param scale scale by which the coordinates where multiplied,
     * so the inverse is taken to get the real coordinates
     * @return
     */
    static from_point_t to_point(key_t key, double cell_size, double scale = 1.) {
        coordinates coords = unpack_key(key);
        coordinates_t row = coords.row;
        coordinates_t column = coords.column;
        auto scale_inv = 1. / scale;
        double tmp_x, tmp_y;
        if (row < 0) {
            tmp_x = ((row * cell_size) + .5 * cell_size) * scale_inv;
        } else if(row > 0){
            tmp_x =  (((row - 1) * cell_size) + .5 * cell_size) * scale_inv;
        } else {
            tmp_x = 0;
        }

        if (column < 0) {
            tmp_y = ((column * cell_size) + .5 * cell_size) * scale_inv;
        } else if(column > 0){
            tmp_y = (((column - 1) * cell_size) + .5 * cell_size) * scale_inv;
        } else {
            tmp_y = 0;
        }

        return from_point_t(tmp_x, tmp_y);

    }
    /**
     * @brief has_same_grid_position
     * @param other_point
     * @return true if it has same row and column
     */
    bool has_same_grid_position(const GridPoint &other_point) const {
        return row == other_point.row &&
                column == other_point.column;
    }

    /**
     * @brief GridPoint construct gridpoint from vector and cell size
     * @param point
     * @param cell_size
     */
    GridPoint(from_point_t point, double cell_size, double scale = 1.) {
        auto coord = to_coordinates(point, cell_size, scale);
        row = coord.row;
        column = coord.column;
    }

    GridPoint(coordinates_t row, coordinates_t column) :
        row(row), column(column) {

    }

    GridPoint() = default;

    /**
     * @brief operator + adds coordinates of another point by adding row and column
     * @param other_point
     * @return
     */
    GridPoint operator+(const GridPoint &other_point) const {
        GridPoint tmp(other_point);
        tmp.row += row;
        tmp.column += column;
        return tmp;
    }

    /**
     * @brief to string
     * @return
     */
    virtual string to_string() const {
        return "point: row: " + std::to_string(row) +
                ", column: " + std::to_string(column) +
                ", hit_cnt: " + std::to_string(hit_cnt);
    }

    virtual ~GridPoint() {}

protected:
    coordinates_t row = 0;
    coordinates_t column = 0;
};

} // end namespace fsd

#endif //PROJECT_GRIDPOINT_HPP
