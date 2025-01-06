
```markdown
# Learned Spatial Indexes

## Overview
Learned Spatial Indexes is a C++ project focused on implementing and experimenting with spatial indexing techniques enhanced with search based on machine learning models.

## Prerequisites
- C++17 or later
- CMake 3.10 or later
- A compatible C++ compiler (e.g., GCC, Clang, MSVC)

## Usage
1. Clone the repository:
    ```sh
    git clone https://github.com/yourusername/learnedspatial.git
    cd learnedspatial
    ```

2. Create a build directory and navigate into it:
    ```sh
    mkdir build
    cd build
    ```

3. Configure the project using CMake:
    ```sh
    cmake ..
    ```

4. Build the project:
    ```sh
    make
    ```

## Usage
After building the project, you can run the executable with the following command:
```sh
./learnedspatial <points_file> <rectangles_file>
```

### Example
```sh
./learnedspatial data/points.csv data/rectangles.csv
```

## Configuration
We use several preprocessor directives to configure the benchmark type and other parameters. You can modify these in the source code or pass them as compiler flags.

### Benchmark Types
- `BENCHMARK_MATERIALIZE` - This option benchmarks the range queries. The results are materialized.
- `BENCHMARK_COUNT` - This option benchmarks the range queries but only reports the count of the result without materializeing the resulting points.
- `BENCHMARK_DISTANCE` - This options runs the distance queries. Instead of using a rectangles file as second arguement, this benchmark expects a file with each row in the format `lat,long,distance(in meters)`.
- `BENCHMARK_POINT` - This options runs the distance queries. Instead of using a rectangles file as second arguement, this benchmark expects a points file with each row in the format `lat,long`.
- `BENCHMARK_JOIN` - This options runs the point-in-polygon join queries. Instead of using a rectangles file as second arguement, this benchmark expects a polygon dataset as second argument in the format `field_a,field_b,polygon_wkt(the points int the polygon should be separated with a skipped comma i.e., '\,')`.

Examples:
```sh
g++ -DBENCHMARK_TYPE=BENCHMARK_MATERIALIZE -o learnedspatial_benchmark_materialize src/main.cpp
```

### Partition Size
By deafult the `PARTITION_SIZE` is set to `1000`. You can set the partition size using the `PARTITION_SIZE` directive. We plan to provide a yaml file with the best configurations that we found during experiments for all experiments in the paper soon. Please note that the best configurations are limited to the datasets that were considered in the paper. In case you need to use your own datasets, you will have to find those by experimenting yourself.

Example:
```sh
g++ -DPARTITION_SIZE=8000 -o learnedspatial_benchmark_materialize src/main.cpp
```

## Folder Structure
The following is the folder structure to navigate the code. Each partitioning technique is implemented in the parititoning_techniqeues folder. These techniques expect a partition or cell type. These could be one of those implemented in the the partition_cells. These partition cells implement the partition and a particular type of search within them. For instance, Spline.hpp implements the search based on Radix Spline. The partition cells also implement the logic for scan within partition as well as searching for a point within the partition.
```
learnedspatial/
├── CMakeLists.txt
├── README.md
├── include/
│   ├── ds/
│   │   ├── geography/
│   │   │   ├── DataTypes.hpp
│   │   │   └── ...
│   │   └── ...
│   ├── partition_cells/
│   │   ├── Spline.hpp
│   │   ├── BinarySearchY.hpp
│   │   ├── BinarySearchX.hpp
│   │   └── FullScan.hpp
│   ├── partitioning_techniques/
│   │   ├── AdaptiveGrid.hpp
│   │   ├── FixedGrid.hpp
│   │   ├── KdTreePartitioning.hpp
│   │   ├── QuadtreePartitioning.hpp
│   │   ├── SinglePartition.hpp
│   │   └── STRPartitioning.hpp
│   ├── queries/
│   │   ├── distance.hpp
│   │   └── join.hpp
│   ├── utils/
│   │   ├── IO.hpp
│   │   ├── SplineUtil.h
│   │   ├── Utils.hpp
│   │   └── ...
│   └── ...
├── src/
│   ├── main.cpp
│   └── ...
```

## License and Disclaimer
This project is licensed under the MIT License and provided as is.

## Contact
For any questions or suggestions, please open an issue.

## Cite

If you liked the work, please cite our papers in your work:

```
@inproceedings{learned_spatial_1,
  author       = {Varun Pandey and
                  Alexander van Renen and
                  Andreas Kipf and
                  Jialin Ding and
                  Ibrahim Sabek and
                  Alfons Kemper},
  title        = {The Case for Learned Spatial Indexes},
  booktitle    = {AIDB@VLDB 2020, 2nd International Workshop on Applied {AI} for Database
                  Systems and Applications, Held with {VLDB} 2020, Monday, August 31,
                  2020, Online Event / Tokyo, Japan},
  year         = {2020},
  url          = {https://drive.google.com/file/d/1Q\_kmSPhM86FeeZb8Kz196eNIn7uWcu8L/view?usp=sharing}
}
```

```
@article{learned_spatial_2,
  author       = {Varun Pandey and
                  Alexander van Renen and
                  Eleni Tzirita Zacharatou and
                  Andreas Kipf and
                  Ibrahim Sabek and
                  Jialin Ding and
                  Volker Markl and
                  Alfons Kemper},
  title        = {Enhancing In-Memory Spatial Indexing with Learned Search},
  journal      = {CoRR},
  volume       = {abs/2309.06354},
  year         = {2023},
  url          = {https://doi.org/10.48550/arXiv.2309.06354},
  doi          = {10.48550/ARXIV.2309.06354},
  eprinttype    = {arXiv},
  eprint       = {2309.06354}
}
```
