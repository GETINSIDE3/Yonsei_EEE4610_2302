# Yonsei_EEE4610_2302
Graph similarity using CUDA

 
---

 
**Compiling your code**

    make clean; make


**Generating your custom graph for input**

    cd ./inputGen

    ./graphgen [num_vertices] [output_suffix]

    (or you can simply ./gen_dataset.sh)


**Testing the code**

    ./matmul [number of vertices to test] [input file name #1 (optional)] [input file name #2 (optional)]
