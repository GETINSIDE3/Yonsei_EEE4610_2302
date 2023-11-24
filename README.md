# Yonsei_EEE4610_2302
Graph similarity using CUDA

 
---

**Setting up the Repo**

    git clone https://github.com/GETINSIDE3/Yonsei_EEE4610_2302.git

    cd Yonsei_EEE4610_2302/

**Compiling your code**

    make clean; make


**Generating your custom graph for input**

    cd ./inputGen

    ./graphgen [num_vertices] [output_file_suffix]

    (or you can simply ./gen_dataset.sh)


**Testing the code**

    ./matmul [# of vertices to test] [input file name #1] [input file name #2]

    (input file names are optional)
