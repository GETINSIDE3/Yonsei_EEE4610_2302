FLAGS :=

matmul: matmul.cu
	nvcc $(FLAGS) -o $@ $<

clean: 
	rm -f matmul
