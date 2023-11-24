#include <iostream>
#include <fstream>
#include <limits>
#include <vector>
#include <cmath>
#include <cuda.h>

// Function to calculate the Frobenius norm of a matrix
double frobeniusNormCPU(const std::vector<std::vector<int>>& matrix) {
    double sum = 0.0;
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            int element = matrix[i][j];
            sum += element * element;
        }
    }
    return std::sqrt(sum);
}

// Function to calculate the similarity score using Cosine similarity
double calculateSimilarityScoreCPU(const std::vector<std::vector<int>>& graph1, const std::vector<std::vector<int>>& graph2) {
    // Check if the graphs have the same size
    if (graph1.size() != graph2.size() || graph1[0].size() != graph2[0].size()) {
        std::cerr << "Error: Graphs must have the same size." << std::endl;
        return -1;
    }

    // Calculate the Frobenius norms of the adjacency matrices
    double normGraph1 = frobeniusNormCPU(graph1);
    double normGraph2 = frobeniusNormCPU(graph2);

    // Calculate the similarity score
    double similarityScore = 0.0;

    for (size_t i = 0; i < graph1.size(); ++i) {
        for (size_t j = 0; j < graph1[0].size(); ++j) {
            similarityScore += (graph1[i][j] * graph2[i][j]) / (normGraph1 * normGraph2);
        }
    }

    return similarityScore;
}

// Function to print a graph (adjacency matrix)
void printMatrix(const std::vector<std::vector<int>>& graph, const std::string& label) {
    std::cout << "Graph " << label << ":" << std::endl;
    for (int i = 0; i < graph.size(); i++) {
        for (int j = 0; j < graph[i].size(); j++) {
            std::cout << graph[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}


std::vector<std::vector<int>> generateAdjacencyMatrix(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return {};
    }

    // Read the number of vertices from the first line
    unsigned long numVertices;
    file >> numVertices;

    // Initialize an adjacency matrix of the appropriate size with all zeros
    std::vector<std::vector<int>> adjacencyMatrix(numVertices, std::vector<int>(numVertices, 0));
    unsigned long edgepernode[numVertices];

    // Skip the next numVertices lines
    unsigned long totalEdges, numEdges;
    for (unsigned long i = 0; i < numVertices; ++i) {
        file >> totalEdges >> numEdges;
        edgepernode[i] = numEdges;
    }

    // Skip the next blank line
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Skip the next line containing a single number
    unsigned long skip;
    file >> skip;

    // Skip the next blank line
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Read the total number of edges from the next line
    file >> totalEdges;

    // Read the edges and fill in the adjacency matrix
    unsigned long dest;
    int weight;
    for (unsigned long i = 0; i < numVertices; ++i) {
        for(int j = 0; j < edgepernode[i]; ++j){
            // printf("Reading edge %d of vertex %lu\n", j, i);
            file >> dest >> weight;
            adjacencyMatrix[i][dest] = weight;
        }
    }

    file.close();
    return adjacencyMatrix;
}

std::vector<std::vector<int>> oppositeMatrix(const std::vector<std::vector<int>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();
    std::vector<std::vector<int>> opposite(cols, std::vector<int>(rows));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            opposite[i][j] = !matrix[i][j];
        }
    }

    return opposite;
}

std::vector<std::vector<int>> halfsameMatrix(const std::vector<std::vector<int>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();
    std::vector<std::vector<int>> halfsame(rows, std::vector<int>(cols));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols / 2; j++) {
            halfsame[i][j] = matrix[i][j];
        }
        for (int j = cols / 2; j < cols; j++) {
            halfsame[i][j] = !matrix[i][j];  // Generate random values
        }
    }
    
    return halfsame;
}

std::vector<std::vector<int>> generateRandomAdjacencyMatrix(int vertices, int weight_range) {
    std::vector<std::vector<int>> matrix(vertices, std::vector<int>(vertices));

    for (int i = 0; i < vertices; ++i) {
        for (int j = i + 1; j < vertices; ++j) {
            matrix[i][j] = rand() % weight_range;
            matrix[j][i] = matrix[i][j];  // For undirected graph
        }
    }

    return matrix;
}


// CUDA kernel to calculate the Frobenius norm of a matrix
__global__ void frobeniusNormKernel(int* matrix, int rows, int cols, double* result) {
    // TODO: Implement the kernel to calculate the Frobenius norm of a matrix
    return;
}

// Function to calculate the Frobenius norm of a matrix using CUDA
double frobeniusNormGPU(const std::vector<std::vector<int>>& matrix) {
    // TODO: Allocate memory on the GPU, copy the matrix to the GPU, launch the kernel,
    // copy the result back to the CPU, and free the GPU memory
    return -1.0;
}


// CUDA kernel to calculate the similarity score using
__global__ void calculateSimilarityScoreKernel(int* graph1, int* graph2, int rows, int cols, double normGraph1, double normGraph2, double* result) {
    // TODO: Implement the kernel to calculate the similarity score
    return;
}

// Function to calculate the similarity score using  with CUDA
double calculateSimilarityScoreGPU(const std::vector<std::vector<int>>& graph1, const std::vector<std::vector<int>>& graph2) {
    // TODO: Allocate memory on the GPU, copy the graphs to the GPU, launch the kernel,
    // copy the result back to the CPU, and free the GPU memory
    return -1.0;
}


int main(int argc, char* argv[]) {

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " [num_vertices_test] [input_file_1(optional)] [input_file_2(optional)]" << std::endl;
        return 1;
    }

    int TEST_VERTICES = std::stoi(argv[1]);

    // Test Case 1: Identical matrices, should result in score 1
    std::vector<std::vector<int>> identicalMatrix1 = generateRandomAdjacencyMatrix(TEST_VERTICES, 10);
    std::vector<std::vector<int>> identicalMatrix2 = identicalMatrix1;

    // // Test Case 2: Completely different matrices, should result in score 0
    std::vector<std::vector<int>> oppositeMatrix1 = generateRandomAdjacencyMatrix(TEST_VERTICES, 10);
    std::vector<std::vector<int>> oppositeMatrix2 = oppositeMatrix(oppositeMatrix1);

    // Test Case 3: Half of the elements are the same, should result in score 0.5
    std::vector<std::vector<int>> halfsameMatrix1 = generateRandomAdjacencyMatrix(TEST_VERTICES, 2);
    std::vector<std::vector<int>> halfsameMatrix2 = halfsameMatrix(halfsameMatrix1);

    double similarityScoreCPU = -1.0;
    double similarityScoreGPU = -1.0;

    // Test #1
    printf("Test 1: Identical matrices (expected score 1)\n");
    similarityScoreCPU = calculateSimilarityScoreCPU(identicalMatrix1, identicalMatrix2);
    similarityScoreGPU = calculateSimilarityScoreGPU(identicalMatrix1, identicalMatrix2);
    printf("CPU Similarity score: %.2f\n", similarityScoreCPU);
    printf("GPU Similarity score: %.2f\n", similarityScoreGPU);
    if (similarityScoreCPU > 0.99 && abs(similarityScoreCPU - similarityScoreGPU) < 0.01) {
        std::cout << "Test Case 1: Passed" << std::endl << std::endl;
    } else {
        printMatrix(identicalMatrix1, "1");
        printMatrix(identicalMatrix2, "2");
        std::cout << "Test Case 1: Failed" << std::endl << std::endl;
        exit(1);
    }

    // Test #2
    printf("Test 2: Opposite matrices (expected score 0)\n");
    similarityScoreCPU = calculateSimilarityScoreCPU(oppositeMatrix1, oppositeMatrix2);
    similarityScoreGPU = calculateSimilarityScoreGPU(oppositeMatrix1, oppositeMatrix2);
    printf("CPU Similarity score: %.2f\n", similarityScoreCPU);
    printf("GPU Similarity score: %.2f\n", similarityScoreGPU);
    if (similarityScoreCPU == 0.0 && abs(similarityScoreCPU - similarityScoreGPU) < 0.01) {
        std::cout << "Test Case 2: Passed" << std::endl << std::endl;
    } else {
        printMatrix(oppositeMatrix1, "1");
        printMatrix(oppositeMatrix2, "2");
        std::cout << "Test Case 2: Failed" << std::endl << std::endl; 
        exit(1);
    }

    // Test #3
    printf("Test 3: Half same matrices (expected score 0.5)\n");
    similarityScoreCPU = calculateSimilarityScoreCPU(halfsameMatrix1, halfsameMatrix2);
    similarityScoreGPU = calculateSimilarityScoreGPU(halfsameMatrix1, halfsameMatrix2);
    printf("CPU Similarity score: %.2f\n", similarityScoreCPU);
    printf("GPU Similarity score: %.2f\n", similarityScoreGPU);
    if (similarityScoreCPU > 0.49 && similarityScoreCPU < 0.51 && abs(similarityScoreCPU - similarityScoreGPU) < 0.01) {
        std::cout << "Test Case 3: Passed" << std::endl << std::endl;
    } else {
        printMatrix(halfsameMatrix1, "1");
        printMatrix(halfsameMatrix2, "2");
        std::cout << "Test Case 3: Failed" << std::endl << std::endl;
        exit(1);
    }

    if(argc > 2) {
        std::string inputFile1 = argv[2];
        std::string inputFile2 = argv[3];
        std::string inputDir   = "./inputGen/";

        printf("Test 4: Custom graphs\n");
    
        // Example adjacency matrices for two graphs
        std::vector<std::vector<int>> adjacencyMatrix1 = generateAdjacencyMatrix(inputDir+inputFile1);
        std::vector<std::vector<int>> adjacencyMatrix2 = generateAdjacencyMatrix(inputDir+inputFile2);

        // Calculate the similarity matrix
        similarityScoreCPU = calculateSimilarityScoreCPU(adjacencyMatrix1, adjacencyMatrix2);
        similarityScoreGPU = calculateSimilarityScoreGPU(adjacencyMatrix1, adjacencyMatrix2);

        if (abs(similarityScoreCPU - similarityScoreGPU) < 0.01) {
            std::cout << "Custom Graph Test: Passed" << std::endl << std::endl;
        } else {
            printMatrix(adjacencyMatrix1, "1");
            printMatrix(adjacencyMatrix2, "2");
            std::cout << "Custom Graph Test: Failed" << std::endl << std::endl;
            exit(1);
        }
    }

    return 0;
}
