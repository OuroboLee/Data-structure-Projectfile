#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define MAX_QUEUE_SIZE 500
#define DEPTH 2

void printArray(int* target, int size) {
    for (int i = 0; i < size; i++) 
        printf("%d ", target[i]);
    printf("\n");
}

void printArray(double* target, int size) {
    for (int i = 0; i < size; i++) 
        printf("%f ", target[i]);
    printf("\n");
}


typedef struct slice {
    // int depth;      // = (the size of rateOfChangeArray) - 1
    int sliceSize; // the size of slice.
    int start, end; // slicing start idx / end idx
    int* baseDemandArray; // 기초가 되는 수요 배열 주소
    int** rateOfChangeArray;
};

void constructRateOfChangeArray(slice* target) {
    //printf("==============================\n");
    //printf("start : %d, end : %d\n", target -> start, target -> end);
    for (int i = 0; i < target -> sliceSize; i++) {
        (target -> rateOfChangeArray)[0][i] = target -> baseDemandArray[(target -> start) + i];
        //printf("base[(target -> start) + %d] : %d\n", i, base[(target -> start) + i]);
    }
    
    for (int i = 1; i < DEPTH + 1; i++) {
        for (int j = 0; j < (target -> sliceSize) - i; j++) {
            (target -> rateOfChangeArray)[i][j] = (target -> rateOfChangeArray)[i-1][j+1] - (target -> rateOfChangeArray)[i-1][j];
        }
    }

    /*printf("RateOfChangeArray: \n");
    for (int i = 0; i < (target -> depth) + 1; i++) {
        printArray((target -> rateOfChangeArray)[i], target -> sliceSize - i);
    }*/
}

void destructRateOfChangeArray(slice* target) {
    for (int j = 0; j < DEPTH + 1; j++) {
        free((target -> rateOfChangeArray)[j]);
    }
    free(target -> rateOfChangeArray);
}

slice* constructSlice(int sliceSize, int start, int end, int* baseDemandArray) {
    //printf("===== Slice Ctor =====\n");
    slice* result = (slice*) malloc(sizeof(slice));
   // printf("Passed result malloc\n");

    result -> sliceSize = sliceSize;
    result -> start = start;
    result -> end = end;
    result -> baseDemandArray = baseDemandArray;

    result -> rateOfChangeArray = (int**) malloc(sizeof(int*) * (DEPTH + 1));
    for (int j = 0; j < DEPTH + 1; j++) {
        (result -> rateOfChangeArray)[j] = (int*) malloc(sizeof(int) * (sliceSize - j));
    }
    //printf("Passed rateOfChangeArray malloc\n");
    constructRateOfChangeArray(result);
    //printf("Passed RateOfChangeArray ctor\n");

    return result;
}


void destructSlice(slice* target) {
    destructRateOfChangeArray(target);
    free(target); 
}

void setSliceSize(slice* target, int newSize) {
    destructRateOfChangeArray(target);
    target -> sliceSize = newSize;
    constructRateOfChangeArray(target);
}

void setStart(slice* target, int newStart) {
    destructRateOfChangeArray(target);
    target -> start = newStart;
    target -> end = newStart + (target -> end);
    constructRateOfChangeArray(target);
}

void printSlice(slice* target) {
    printf("Depth : %d\n", DEPTH);
    printf("Slice Size : %d\n", target -> sliceSize);
    printf("Start : %d, End : %d\n", target -> start, target -> end);
    for (int i = 0; i < DEPTH + 1; i++) {
        printf("Depth %d : \n\t", i);
        printArray((target -> rateOfChangeArray)[i], (target -> sliceSize) - i);
    }
}

void printVector(double* vector, int vectorSize) {
    printf("{ ");
    for (int i = 0; i < vectorSize; i++) {
        printf("%f", vector[i]);
        if (i != vectorSize - 1) printf(", ");
    }
    printf(" }");
}

double distance(int* vector1, int* vector2, int vectorSize) {
    double result;
    for (int i = 0 ; i < vectorSize; i++) {
        result += pow(vector2[i] - vector1[i], 2.0);
    }
    return sqrt(result);
}

double distance(double* vector1, double* vector2, int vectorSize) {
    double result;
    for (int i = 0 ; i < vectorSize; i++) {
        result += pow(vector2[i] - vector1[i], 2.0);
    }
    return sqrt(result);
}

typedef struct KDTreeNode2 {
    double vector[3];
    slice* origin;
    struct KDTreeNode2* left;
    struct KDTreeNode2* right;
} KDTreeNode2;

void printTreeNode(KDTreeNode2* node) {
    printf("Vector : \n\t");
    printArray(node -> vector, 3);
    printf("Origin slice : \n");
    printSlice(node -> origin);
}

// Create a new node

KDTreeNode2* createNode(double vector[], slice* origin) {
    KDTreeNode2* newNode = (KDTreeNode2*) malloc(sizeof(KDTreeNode2));
    for (int i = 0; i < 3; i++) (newNode -> vector)[i] = vector[i];
    printArray(newNode -> vector, 3);
    newNode -> origin = origin;
    newNode -> left = newNode -> right = NULL;
    return newNode;
}

// Insert a point into the k-d tree

KDTreeNode2* insert(KDTreeNode2* root, double vector[], slice* origin, int depth) {
    //printVector(vector, 3);

    if (root == NULL) {
        return createNode(vector, origin);
    }

    int cd = depth % 3;  // Current dimension (0=x, 1=y, 2=z)
    bool judgement = (cd == 0 && vector[0] < (root -> vector)[0]);
    for (int i = 1; i < 3; i++) 
        judgement = judgement || (cd == i && vector[i] < (root -> vector)[i]);
    
    // Compare the point to insert with the root point based on the current dimension
    if (judgement) {
       //printf("LEFT\n");
        root->left = insert(root->left, vector, origin, depth + 1);
    } else {
        //printf("RIGHT\n");
        root->right = insert(root->right, vector, origin, depth + 1);
    }
    //printf("!!!!!HERE!!!!!\n");

    return root;
}

// Function to find the closest point in the k-d tree to a given query point

KDTreeNode2* findClosestPoint(KDTreeNode2* root, double targetVector[], int depth, KDTreeNode2* best, double* bestDist) {
    if (root == NULL) {
        return best;
    }

    // Calculate the distance between the current node and the query point
    double dist = distance(root -> vector, targetVector, 3);
    if (dist < *bestDist) {
        *bestDist = dist;
        best = root;
    }

    // Calculate current dimension (0=x, 1=y, 2=z)
    int cd = depth % 3;
    bool judgement = (cd == 0 && targetVector[0] < (root -> vector)[0]);
    for (int i = 1; i < 3; i++) 
        judgement = judgement || (cd == i && targetVector[i] < (root -> vector)[i]);

    // Decide whether to go left or right based on the current dimension
    KDTreeNode2* nextBest = NULL;
    KDTreeNode2* otherSide = NULL;

    if (judgement) {
        nextBest = root->left;
        otherSide = root->right;
    } else {
        nextBest = root->right;
        otherSide = root->left;
    }

    // Recursively search in the selected subtree
    best = findClosestPoint(nextBest, targetVector, depth + 1, best, bestDist);

    // Check whether we need to search the other side of the tree
    // If the distance to the splitting plane is smaller than the best distance found so far
    double planeDist;
    for (int i = 0; i < 3; i++) {
        if (cd = i) {
            planeDist = fabs(targetVector[i] - (root -> vector)[i]);
            break;
        }
    }

    if (planeDist < *bestDist) {
        best = findClosestPoint(otherSide, targetVector, depth + 1, best, bestDist);
    }

    //printTreeNode(best);

    return best;
}

void freeKDTree(KDTreeNode2* root) {
    if (root == NULL) return;

    freeKDTree(root -> left);
    freeKDTree(root -> right);
    free(root);
}

void inorder(KDTreeNode2* root) { // LVR
    if (root == NULL) return;

    inorder(root -> left);
    printVector(root -> vector, 3);
    printf("\n");
    inorder(root -> right);
}

void preorder(KDTreeNode2* root) { // VLR
    if (root == NULL) return;

    printVector(root -> vector, 3);
    printf("\n");
    preorder(root -> left);
    preorder(root -> right);
}

// 통합 자료구조
typedef struct model {
    int* demandArray;
    int demandArraySize;
    int sliceSize;
    int sliceArraySize;
    slice** sliceArray;
    KDTreeNode2* KDTreeRoot;
} model;

model* constructModel(int* demandArray, int demandArraySize, int sliceSize) {
    model* result = (model*) malloc(sizeof(model));
    //printf("Passed result malloc\n");

    result -> demandArray = (int*) malloc(sizeof(int) * demandArraySize);
    result -> demandArraySize = demandArraySize;
    for (int i = 0; i < demandArraySize; i++) {
        result -> demandArray[i] = demandArray[i];
    }
    //printf("Passed demandArray malloc\n");

    result -> sliceSize = sliceSize;
    result -> sliceArraySize = demandArraySize - sliceSize + 1;
    result -> sliceArray = (slice**) malloc(sizeof(slice*) * (result -> sliceArraySize));
    /*for (int i = 0; i < result -> sliceArraySize; i++) {
        result -> sliceArray[i] = (slice*) malloc(sizeof(slice));
    }*/
    //printf("Passed sliceArray malloc\n");
    
    for (int i = 0; i < result -> sliceArraySize; i++) {
        //printf("Start Slice ctor idx : %d\n", i);
        (result -> sliceArray)[i] = constructSlice(sliceSize, i, i + sliceSize - 1, demandArray);
        //printSlice((result -> sliceArray)[i]);
        //printf("End Slice ctor idx : %d\n", i);
    }
    //printf("Passed sliceArray ctors\n");

    // 유사도 배열 생성
    // 비교 대상 : slice_array[size - sliceSize] (맨 끝 개체)

    int similarity_array_size = demandArraySize - sliceSize;
    double** similarity_array = (double**) malloc(sizeof(double*) * similarity_array_size);
    for (int i = 0; i < similarity_array_size; i++) {
        similarity_array[i] = (double*) malloc(sizeof(double) * 3);
    }
    //printf("Passed similarity_array malloc\n");

    for (int i = 0; i < similarity_array_size; i++) {
        for (int j = 0; j < 3; j++) {
            //printf("===== %d =====\n", i);
            //printSlice((result -> sliceArray)[i]);
            //printArray((result -> sliceArray)[i] -> rateOfChangeArray[j], sliceSize - j); 
            printArray((result -> sliceArray)[similarity_array_size] -> rateOfChangeArray[j], sliceSize - j); // <- 없으면 이상한 값이 붙는데.. 아마 메모리 할당 문제?
            similarity_array[i][j] = distance((result -> sliceArray)[i] -> rateOfChangeArray[j], (result -> sliceArray)[similarity_array_size] -> rateOfChangeArray[j], sliceSize - j);
        }
    }
    //printf("Similarity Array:\n");
    //printf("==============================\n");
    for (int i = 0; i < similarity_array_size; i++) {
        //printf("Idx %d : ", i);
        //printArray(similarity_array[i], 3);
    }

    //K-D 트리 KDTreeRoot
    result -> KDTreeRoot = NULL;
    for (int i = 0; i < similarity_array_size; i++) {
        //printArray(similarity_array[i], depth + 1);
        result -> KDTreeRoot = insert(result -> KDTreeRoot, similarity_array[i], result -> sliceArray[i], 0);
    }
    //printf("Inorder : \n");
    //inorder(result -> KDTreeRoot);
    //printf("Preorder : \n");
    //preorder(result -> KDTreeRoot);

    //printf("Construction Success\n");
    return result;
}

void destructModel(model* target) {
    //printf("===== destructModel =====\n");
    freeKDTree(target -> KDTreeRoot);

    //printf("Passed freeKDTree\n");
    for (int i = 0; i < target -> sliceArraySize; i++) {
        destructSlice(target -> sliceArray[i]);
    }
    //printf("Passed Slice dtors\n");

    free(target -> sliceArray);
    free(target -> demandArray);
    free(target);

    target = NULL;
}

void pushDemandNumber(model* target, int number) { // 맨 앞에 대상 추가
    int newDemandArraySize = (target -> demandArraySize) + 1;
    int sliceSize = target -> sliceSize;
    int* newDemandArray = (int*) malloc(sizeof(int) * newDemandArraySize);
    for (int i = 0 ; i < newDemandArraySize - 1; i++) {
        newDemandArray[i] = (target -> demandArray)[i];
    }
    newDemandArray[newDemandArraySize - 1] = number;
    //printf("Passed push before dtor\n");
    destructModel(target); // <- 트리 할당 해제 및 재할당 한 후 최근접 이웃 탐색 알고리즘에서 버그 발생하는 듯. 
    //printf("Passed push dtor\n");
    target = constructModel(newDemandArray, newDemandArraySize, sliceSize);
}

void deleteLastDemandNumber(model* target) { // 맨 끝 (0번 인덱스) 대상 삭제
    int newDemandArraySize = (target -> demandArraySize) - 1;
    int sliceSize = target -> sliceSize;
    int* newDemandArray = (int*) malloc(sizeof(int) * newDemandArraySize);
    for (int i = 1 ; i < newDemandArraySize + 1; i++) {
        newDemandArray[i] = (target -> demandArray)[i];
    }
    //printf("Passed delete before dtor\n");
    // destructModel(target);
   // printf("Passed delete dtor\n");
    target = constructModel(newDemandArray, newDemandArraySize, sliceSize);
}

void setSliceSize(model* target, int newSliceSize) {
    int tempDemandArraySize = (target -> demandArraySize);
    int* tempDemandArray = (int*) malloc(sizeof(int) * tempDemandArraySize);
    for (int i = 0; i < tempDemandArraySize; i++) {
        tempDemandArray[i] = (target -> demandArray)[i];
    }
    // destructModel(target);
    target = constructModel(tempDemandArray, tempDemandArraySize, newSliceSize);
}

int demandPrediction(model* target) {
    // 최단 거리 | 가장 가까운 이웃 탐색
    double targetPoint[3] = {0, 0, 0};

    KDTreeNode2* closest = NULL;
    double bestDistance = INFINITY;
    closest = findClosestPoint(target -> KDTreeRoot, targetPoint, 0, closest, &bestDistance);

    //printf("Closest Neighbor : \n");
    //printf("%d\n", closest);
    //printArray(closest -> vector, 3);


    // 최종 수요 예측

    slice* target_slice = closest -> origin;
    //printf("Target Slice : \n");
    //printSlice(target_slice);

    int result = 
        (target -> demandArray)[target -> demandArraySize - 1] + ((target -> demandArray)[(target_slice -> end) + 1] - (target -> demandArray)[target_slice -> end]);

    return result;
    
}

// 유클리드 거리가 0에 가까울수록 유사도 높은 것으로 상정.
// 각 depth에 대응하는 slice의 rateOfChangeArray의 인덱스에 접근, 

int main() {
    
    int size1 = 90;
    int demand1[90] = {
        123, 135, 141, 152, 168, 176, 188, 207, 199, 215,   // 1~10일 (주간 수요)
        242, 251, 237, 273, 263, 257, 291, 299, 308, 317,   // 11~20일 (중순 수요)
        332, 324, 314, 331, 359, 374, 386, 400, 417, 429,   // 21~30일 (월말 수요)
        97, 108, 113, 127, 134, 141, 158, 167, 154, 169,   // 1~10일 (주간 수요)
        172, 186, 178, 199, 211, 225, 242, 235, 248, 257,   // 11~20일 (중순 수요)
        271, 283, 295, 308, 320, 335, 348, 356, 367, 380,   // 21~30일 (월말 수요)
        105, 115, 120, 134, 145, 158, 165, 177, 189, 195,   // 1~10일 (주간 수요)
        201, 214, 220, 234, 246, 258, 272, 265, 280, 293,   // 11~20일 (중순 수요)
        305, 319, 328, 340, 355, 367, 380, 390, 401, 412,   // 21~30일 (월말 수요)
    };

    int size2 = 180;
    int demand2[180] = {
        123, 135, 141, 152, 168, 176, 188, 207, 199, 215,   // 1~10일 (주간 수요)
        242, 251, 237, 273, 263, 257, 291, 299, 308, 317,   // 11~20일 (중순 수요)
        332, 324, 314, 331, 359, 374, 386, 400, 417, 429,   // 21~30일 (월말 수요)
        97, 108, 113, 127, 134, 141, 158, 167, 154, 169,   // 1~10일 (주간 수요)
        172, 186, 178, 199, 211, 225, 242, 235, 248, 257,   // 11~20일 (중순 수요)
        271, 283, 295, 308, 320, 335, 348, 356, 367, 380,   // 21~30일 (월말 수요)
        105, 115, 120, 134, 145, 158, 165, 177, 189, 195,   // 1~10일 (주간 수요)
        201, 214, 220, 234, 246, 258, 272, 265, 280, 293,   // 11~20일 (중순 수요)
        305, 319, 328, 340, 355, 367, 380, 390, 401, 412,   // 21~30일 (월말 수요)
        98, 105, 112, 127, 135, 148, 160, 178, 169, 183,   // 1~10일 (주간 수요)
        211, 220, 206, 240, 229, 235, 265, 274, 285, 294,   // 11~20일 (중순 수요)
        310, 320, 309, 326, 340, 358, 373, 385, 398, 405,   // 21~30일 (월말 수요)
        89, 95, 101, 115, 124, 131, 148, 156, 145, 160,     // 31~40일 (주간 수요)
        167, 181, 173, 191, 204, 219, 233, 239, 252, 261,   // 41~50일 (중순 수요)
        274, 288, 302, 313, 325, 339, 354, 363, 375, 387,   // 51~60일 (월말 수요)
        92, 100, 108, 120, 134, 141, 153, 165, 179, 185,     // 61~70일 (주간 수요)
        193, 207, 212, 225, 237, 249, 261, 275, 286, 298,   // 71~80일 (중순 수요)
        312, 325, 336, 345, 359, 374, 385, 398, 410, 420    // 81~90일 (월말 수요)
    };

    int size3 = 360;
    int demand3[360] ={
         123, 135, 141, 152, 168, 176, 188, 207, 199, 215,   // 1~10일 (주간 수요)
        242, 251, 237, 273, 263, 257, 291, 299, 308, 317,   // 11~20일 (중순 수요)
        332, 324, 314, 331, 359, 374, 386, 400, 417, 429,   // 21~30일 (월말 수요)
        97, 108, 113, 127, 134, 141, 158, 167, 154, 169,   // 1~10일 (주간 수요)
        172, 186, 178, 199, 211, 225, 242, 235, 248, 257,   // 11~20일 (중순 수요)
        271, 283, 295, 308, 320, 335, 348, 356, 367, 380,   // 21~30일 (월말 수요)
        105, 115, 120, 134, 145, 158, 165, 177, 189, 195,   // 1~10일 (주간 수요)
        201, 214, 220, 234, 246, 258, 272, 265, 280, 293,   // 11~20일 (중순 수요)
        305, 319, 328, 340, 355, 367, 380, 390, 401, 412,   // 21~30일 (월말 수요)
        98, 105, 112, 127, 135, 148, 160, 178, 169, 183,   // 1~10일 (주간 수요)
        211, 220, 206, 240, 229, 235, 265, 274, 285, 294,   // 11~20일 (중순 수요)
        310, 320, 309, 326, 340, 358, 373, 385, 398, 405,   // 21~30일 (월말 수요)
        89, 95, 101, 115, 124, 131, 148, 156, 145, 160,     // 31~40일 (주간 수요)
        167, 181, 173, 191, 204, 219, 233, 239, 252, 261,   // 41~50일 (중순 수요)
        274, 288, 302, 313, 325, 339, 354, 363, 375, 387,   // 51~60일 (월말 수요)
        92, 100, 108, 120, 134, 141, 153, 165, 179, 185,     // 61~70일 (주간 수요)
        193, 207, 212, 225, 237, 249, 261, 275, 286, 298,   // 71~80일 (중순 수요)
        312, 325, 336, 345, 359, 374, 385, 398, 410, 420,    // 81~90일 (월말 수요)
       113, 127, 134, 145, 158, 170, 185, 198, 206, 215,   // 1~10일 (주간 수요)
        231, 240, 225, 265, 255, 268, 290, 298, 305, 318,   // 11~20일 (중순 수요)
        325, 340, 318, 330, 350, 370, 380, 390, 400, 415,   // 21~30일 (월말 수요)
        98, 112, 120, 133, 140, 150, 165, 175, 160, 175,     // 31~40일 (주간 수요)
        182, 193, 185, 198, 210, 223, 230, 240, 255, 265,   // 41~50일 (중순 수요)
        275, 290, 310, 315, 325, 335, 350, 360, 370, 385,   // 51~60일 (월말 수요)
        105, 115, 125, 137, 145, 155, 165, 175, 190, 200,    // 61~70일 (주간 수요)
        210, 225, 240, 250, 265, 275, 285, 295, 305, 315,   // 71~80일 (중순 수요)
        325, 335, 340, 350, 365, 375, 385, 395, 405, 415,   // 81~90일 (월말 수요)
        76, 88, 101, 110, 122, 135, 148, 161, 173, 181,   // 1~10일 (주간 수요)
        195, 208, 198, 225, 217, 228, 240, 253, 265, 278,   // 11~20일 (중순 수요)
        285, 295, 275, 305, 320, 335, 345, 355, 365, 375,   // 21~30일 (월말 수요)
        85, 98, 107, 116, 125, 132, 148, 159, 145, 157,     // 31~40일 (주간 수요)
        165, 175, 185, 197, 210, 222, 235, 240, 255, 262,   // 41~50일 (중순 수요)
        270, 285, 295, 305, 315, 325, 335, 345, 355, 365,   // 51~60일 (월말 수요)
        92, 103, 110, 121, 135, 142, 154, 165, 176, 185,    // 61~70일 (주간 수요)
        193, 207, 218, 230, 242, 250, 263, 278, 287, 295,   // 71~80일 (중순 수요)
        310, 320, 330, 340, 350, 360, 370, 380, 390, 400    // 81~90일 (월말 수요)
    };

    int bigDemand[5000];

    srand(time(NULL));

    for (int i = 0; i < 5000; i++) {
        bigDemand[i] = rand() % 100;  // 0부터 99 사이의 무작위 정수 생성
    }
    


    int sliceSize = 33;

    clock_t start, end;

    // bigDemand 

    // 자료구조 생성 -> 47.793s
    
    
    // start = clock();

    // model* bigDemandPredictionModel = constructModel(bigDemand, 5000, sliceSize);

    // end = clock();

    // double time_taken = (double)(end - start) / CLOCKS_PER_SEC;
    // printf("Time : %f\n", time_taken);
    
    

    // 데이터 삽입 / 삭제 -> 71.159s

    
    // model* bigDemandPredictionModel = constructModel(bigDemand, 5000, sliceSize);
    
    // start = clock();

    // pushDemandNumber(bigDemandPredictionModel, 99);

    // end = clock();

    // double time_taken = (double)(end - start) / CLOCKS_PER_SEC;
    // printf("Time : %f\n", time_taken);
    
    

    // 수요 예측 결과 산출 -> 0s
    
    model* bigDemandPredictionModel = constructModel(bigDemand, 5000, sliceSize);
    
    start = clock();

    int result = demandPrediction(bigDemandPredictionModel);
    printf("Result: %d\n", result);

    end = clock();

    double time_taken = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time : %f\n", end - start);
    

    // demand1 (size1 = 90)

    // 자료구조 생성 -> 0.333000s
    
    /*
    start = clock();

    model* demand1PredictionModel = constructModel(demand1, size1, sliceSize);

    end = clock();

    double time_taken = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time : %f\n", time_taken);
    */
    

    // 데이터 삽입 / 삭제 -> 0.315000s

    /*
    model* demand1PredictionModel = constructModel(demand1, size1, sliceSize);
    
    start = clock();

    pushDemandNumber(demand1PredictionModel, 99);

    end = clock();

    double time_taken = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time : %f\n", time_taken);
    */
    

    // 수요 예측 결과 산출 -> 0.000000s
    /*
    model* demand1PredictionModel = constructModel(demand1, size1, sliceSize);
    
    start = clock();

    int result = demandPrediction(demand1PredictionModel);
    printf("Result: %d\n", result);

    end = clock();

    double time_taken = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time : %f\n", end - start);
    */



    // demand2 (size2 = 180)

    // 자료구조 생성 -> 0.687000s

    
    // start = clock();    

    // model* demand2PredictionModel = constructModel(demand2, size2, sliceSize);

    // end = clock();

    // double time_taken = (double)(end - start) / CLOCKS_PER_SEC;
    // printf("Time : %f\n", time_taken);
    

    // 데이터 삽입 / 삭제 -> 0.585000s

    // model* demand2PredictionModel = constructModel(demand2, size2, sliceSize);

    // start = clock();

    // pushDemandNumber(demand2PredictionModel, 99);

    // end = clock();

    // double time_taken = (double)(end - start) / CLOCKS_PER_SEC;
    // printf("Time : %f\n", time_taken);

    // 수요 예측 결과 산출 -> 0.000000s

    // model* demand2PredictionModel = constructModel(demand2, size2, sliceSize);

    // start = clock(); 

    // int result = demandPrediction(demand2PredictionModel);
    // printf("Result: %d\n", result);

    // end = clock();

    // double time_taken = (double)(end - start) / CLOCKS_PER_SEC;
    // printf("Time : %f\n", time_taken);



    // demand3 (size1 = 360)

    // 자료구조 생성 -> 1.528000s

    // start = clock();

    // model* demand3PredictionModel = constructModel(demand3, size3, sliceSize);

    // end = clock();

    // double time_taken = (double)(end - start) / CLOCKS_PER_SEC;
    // printf("Time : %f\n", time_taken);

    // 데이터 삽입 / 삭제 -> 1.248000s

    // model* demand3PredictionModel = constructModel(demand3, size3, sliceSize);

    // start = clock();

    // pushDemandNumber(demand3PredictionModel, 99);

    // end = clock();

    // double time_taken = (double)(end - start) / CLOCKS_PER_SEC;
    // printf("Time : %f\n", time_taken);

    // 수요 예측 결과 산출 -> 0.0000000s

    // model* demand3PredictionModel = constructModel(demand3, size3, sliceSize);
    
    // start = clock();

    // int result = demandPrediction(demand3PredictionModel);
    // printf("Result: %d\n", result);

    // end = clock();

    // double time_taken = (double)(end - start) / CLOCKS_PER_SEC;
    // printf("Time : %f\n", end - start);

    /*int example[14] = {
        17, 30, 21, 37, 44, 59, 61, 21, 16, 20, 21, 26, 30, 31
    };

    int ex_sliceSize = 6;

    model* examplePredictionModel = constructModel(example, 14, ex_sliceSize);

    //printArray(examplePredictionModel -> KDTreeRoot -> vector, 3); 

    int example_result = demandPrediction(examplePredictionModel);
    printf("Prediction Result: %d\n", example_result);*/

    
    return 0;
}


