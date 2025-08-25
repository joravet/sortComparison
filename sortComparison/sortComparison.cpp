//Jackson Oravetz; Ahmed, Kishwar. Code to run three different sorting algorithms: insetion, merge, and radix

#include <iostream>
#include <vector>
#include <chrono>
#include <ctime>
#include <ratio>
#include <math.h>
using namespace std;

void insertionSort(std::vector<int>& numbers) {
    for (int i = 1; i < numbers.size(); i++) {
        int key = numbers.at(i);
        int j = i - 1;  //inset numbers at i into the sorted subvector numbers 0:i-1
            while (j >= 0 && numbers.at(j) > key) {
                numbers.at(j + 1) = numbers.at(j);
                j = j - 1;
            }
            numbers.at(j + 1) = key;
    }

}

void merge(std::vector<int>& numbers, int p, int q, int r) {
    int n_left = q - p + 1; //length of numbers p:q
    int n_right = r - q;    //length of numbers q+1:r
    vector<int> leftVector (n_left);
    vector<int> rightVector (n_right);
    for (int i = 0; i <= n_left - 1; i++)   //copy numbers p:q into leftVector 0:n_left-1
        leftVector.at(i) = numbers.at(p + i);
    for (int j = 0; j <= n_right - 1; j++)  //copy numbers q+1:r into rightVector 0:n_right-1
        rightVector.at(j) = numbers.at(q + j + 1);
    int i = 0;  //i indexes the smallest remaining element in leftVector
    int j = 0;  //j indexes the smallest remaining element in rightVector
    int k = p;  //k indexes the location in numbers to fill

    while (i < n_left && j < n_right) { //As long leftVector and rightVector contain unmerged elements
                                        //copy the smallest unmerged element back into elements p:r
        if (leftVector.at(i) <= rightVector.at(j)) {
            numbers.at(k) = leftVector.at(i);
            i = i + 1;
        }
        else {
            numbers.at(k) = rightVector.at(j);
            j = j + 1;
        }
        k = k + 1;
    }

    while (i < n_left) {    //Having gone through one of leftVector and rightVector entirely
                            //copy the remainder of the other to the end of numbers p:r
        numbers.at(k) = leftVector.at(i);
        i = i + 1;
        k = k + 1;
    }

    while (j < n_right) {
        numbers.at(k) = rightVector.at(j);
        j = j + 1;
        j = k + 1;
    }
}

void mergeSort(std::vector<int>& numbers, int p, int r) {
    if (p >= r)     //zero or one element
        return;
    int q = floor((p + r) / 2); //midpoint of numbers from p:r
    mergeSort(numbers, p, q);   //recursively sort numbers from p:q
    mergeSort(numbers, q + 1, r);   //recursively sort numbers from q+1:r
    merge(numbers, p, q, r);    //merge numbers p:q and numbers q+1:r into numbers p:r
}

int radixGetLength(int value) { //simple length check
    if (value == 0) 
        return 1;
    int digits = 0;
    while (value != 0) {
        digits = digits + 1;
        value = value / 10;
    }
    return digits;
}

int radixGetMaxLength(std::vector<int>& numbers) {
    int maxDigits = 0;
    for (int i = 0; i < numbers.size(); i++) {  //loop for entire vector
        int digitCount = radixGetLength(numbers.at(i)); //get vector length
        if (digitCount > maxDigits) //if bigger than maxDigits, make new max length
            maxDigits = digitCount;
    }

    return maxDigits;
}

void radixSort(std::vector<int>& numbers) {
    vector<vector<int>> buckets(10);    //vector of buckets
    int maxDigits = radixGetMaxLength(numbers); //find the mx length, in number of digits
    int pow10 = 1;  //start with least significant digit

    for (int digitIndex = 0; digitIndex < maxDigits; digitIndex++) {    //loop for each digit
        for (int i = 0; i < numbers.size(); i++) {  //loop for each index of vector
            int bucketIndex = abs(numbers.at(i) / pow10) % 10;
            buckets.at(bucketIndex).push_back(numbers.at(i));   //append numbers at i to bucketIndex
        }
        int vectorIndex = 0;
        for (int i = 0; i < 10; i++) {  //nested loop
            for (int j = 0; j < buckets.at(i).size(); j++)  
               numbers.at(vectorIndex++) = buckets.at(i).at(j); //access bucket i, index j
        }
        pow10 = 10 * pow10;
        for (int i = 0; i < 10; i++) //clear all buckets
            buckets.at(i).clear();
    }
}


int main()
{
    vector<int> test;   //test vector
    int vecSize = 10;   //size of vector
    int testNum = 0;    //trial run
    double time = 0;    //timevar
    //insert portion
    while (testNum != 4) {  //loop 4 times
        time = 0;   //reset time variable
        for (int i = 0; i < 10; i++) {  //loop for ten trials
            for (int j = 0; j < vecSize; j++) { //fill vector to desired size with int from 0-999
                test.push_back(rand() % 1000);  
            }
            chrono::steady_clock::time_point t1 = chrono::steady_clock::now();  //start clock
            insertionSort(test);    //run sort
            chrono::steady_clock::time_point t2 = chrono::steady_clock::now();  //stop clock
            chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
            time += time_span.count();  //track total time for all 10 trials
            test.clear();   //clear vector for next trial
        }
        time /= 10.0;   //divide by 10 to find avg time
        cout << "The insertion sort wall clock time at N = " <<  vecSize << " is: " << time << endl;    //ouput avg wall clock time
        vecSize *= 10;  //increase vecSize for next trial
        testNum++;  //next trial
    }

    //code works identically to insert portion, just with merge sort
    cout << endl;
    vecSize = 10;
    testNum = 0;
    time = 0;
    while (testNum != 4) {
        time = 0;
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < vecSize; j++) {
                test.push_back(rand() % 1000);
            }
            chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
            mergeSort(test, 0, vecSize - 1);
            chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
            chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
            time += time_span.count();
            test.clear();
        }
        time /= 10.0;
        cout << "The merge sort wall clock time at N = " << vecSize << " is: " << time << endl;
        vecSize *= 10;
        testNum++;
    }

    //code works identically to insetion portion, just with radix sort
    cout << endl;
    vecSize = 10;
    testNum = 0;
    time = 0;
    while (testNum != 4) {
        time = 0;
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < vecSize; j++) {
                test.push_back(rand() % 1000);
            }
            chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
            radixSort(test);
            chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
            chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
            time += time_span.count();
            test.clear();
        }
        time /= 10.0;
        cout << "The radix sort wall clock time at N = " << vecSize << " is: " << time << endl;
        vecSize *= 10;
        testNum++;
    }
 }