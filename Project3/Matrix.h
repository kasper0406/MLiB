#pragma once

#include <vector>
#include <functional>
#include <cassert>

using namespace std;

template<class T>
class Matrix
{
public:
    Matrix(size_t n, size_t m, const T& value) : n(n), m(m), elements(vector<T>(n * m, value))
    { }
    
    inline T operator()(size_t row, size_t column) const {
        assert(row < n && column < m);
        return elements[m * row + column];
    }
    
    inline T& operator()(size_t row, size_t column) {
        assert(row < n && column < m);
        return elements[m * row + column];
    }
    
    void map(function<T(T)> op) {
        for (size_t i = 0; i < n * m; i++)
            elements[i] = op(elements[i]);
    }
    
private:
    const size_t n, m;
    vector<T> elements;
};
