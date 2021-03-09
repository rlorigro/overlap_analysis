#ifndef OVERLAP_ANALYSIS_SPAN_HPP
#define OVERLAP_ANALYSIS_SPAN_HPP

#include <stdexcept>
#include <string>

// Gcc (for backtraces).
#include "execinfo.h"

void writeBackTrace();

#define ASSERT(expression) ((expression) ? (static_cast<void>(0)) : \
    (/*writeBackTrace(),*/ throw std::runtime_error(std::string("Assertion failed: ") + #expression + " at " + __PRETTY_FUNCTION__ + " in " +  __FILE__ + " line " + std::to_string(__LINE__))))


// A span class similar to std::span in C++20.
// We currently compile using the C++14 standard, so we cannot use std::span.

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

using std::ostream_iterator;
using std::ostream;
using std::runtime_error;


template<class T> class span;

// Output a span<const char> as a string.
inline ostream& operator<<(ostream&, const span<const char>&);

// Convert a span<const char> to an integer.
// This cannot be done using std::atol because the
// span<const char> is not null terminated.
uint64_t atoul(const span<const char>&);

// Convert a span<const char> to an std::string.
string convertToString(const span<const char>&);


template<class T> class span {
public:

    span(T* begin, T* end) :
            dataBegin(begin),
            dataEnd(end)
    {
    }

    span(vector<T>& v) :
            dataBegin(&v[0]),
            dataEnd(dataBegin + v.size())
    {
    }

    span() : dataBegin(0), dataEnd(0) {}

    size_t size() const
    {
        return dataEnd - dataBegin;
    }
    bool empty() const
    {
        return dataBegin == dataEnd;
    }
    T* begin() const
    {
        return dataBegin;
    }
    T* end() const
    {
        return dataEnd;
    }
    T& operator[](size_t i) const
    {
        return dataBegin[i];
    }

    T& front()
    {
        ASSERT(dataBegin);
        ASSERT(dataEnd);
        ASSERT(!empty());
        return *dataBegin;
    }
    const T& front() const
    {
        ASSERT(dataBegin);
        ASSERT(dataEnd);
        ASSERT(!empty());
        return *dataBegin;
    }

    T& back()
    {
        ASSERT(dataBegin);
        ASSERT(dataEnd);
        ASSERT(!empty());
        return *(dataEnd - 1);
    }

    const T& back() const
    {
        ASSERT(dataBegin);
        ASSERT(dataEnd);
        ASSERT(!empty());
        return *(dataEnd - 1);
    }

    bool operator==(const span<T>& that) const
    {
        if(size() == that.size()) {
            return std::equal(begin(), end(), that.begin());
        } else {
            return false;
        }
    }
    bool operator!=(const span<T>& that) const
    {
        return not((*this) == that);
    }
private:
    T* dataBegin;
    T* dataEnd;
};



// Output a span<const char> as a string.
inline std::ostream& operator<<(
        std::ostream& s,
        const span<const char>&  m)
{
    copy(m.begin(), m.end(), ostream_iterator<char>(s));
    return s;
}


// Convert a span<const char> to an integer.
// This cannot be done using std::atol because the
// span<const char> is not null terminated.
// The string can only contain numeric characters.
inline uint64_t atoul(const span<const char>& s)
{
    uint64_t n = 0;
    for(uint64_t i=0; ; i++) {
        const char c = s[i];
        if(not std::isdigit(c)) {
            throw runtime_error("Non-digit found in " + convertToString(s));
        }

        n += (c- '0');

        if(i ==  s.size()-1) {
            return n;
        }

        n *= 10;
    }
}



// Convert a span<const char> to an std::string.
inline std::string convertToString(const span<const char>& m)
{
    return string(m.begin(), m.size());
}



#endif //OVERLAP_ANALYSIS_SPAN_HPP
